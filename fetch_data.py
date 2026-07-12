#!/usr/bin/env python3
"""Build the person-temperature-hours dataset for the humantemp site.

Combines three freely-accessible, no-authentication sources onto one regular
lat/lon grid, then writes compact JSON the web page renders client-side:

  * Temperature  ERA5 hourly 2m air temperature, streamed from the anonymous
                 ARCO-ERA5 Zarr store on Google Cloud (0.25deg, coarsened here).
  * Population   GHS-POP 2020, WGS84 30 arcsec (=1/120 deg), block-summed to
                 the grid (population counts add exactly, no resampling).
  * Continents   Natural Earth 110m country polygons dissolved by continent,
                 rasterized onto the grid (this also serves as the land mask).

For every land (or populated) cell we accumulate, over the whole year, the
number of hours spent in each 1 deg-C temperature bin. The page multiplies
those hour counts by population (person-hours) or cell area, and re-aggregates
under whatever region filter is active.

Nothing raw is kept on disk: each block of timesteps is downloaded, coarsened,
histogrammed, and discarded.

Usage:
  pip install -r requirements.txt
  python fetch_data.py                 # 1deg grid, full year 2024, hourly
  python fetch_data.py --days 7        # quick smoke test (first week only)
  python fetch_data.py --deg 0.5 --stride 3
"""

import argparse
import io
import json
import os
import time
import zipfile

import numpy as np
import requests

ARCO = "gs://gcp-public-data-arco-era5/ar/full_37-1h-0p25deg-chunk-1.zarr-v3"
GHSPOP_URL = (
    "https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_POP_GLOBE_R2023A/"
    "GHS_POP_E2020_GLOBE_R2023A_4326_30ss/V1-0/"
    "GHS_POP_E2020_GLOBE_R2023A_4326_30ss_V1_0.zip"
)
CONTINENTS_URL = (
    "https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/"
    "geojson/ne_110m_admin_0_countries.geojson"
)
OCEAN_LABEL = "Seven seas (open ocean)"

# Temperature binning (degrees Celsius). Covers Earth's inhabited-plus range.
TMIN, TMAX, BINW = -70.0, 56.0, 1.0


def log(msg):
    print(msg, flush=True)


# --------------------------------------------------------------------------- #
# Population: GHS-POP 30 arcsec -> block-sum onto the target grid.
# --------------------------------------------------------------------------- #
def population_grid(deg, cache_dir):
    import rasterio

    os.makedirs(cache_dir, exist_ok=True)
    tif = os.path.join(cache_dir, "ghs_pop_2020_30ss.tif")
    if not os.path.exists(tif):
        zpath = os.path.join(cache_dir, "ghs_pop.zip")
        if not os.path.exists(zpath):
            log("Downloading GHS-POP 2020 30 arcsec (~460 MB, one-time)...")
            with requests.get(GHSPOP_URL, stream=True, timeout=120) as r:
                r.raise_for_status()
                with open(zpath, "wb") as f:
                    for chunk in r.iter_content(1 << 20):
                        f.write(chunk)
        log("Extracting GHS-POP GeoTIFF...")
        with zipfile.ZipFile(zpath) as z:
            name = next(n for n in z.namelist() if n.lower().endswith(".tif"))
            with z.open(name) as src, open(tif, "wb") as dst:
                dst.write(src.read())

    from rasterio.windows import Window

    nlat, nlon = int(round(180 / deg)), int(round(360 / deg))
    grid = np.zeros(nlat * nlon, dtype=np.float64)
    with rasterio.open(tif) as ds:
        H, W = ds.height, ds.width
        t = ds.transform  # a=pixel width, e=pixel height (neg), c/f=origin x/y
        log(f"GHS-POP raster {W}x{H}, origin=({t.c:.3f},{t.f:.3f}), "
            f"pixel={t.a:.6f} deg")
        # The raster is offset from a clean grid, so bin each source pixel into
        # the target cell that contains its centre (exact for any origin/extent).
        lon = t.c + (np.arange(W) + 0.5) * t.a
        tj = np.mod(np.floor((lon + 180.0) / deg).astype(np.int64), nlon)
        for r0 in range(0, H, 2048):
            n = min(2048, H - r0)
            band = ds.read(1, window=Window(0, r0, W, n)).astype(np.float64)
            band[band < 0] = 0.0  # nodata (-200)
            lat = t.f + (np.arange(r0, r0 + n) + 0.5) * t.e
            ti = np.clip(np.floor((90.0 - lat) / deg).astype(np.int64), 0, nlat - 1)
            lin = (ti[:, None] * nlon + tj[None, :]).ravel()
            grid += np.bincount(lin, weights=band.ravel(), minlength=nlat * nlon)
    grid = grid.reshape(nlat, nlon)
    log(f"Population grid {grid.shape}, total = {grid.sum() / 1e9:.2f} billion")
    return grid  # (nlat, nlon), persons per cell


# --------------------------------------------------------------------------- #
# Continents + land mask: rasterize NE polygons onto the grid.
# --------------------------------------------------------------------------- #
def continent_grid(lat_c, lon_c, cache_dir):
    import geopandas as gpd
    import regionmask

    gj = os.path.join(cache_dir, "ne_110m_countries.geojson")
    if not os.path.exists(gj):
        log("Downloading Natural Earth continents...")
        r = requests.get(CONTINENTS_URL, timeout=60)
        r.raise_for_status()
        with open(gj, "wb") as fh:
            fh.write(r.content)
    world = gpd.read_file(gj)
    col = next(c for c in world.columns if c.upper() == "CONTINENT")
    world = world[world[col] != OCEAN_LABEL]
    cont = world.dissolve(by=col).reset_index()
    regions = regionmask.from_geopandas(cont, names=col, name="continent")
    names = list(regions.names)
    idx = regions.mask(lon_c, lat_c).values  # (nlat,nlon) float codes / NaN
    log(f"Continents: {names}")
    return idx, names


def fill_missing_continents(cont_idx, include):
    """Assign a continent to included cells that fell on no polygon (coastline)
    by copying the nearest cell that has one."""
    have = np.isfinite(cont_idx)
    need = include & ~have
    if not need.any():
        return cont_idx
    hi, hj = np.where(have)
    ni, nj = np.where(need)
    for i, j in zip(ni, nj):
        d = (hi - i) ** 2 + (hj - j) ** 2
        cont_idx[i, j] = cont_idx[hi[d.argmin()], hj[d.argmin()]]
    return cont_idx


# --------------------------------------------------------------------------- #
# Temperature: stream ERA5 hourly, coarsen, accumulate per-cell histograms.
# --------------------------------------------------------------------------- #
def temperature_hist(deg, year, stride, days, nbins):
    import xarray as xr

    ds = xr.open_zarr(ARCO, storage_options={"token": "anon"}, chunks=None)
    da = ds["2m_temperature"]
    sub = da.sel(time=slice(f"{year}-01-01T00:00", f"{year}-12-31T23:00"))
    if stride > 1:
        sub = sub.isel(time=slice(None, None, stride))
    nT = sub.sizes["time"]
    if days:
        nT = min(nT, days * 24 // stride)

    f = int(round(deg / 0.25))          # native cells per output cell
    nlat, nlon = 720 // f, 1440 // f
    ncell = nlat * nlon
    roll = 1440 // 2                     # shift lon 0..360 -> -180..180

    acc = np.zeros((ncell, nbins), dtype=np.int32)   # hours per (cell,bin)
    tsum = np.zeros(ncell, dtype=np.float64)
    hours_per_sample = stride
    B = 48
    t0 = time.time()
    log(f"Streaming {nT} timesteps of ERA5 (stride={stride}h) -> {nlat}x{nlon} grid")
    for start in range(0, nT, B):
        n = min(B, nT - start)
        block = sub.isel(time=slice(start, start + n)).values  # (n,721,1440) K
        block = np.roll(block, roll, axis=2)[:, :720, :]
        coarse = block.reshape(n, nlat, f, nlon, f).mean(axis=(2, 4))  # (n,nlat,nlon) K
        coarse = coarse.reshape(n, ncell) - 273.15                     # -> Celsius
        tsum += coarse.sum(axis=0)
        bins = np.clip(((coarse - TMIN) / BINW).astype(np.int64), 0, nbins - 1)
        cell = np.arange(ncell, dtype=np.int64)[None, :]
        lin = (cell * nbins + bins).ravel()
        acc += np.bincount(lin, minlength=ncell * nbins).reshape(ncell, nbins).astype(np.int32)
        done = start + n
        if (start // B) % 10 == 0 or done == nT:
            el = time.time() - t0
            log(f"  {done}/{nT} timesteps  ({el:.0f}s, eta {el/done*(nT-done):.0f}s)")

    acc *= hours_per_sample             # samples -> hours
    tmean = tsum / nT                   # mean Celsius per cell
    return acc.reshape(nlat, nlon, nbins), tmean.reshape(nlat, nlon), nT * hours_per_sample


# --------------------------------------------------------------------------- #
# Assemble + write JSON.
# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--deg", type=float, default=1.0, help="grid resolution (deg)")
    ap.add_argument("--year", type=int, default=2024)
    ap.add_argument("--stride", type=int, default=1, help="hours between samples")
    ap.add_argument("--days", type=int, default=None, help="limit to first N days (testing)")
    ap.add_argument("--out", default="data")
    ap.add_argument("--cache", default=".cache")
    args = ap.parse_args()

    deg = args.deg
    nbins = int(round((TMAX - TMIN) / BINW))
    nlat, nlon = int(round(180 / deg)), int(round(360 / deg))
    lat_c = 90 - (np.arange(nlat) + 0.5) * deg   # 89.5 .. -89.5
    lon_c = -180 + (np.arange(nlon) + 0.5) * deg  # -179.5 .. 179.5

    pop = population_grid(deg, args.cache)
    cont_idx, cont_names = continent_grid(lat_c, lon_c, args.cache)
    hist, tmean, hours_total = temperature_hist(deg, args.year, args.stride, args.days, nbins)

    land = np.isfinite(cont_idx)
    include = land | (pop > 0)
    cont_idx = fill_missing_continents(cont_idx.copy(), include)

    # Emit cells columnar + run-length histograms (first non-zero bin + counts).
    lat_i, lon_i = np.where(include)
    lats, lons, pops, conts, tmeans = [], [], [], [], []
    t0s, hists = [], []
    for i, j in zip(lat_i.tolist(), lon_i.tolist()):
        h = hist[i, j]
        nz = np.nonzero(h)[0]
        if nz.size == 0:
            continue
        a, b = int(nz[0]), int(nz[-1]) + 1
        lats.append(round(float(lat_c[i]), 3))
        lons.append(round(float(lon_c[j]), 3))
        pops.append(int(round(float(pop[i, j]))))
        c = cont_idx[i, j]
        conts.append(int(c) if np.isfinite(c) else -1)
        tmeans.append(round(float(tmean[i, j]), 1))
        t0s.append(a)
        hists.append(h[a:b].astype(int).tolist())

    os.makedirs(args.out, exist_ok=True)
    meta = {
        "year": args.year,
        "deg": deg,
        "stride_hours": args.stride,
        "hours_total": int(hours_total),
        "tmin": TMIN, "tmax": TMAX, "binw": BINW, "nbins": nbins,
        "continents": cont_names,
        "ncells": len(lats),
        "source": "ERA5 (ARCO-ERA5, ECMWF/Copernicus) + GHS-POP 2020 (EU JRC) "
                  "+ Natural Earth",
    }
    with open(os.path.join(args.out, "meta.json"), "w") as f:
        json.dump(meta, f)
    cells = {
        "lat": lats, "lon": lons, "pop": pops, "cont": conts,
        "tmean": tmeans, "t0": t0s, "hist": hists,
    }
    with open(os.path.join(args.out, "cells.json"), "w") as f:
        json.dump(cells, f, separators=(",", ":"))

    size = os.path.getsize(os.path.join(args.out, "cells.json"))
    log(f"Wrote {len(lats)} cells, cells.json = {size/1e6:.1f} MB")
    log(f"Population represented: {sum(pops)/1e9:.2f} billion")


if __name__ == "__main__":
    main()
