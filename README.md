# The temperatures humans live in

A static, client-side visualization of **person-hours by temperature**: for every hour of a
year, how many people worldwide were experiencing each air temperature. Combines hourly ERA5
temperature with gridded population on a common lat/lon grid.

- **Map** of the globe with selectable views: mean temperature, **mean daily high**, **mean daily
  low**, **annual maximum**, **annual minimum**, **mean daily range**, or population. Drag to
  select a box; hover a chart bar to see where that temperature is felt.
- **Distribution** of person-hours across temperature, with a slider that draws symmetric
  **percentile lines** (e.g. 95 → the 5th and 95th percentile temperatures of the selection).
- **Daily-swing heatmap**: a 2D histogram of each place-day's overnight low vs daytime high,
  shaded (white→black) by people-days — points far above the diagonal are big daily swings.
- Filters applied to all panels: weight by **population** (person-hours) or **land area**
  (km²-hours), and restrict to **continents** or a **lat/lon box**.

Both **°C and °F** are shown at once — dual axes on the charts, and paired labels everywhere
else. Everything the page needs is precomputed into `data/*.json`; all aggregation runs in the
browser. Colour is reserved for absolute temperature — a diverging scale **centered on 0 °C**
(blue below freezing, red above) — while every other magnitude (population, daily range,
people-days) is grayscale so it never reads as a temperature.

## Data sources (all free, no account required)

| Layer | Source | Access |
|------|--------|--------|
| Temperature | ERA5 hourly 2 m air temperature (ECMWF / Copernicus) | [ARCO-ERA5](https://cloud.google.com/storage/docs/public-datasets/era5) Zarr on Google Cloud, anonymous |
| Population | [GHS-POP 2020](https://human-settlement.emergency.copernicus.eu/) 30 arcsec (EU JRC) | direct HTTP download |
| Land & continents | [Natural Earth](https://www.naturalearthdata.com/) 110m | GeoJSON |

## Regenerate the data

```bash
pip install -r requirements.txt
python fetch_data.py                 # 1° grid, full year 2024, hourly (~7 min, streams ~36 GB)
```

The script streams ERA5 one block of timesteps at a time, coarsens each to the grid, and
accumulates per-cell temperature histograms — nothing raw is stored. Options:

```bash
python fetch_data.py --days 7        # quick smoke test (first week only)
python fetch_data.py --deg 0.5       # finer grid (larger data file, longer run)
python fetch_data.py --stride 3      # sample every 3rd hour (smaller download)
python fetch_data.py --year 2023
```

Outputs three files: `data/meta.json` (grid + bin definitions), `data/cells.json` (per-cell
population, continent, mean temperature, mean daily range, and run-length
hour-per-temperature-bin histogram), and `data/daily.json` (per-cell sparse 2D histogram of
days by (daily-min, daily-max) 2 °C bins, aligned to the same cell order).

## Run locally

```bash
python -m http.server 8000
# open http://localhost:8000
```

## Publish on GitHub Pages

Commit everything (including `data/`) and enable Pages for the branch, serving from the repo
root. No build step. `.nojekyll` is included so all files are served as-is.

## Notes on method

- Grid is 1°; ERA5 (0.25°) is averaged to it, GHS-POP (30 arcsec) is summed to it.
- A cell is included if it is land or has population. "Land area" weighting uses true cell area
  (∝ cos lat); "population" weighting uses GHS-POP counts.
- Temperature bins are 1 °C wide. Each cell stores hours spent in each bin over the year, so
  `person-hours(bin) = Σ population × hours` and `area-hours(bin) = Σ area × hours` are
  recomputed for any region filter in the browser.
- Daily min/max are computed per calendar day from the same hourly stream, then binned into a
  2 °C × 2 °C grid per cell (people-days = Σ population × days). Mean daily range is the
  per-cell mean of (daily max − daily min). Temperature *differences* convert to °F as ×9/5
  (no +32 offset); absolute temperatures as ×9/5 + 32.
