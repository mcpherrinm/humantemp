"use strict";

// ---------------------------------------------------------------------------
// Color scales. Temperature uses one diverging blue->red ramp, shared by the
// map and the chart bars, so a color means the same thing everywhere. Anything
// that is not temperature stays grayscale.
// ---------------------------------------------------------------------------
const TEMP_STOPS = [
  [49, 54, 149], [69, 117, 180], [116, 173, 209], [171, 217, 233],
  [224, 243, 248], [255, 255, 191], [254, 224, 144], [253, 174, 97],
  [244, 109, 67], [215, 48, 39], [165, 0, 38],
];
const TEMP_DOMAIN = [-30, 45];

function lerp(a, b, t) { return a + (b - a) * t; }

function tempColor(tc) {
  let f = (tc - TEMP_DOMAIN[0]) / (TEMP_DOMAIN[1] - TEMP_DOMAIN[0]);
  f = Math.max(0, Math.min(1, f));
  const x = f * (TEMP_STOPS.length - 1);
  const i = Math.min(TEMP_STOPS.length - 2, Math.floor(x));
  const u = x - i, a = TEMP_STOPS[i], b = TEMP_STOPS[i + 1];
  return [lerp(a[0], b[0], u), lerp(a[1], b[1], u), lerp(a[2], b[2], u)];
}
const rgb = (c, alpha) => `rgba(${c[0] | 0},${c[1] | 0},${c[2] | 0},${alpha == null ? 1 : alpha})`;

function popColor(frac) { // grayscale, light -> dark
  const g = Math.round(lerp(238, 24, Math.max(0, Math.min(1, frac))));
  return `rgb(${g},${g},${g})`;
}

// ---------------------------------------------------------------------------
// State
// ---------------------------------------------------------------------------
let META, C, N;                 // metadata, cells (columnar), cell count
let AREA;                       // per-cell area weight (km^2)
let state = {
  weight: "pop",                // 'pop' | 'area'
  mapview: "temp",              // 'temp' | 'pop'
  continents: new Set(),        // active continent indices
  box: null,                    // {latMin,latMax,lonMin,lonMax}
  hoverBin: null,
};

const $ = (id) => document.getElementById(id);

// ---------------------------------------------------------------------------
// Load
// ---------------------------------------------------------------------------
async function load() {
  [META, C] = await Promise.all([
    fetch("data/meta.json").then((r) => r.json()),
    fetch("data/cells.json").then((r) => r.json()),
  ]);
  N = C.lat.length;
  const R2 = 6371 * 6371, d = (META.deg * Math.PI) / 180;
  AREA = new Float64Array(N);
  for (let k = 0; k < N; k++) {
    AREA[k] = R2 * d * d * Math.cos((C.lat[k] * Math.PI) / 180);
  }
  META.continents.forEach((_, i) => state.continents.add(i));
  $("src").textContent =
    `${META.ncells.toLocaleString()} cells at ${META.deg}°, ` +
    `${META.hours_total.toLocaleString()} hours of ${META.year}. Source: ${META.source}.`;
  buildControls();
  window.addEventListener("resize", () => { drawMap(); drawChart(); });
  update();
}

// ---------------------------------------------------------------------------
// Filtering + aggregation
// ---------------------------------------------------------------------------
function passes(k) {
  if (!state.continents.has(C.cont[k])) return false;
  const b = state.box;
  if (b) {
    const la = C.lat[k], lo = C.lon[k];
    if (la < b.latMin || la > b.latMax || lo < b.lonMin || lo > b.lonMax) return false;
  }
  return true;
}
const wt = (k) => (state.weight === "pop" ? C.pop[k] : AREA[k]);

function distribution() {
  const nb = META.nbins, bins = new Float64Array(nb);
  for (let k = 0; k < N; k++) {
    if (!passes(k)) continue;
    const w = wt(k);
    if (w <= 0) continue;
    const h = C.hist[k], t0 = C.t0[k];
    for (let j = 0; j < h.length; j++) bins[t0 + j] += w * h[j];
  }
  let total = 0;                        // total exposure (person-hours / km²-hours)
  for (let b = 0; b < nb; b++) total += bins[b];
  return { bins, total };
}

const binTemp = (b) => META.tmin + (b + 0.5) * META.binw;   // bin center, Celsius

function fmtBig(x, unit) {
  const a = Math.abs(x);
  if (a >= 1e12) return (x / 1e12).toFixed(1) + " T " + unit;
  if (a >= 1e9) return (x / 1e9).toFixed(1) + " B " + unit;
  if (a >= 1e6) return (x / 1e6).toFixed(1) + " M " + unit;
  return Math.round(x).toLocaleString() + " " + unit;
}

// ---------------------------------------------------------------------------
// Map
// ---------------------------------------------------------------------------
function fit(cv, cssH) {
  const cssW = cv.clientWidth || cv.parentElement.clientWidth;
  const dpr = window.devicePixelRatio || 1;
  cv.width = Math.round(cssW * dpr);
  cv.height = Math.round(cssH * dpr);
  cv.style.height = cssH + "px";
  const ctx = cv.getContext("2d");
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  return { ctx, w: cssW, h: cssH };
}

const px = (lon, w) => ((lon + 180) / 360) * w;
const py = (lat, h) => ((90 - lat) / 180) * h;

function drawMap() {
  const cv = $("map");
  const H = cv.clientWidth / 2;
  const { ctx, w, h } = fit(cv, H);
  ctx.clearRect(0, 0, w, h);
  ctx.fillStyle = "#fbfbfb"; ctx.fillRect(0, 0, w, h);

  const deg = META.deg, cw = (deg / 360) * w + 0.6, ch = (deg / 180) * h + 0.6;

  // precompute normalisation for the active view
  let maxPop = 1, hi = null, hiMax = 1, hiCol = null;
  if (state.hoverBin != null) {
    hi = new Float64Array(N);
    for (let k = 0; k < N; k++) {
      if (!passes(k)) continue;
      const b = state.hoverBin, t0 = C.t0[k], j = b - t0;
      if (j >= 0 && j < C.hist[k].length) { hi[k] = C.hist[k][j] * wt(k); if (hi[k] > hiMax) hiMax = hi[k]; }
    }
    hiCol = tempColor(binTemp(state.hoverBin));
  } else if (state.mapview === "pop") {
    for (let k = 0; k < N; k++) if (C.pop[k] > maxPop) maxPop = C.pop[k];
  }
  const logMax = Math.log(1 + maxPop);

  for (let k = 0; k < N; k++) {
    const x = px(C.lon[k] - deg / 2, w), y = py(C.lat[k] + deg / 2, h);
    const ok = passes(k);
    if (hi) {
      if (!ok || hi[k] <= 0) { ctx.fillStyle = "#f0f0f0"; }
      else ctx.fillStyle = rgb(hiCol, 0.15 + 0.85 * Math.sqrt(hi[k] / hiMax));
    } else if (!ok) {
      ctx.fillStyle = "#ededed";
    } else if (state.mapview === "pop") {
      ctx.fillStyle = C.pop[k] > 0 ? popColor(Math.log(1 + C.pop[k]) / logMax) : "#f2f2f2";
    } else {
      ctx.fillStyle = rgb(tempColor(C.tmean[k]));
    }
    ctx.fillRect(x, y, cw, ch);
  }

  graticule(ctx, w, h);
  if (state.box) drawBox(ctx, w, h, state.box);
  updateLegend();
}

function graticule(ctx, w, h) {
  ctx.strokeStyle = "rgba(0,0,0,0.06)"; ctx.lineWidth = 1;
  for (let lon = -180; lon <= 180; lon += 30) {
    ctx.beginPath(); ctx.moveTo(px(lon, w), 0); ctx.lineTo(px(lon, w), h); ctx.stroke();
  }
  for (let lat = -60; lat <= 60; lat += 30) {
    ctx.beginPath(); ctx.moveTo(0, py(lat, h)); ctx.lineTo(w, py(lat, h)); ctx.stroke();
  }
  ctx.strokeStyle = "rgba(0,0,0,0.14)";
  ctx.beginPath(); ctx.moveTo(0, py(0, h)); ctx.lineTo(w, py(0, h)); ctx.stroke();
}

function drawBox(ctx, w, h, b) {
  const x = px(b.lonMin, w), y = py(b.latMax, h);
  const x2 = px(b.lonMax, w), y2 = py(b.latMin, h);
  ctx.strokeStyle = "#111"; ctx.lineWidth = 1.5;
  ctx.setLineDash([4, 3]); ctx.strokeRect(x, y, x2 - x, y2 - y); ctx.setLineDash([]);
}

// ---------------------------------------------------------------------------
// Chart
// ---------------------------------------------------------------------------
let chartGeom = null;

function drawChart() {
  const cv = $("chart");
  const { ctx, w, h } = fit(cv, 260);
  ctx.clearRect(0, 0, w, h);
  const { bins, total } = distribution();

  // visible temperature window: trim empty tails, clamp to a sensible range
  let lo = 0, hi = bins.length - 1;
  while (lo < hi && bins[lo] === 0) lo++;
  while (hi > lo && bins[hi] === 0) hi--;
  lo = Math.max(0, lo - 1); hi = Math.min(bins.length - 1, hi + 1);

  const padL = 44, padR = 8, padT = 10, padB = 26;
  const plotW = w - padL - padR, plotH = h - padT - padB;
  let maxShare = 0;
  for (let b = lo; b <= hi; b++) maxShare = Math.max(maxShare, bins[b] / total);
  maxShare = maxShare || 1;

  const bx = (b) => padL + ((b - lo) / (hi - lo + 1)) * plotW;
  const bw = plotW / (hi - lo + 1);
  chartGeom = { lo, hi, padL, padT, plotW, plotH, bw, h, w, padB };

  // y gridlines (share %)
  ctx.fillStyle = "#999"; ctx.font = "11px system-ui, sans-serif";
  ctx.textAlign = "right"; ctx.textBaseline = "middle";
  const ticks = niceTicks(maxShare, 4);
  ctx.strokeStyle = "#eee"; ctx.lineWidth = 1;
  for (const t of ticks) {
    const y = padT + plotH - (t / maxShare) * plotH;
    ctx.beginPath(); ctx.moveTo(padL, y); ctx.lineTo(w - padR, y); ctx.stroke();
    ctx.fillText((t * 100).toFixed(t < 0.01 ? 1 : 0) + "%", padL - 6, y);
  }

  // bars
  for (let b = lo; b <= hi; b++) {
    const share = bins[b] / total;
    const bh = (share / maxShare) * plotH;
    const x = bx(b), y = padT + plotH - bh;
    ctx.fillStyle = rgb(tempColor(binTemp(b)), b === state.hoverBin ? 1 : 0.9);
    ctx.fillRect(x + 0.4, y, Math.max(1, bw - 0.8), bh);
    if (b === state.hoverBin) { ctx.strokeStyle = "#111"; ctx.lineWidth = 1; ctx.strokeRect(x + 0.4, y, Math.max(1, bw - 0.8), bh); }
  }

  // x axis: temperature labels + 0C marker
  ctx.strokeStyle = "#ccc"; ctx.beginPath();
  ctx.moveTo(padL, padT + plotH); ctx.lineTo(w - padR, padT + plotH); ctx.stroke();
  ctx.fillStyle = "#666"; ctx.textAlign = "center"; ctx.textBaseline = "top";
  const step = (hi - lo) > 70 ? 20 : 10;
  for (let t = Math.ceil(binTemp(lo) / step) * step; t <= binTemp(hi); t += step) {
    const b = (t - META.tmin) / META.binw - 0.5;
    const x = bx(b) + bw / 2;
    ctx.fillStyle = t === 0 ? "#111" : "#666";
    ctx.fillText(t + "°", x, padT + plotH + 5);
    if (t === 0) { ctx.strokeStyle = "rgba(0,0,0,0.25)"; ctx.beginPath(); ctx.moveTo(x, padT); ctx.lineTo(x, padT + plotH); ctx.stroke(); }
  }
  ctx.fillStyle = "#999"; ctx.textAlign = "left";
  ctx.fillText("air temperature", padL, padT + plotH + 14);
}

function niceTicks(max, n) {
  const raw = max / n, mag = Math.pow(10, Math.floor(Math.log10(raw)));
  const step = [1, 2, 2.5, 5, 10].map((m) => m * mag).find((s) => s >= raw) || mag;
  const out = [];
  for (let t = step; t <= max * 1.0001; t += step) out.push(t);
  return out;
}

// ---------------------------------------------------------------------------
// Stats
// ---------------------------------------------------------------------------
function update() {
  drawMap();
  drawChart();
  const { bins } = distribution();
  const nb = META.nbins;
  // total exposure = person-hours (or km²-hours); this is the correct denominator
  let totPH = 0;
  for (let b = 0; b < nb; b++) totPH += bins[b];
  // weighted mean and median experienced temperature
  let mean = 0, cum = 0, median = null, above30 = 0, below0 = 0;
  for (let b = 0; b < nb; b++) mean += binTemp(b) * bins[b];
  mean = totPH ? mean / totPH : 0;
  for (let b = 0; b < nb; b++) {
    cum += bins[b];
    if (median == null && cum >= totPH / 2) median = binTemp(b);
    if (META.tmin + (b + 1) * META.binw > 30) above30 += bins[b];
    if (META.tmin + b * META.binw < 0) below0 += bins[b];
  }
  const unit = state.weight === "pop" ? "person-hours" : "km²-hours";
  const totLabel = state.weight === "pop"
    ? (selPop() / 1e9).toFixed(2) + " B people"
    : (selArea() / 1e6).toFixed(1) + " M km²";
  const pct = (x) => totPH ? (100 * x / totPH).toFixed(1) + "%" : "–";
  $("stats").innerHTML = stat(mean.toFixed(1) + "°C", "mean temperature experienced")
    + stat((median == null ? "–" : median.toFixed(0) + "°C"), "median")
    + stat(pct(above30), "of " + unit + " above 30°C")
    + stat(pct(below0), "below 0°C")
    + stat(totLabel, "in selection");
}
function stat(n, l) { return `<div class="stat"><div class="n">${n}</div><div class="l">${l}</div></div>`; }
function selPop() { let s = 0; for (let k = 0; k < N; k++) if (passes(k)) s += C.pop[k]; return s; }
function selArea() { let s = 0; for (let k = 0; k < N; k++) if (passes(k)) s += AREA[k]; return s; }

// ---------------------------------------------------------------------------
// Legend
// ---------------------------------------------------------------------------
function updateLegend() {
  const bar = $("legend-bar"), ticks = $("legend-ticks"), label = $("legend-label");
  if (state.mapview === "pop" && state.hoverBin == null) {
    label.textContent = "Population per cell";
    bar.style.background = "linear-gradient(90deg, #eee, #181818)";
    ticks.innerHTML = "<span>low</span><span>high</span>";
  } else {
    const t = state.hoverBin != null ? binTemp(state.hoverBin) : null;
    label.textContent = t != null
      ? `Where ${(t - META.binw / 2).toFixed(0)}–${(t + META.binw / 2).toFixed(0)}°C is felt`
      : "Mean temperature";
    const stops = [];
    for (let i = 0; i <= 10; i++) {
      const c = tempColor(TEMP_DOMAIN[0] + (i / 10) * (TEMP_DOMAIN[1] - TEMP_DOMAIN[0]));
      stops.push(rgb(c) + " " + i * 10 + "%");
    }
    bar.style.background = `linear-gradient(90deg, ${stops.join(",")})`;
    ticks.innerHTML = `<span>${TEMP_DOMAIN[0]}°</span><span>${(TEMP_DOMAIN[0] + TEMP_DOMAIN[1]) / 2 | 0}°</span><span>${TEMP_DOMAIN[1]}°+</span>`;
  }
}

// ---------------------------------------------------------------------------
// Controls + interaction
// ---------------------------------------------------------------------------
function buildControls() {
  seg("weight", (v) => { state.weight = v; update(); });
  seg("mapview", (v) => { state.mapview = v; drawMap(); });
  const box = $("continents");
  META.continents.forEach((name, i) => {
    const b = document.createElement("button");
    b.textContent = name; b.setAttribute("aria-pressed", "true");
    b.onclick = () => {
      if (state.continents.has(i)) state.continents.delete(i); else state.continents.add(i);
      if (state.continents.size === 0) state.continents.add(i); // never empty
      b.setAttribute("aria-pressed", state.continents.has(i));
      update();
    };
    box.appendChild(b);
  });
  const all = document.createElement("button");
  all.className = "mini"; all.textContent = "all";
  all.onclick = () => {
    META.continents.forEach((_, i) => state.continents.add(i));
    [...box.querySelectorAll("button[aria-pressed]")].forEach((x) => x.setAttribute("aria-pressed", "true"));
    update();
  };
  box.appendChild(all);

  $("clearbox").onclick = (e) => { e.preventDefault(); state.box = null; $("clearbox").style.display = "none"; update(); };
  wireMap();
  wireChart();
}

function seg(id, cb) {
  const el = $(id);
  el.querySelectorAll("button").forEach((b) => {
    b.onclick = () => {
      el.querySelectorAll("button").forEach((x) => x.setAttribute("aria-pressed", "false"));
      b.setAttribute("aria-pressed", "true");
      cb(b.dataset.v);
    };
  });
}

function evLonLat(cv, e) {
  const r = cv.getBoundingClientRect();
  const x = (e.clientX - r.left) / r.width, y = (e.clientY - r.top) / r.height;
  return { lon: x * 360 - 180, lat: 90 - y * 180 };
}

function wireMap() {
  const cv = $("map"), ro = $("readout");
  let drag = null;
  cv.addEventListener("mousedown", (e) => { drag = evLonLat(cv, e); });
  cv.addEventListener("mousemove", (e) => {
    const p = evLonLat(cv, e);
    if (drag) {
      state.box = norm(drag, p);
      drawMap();
    } else {
      const k = nearestCell(p.lon, p.lat);
      if (k >= 0) {
        ro.style.opacity = 1;
        ro.textContent = `${fmtLat(C.lat[k])}, ${fmtLon(C.lon[k])} · ${(C.pop[k] / 1e6).toFixed(2)} M · mean ${C.tmean[k].toFixed(1)}°C`;
      }
    }
  });
  cv.addEventListener("mouseleave", () => { ro.style.opacity = 0; });
  window.addEventListener("mouseup", (e) => {
    if (!drag) return;
    const p = evLonLat(cv, e), b = norm(drag, p);
    drag = null;
    if ((b.latMax - b.latMin) < META.deg || (b.lonMax - b.lonMin) < META.deg) {
      state.box = null; $("clearbox").style.display = "none"; // treated as a click -> clear
    } else {
      state.box = b; $("clearbox").style.display = "inline";
    }
    update();
  });
}
function norm(a, b) {
  return { latMin: Math.min(a.lat, b.lat), latMax: Math.max(a.lat, b.lat),
           lonMin: Math.min(a.lon, b.lon), lonMax: Math.max(a.lon, b.lon) };
}
function nearestCell(lon, lat) {
  let best = -1, bd = Infinity;
  for (let k = 0; k < N; k++) {
    const d = Math.abs(C.lat[k] - lat) + Math.abs(C.lon[k] - lon);
    if (d < bd) { bd = d; best = k; }
  }
  return bd < META.deg ? best : -1;
}
const fmtLat = (v) => Math.abs(v).toFixed(1) + (v >= 0 ? "°N" : "°S");
const fmtLon = (v) => Math.abs(v).toFixed(1) + (v >= 0 ? "°E" : "°W");

function wireChart() {
  const cv = $("chart"), ro = $("chart-readout");
  cv.addEventListener("mousemove", (e) => {
    if (!chartGeom) return;
    const r = cv.getBoundingClientRect();
    const x = e.clientX - r.left;
    const g = chartGeom;
    const b = g.lo + Math.floor((x - g.padL) / g.bw);
    if (b < g.lo || b > g.hi) { clearHover(); return; }
    state.hoverBin = b;
    drawMap(); drawChart();
    const { bins, total } = distribution();
    const share = total ? (100 * bins[b] / total) : 0;
    const abs = fmtBig(bins[b], state.weight === "pop" ? "person-hours" : "km²-hours");
    ro.style.opacity = 1;
    ro.style.left = Math.min(g.w - 150, x + 10) + "px";
    ro.textContent = `${(binTemp(b) - META.binw / 2).toFixed(0)} to ${(binTemp(b) + META.binw / 2).toFixed(0)}°C · ${share.toFixed(1)}% · ${abs}`;
  });
  cv.addEventListener("mouseleave", clearHover);
  function clearHover() {
    if (state.hoverBin == null) return;
    state.hoverBin = null; ro.style.opacity = 0; drawMap(); drawChart();
  }
}

load();
