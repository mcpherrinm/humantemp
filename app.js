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
const TEMP_DOMAIN = [-40, 40]; // symmetric about 0°C: blues negative, reds positive

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
function grayShade(frac) { // white(0) -> near-black(1); for non-temperature magnitudes
  const g = Math.round(255 - 235 * Math.max(0, Math.min(1, frac)));
  return `rgb(${g},${g},${g})`;
}
const RANGE_MAX = 18; // °C, top of the mean-daily-range grayscale

// map views. The temperature-coloured ones map to a per-cell °C field; 'range'
// and 'pop' are grayscale magnitudes handled separately.
const TEMPVIEW = { temp: "tmean", mdhigh: "mdhigh", mdlow: "mdlow", tmax: "tmax", tmin: "tmin" };
const VIEW_LABEL = {
  temp: "Mean temperature", mdhigh: "Mean daily high", mdlow: "Mean daily low",
  tmax: "Annual maximum", tmin: "Annual minimum",
  range: "Mean daily range (high − low)", pop: "Population per cell",
};

// ---------------------------------------------------------------------------
// Units. Data is stored in °C and both scales are shown together. Absolute
// temperatures use toF; a temperature *difference* (daily range) uses toFd
// (the ×9/5 slope, no +32 offset). cFromF maps a °F tick back to °C for
// positioning on the °C-based axes.
// ---------------------------------------------------------------------------
const toF = (c) => c * 9 / 5 + 32;
const toFd = (c) => c * 9 / 5;
const cFromF = (f) => (f - 32) * 5 / 9;
const fmtT = (c, dp = 0) => `${c.toFixed(dp)}°C / ${toF(c).toFixed(dp)}°F`;   // absolute
const fmtD = (c, dp = 0) => `${c.toFixed(dp)}°C / ${toFd(c).toFixed(dp)}°F`;  // difference
const fmtThtml = (c, dp = 0) => `${c.toFixed(dp)}°C <span class="alt">${toF(c).toFixed(dp)}°F</span>`;

// ---------------------------------------------------------------------------
// State
// ---------------------------------------------------------------------------
let META, C, N, DAILY;          // metadata, cells (columnar), cell count, daily 2D
let AREA;                       // per-cell area weight (km^2)
let state = {
  weight: "pop",                // 'pop' | 'area'
  mapview: "temp",              // 'temp' | 'range' | 'pop'
  continents: new Set(),        // active continent indices
  box: null,                    // {latMin,latMax,lonMin,lonMax}
  hoverBin: null,               // highlighted temperature bin (distribution hover)
  heatHover: null,              // {mn,mx} highlighted daily-range cell
  pct: 95,                      // upper percentile; band is (100-pct)..pct
};

const $ = (id) => document.getElementById(id);

// ---------------------------------------------------------------------------
// Load
// ---------------------------------------------------------------------------
async function load() {
  [META, C, DAILY] = await Promise.all([
    fetch("data/meta.json").then((r) => r.json()),
    fetch("data/cells.json").then((r) => r.json()),
    fetch("data/daily.json").then((r) => r.json()),
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
  window.addEventListener("resize", () => { drawMap(); drawChart(); drawHeat(); });
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

// temperature at cumulative fraction `frac` of a weighted histogram (interpolated)
function pctTemp(bins, total, frac) {
  const target = frac * total;
  let cum = 0;
  for (let b = 0; b < bins.length; b++) {
    if (cum + bins[b] >= target) {
      const pos = bins[b] > 0 ? (target - cum) / bins[b] : 0;
      return META.tmin + (b + pos) * META.binw;
    }
    cum += bins[b];
  }
  return META.tmax;
}
function ordinal(n) {
  const t = n % 100, u = n % 10;
  const s = (t >= 11 && t <= 13) ? "th" : u === 1 ? "st" : u === 2 ? "nd" : u === 3 ? "rd" : "th";
  return n + s;
}

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
    } else if (state.mapview === "range") {
      ctx.fillStyle = grayShade(C.drange[k] / RANGE_MAX);
    } else {
      ctx.fillStyle = rgb(tempColor(C[TEMPVIEW[state.mapview]][k]));
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

  const padL = 44, padR = 8, padT = 22, padB = 26;
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

  // freezing marker (0°C = 32°F), drawn at its bin position regardless of ticks
  const xpos = (c) => bx((c - META.tmin) / META.binw - 0.5) + bw / 2;
  ctx.strokeStyle = "rgba(0,0,0,0.22)"; ctx.beginPath();
  ctx.moveTo(xpos(0), padT); ctx.lineTo(xpos(0), padT + plotH); ctx.stroke();

  // x axes: °C along the bottom, °F along the top (both shown at once)
  const cLo = binTemp(lo), cHi = binTemp(hi);
  ctx.strokeStyle = "#ccc"; ctx.beginPath();
  ctx.moveTo(padL, padT + plotH); ctx.lineTo(w - padR, padT + plotH); ctx.stroke();
  ctx.fillStyle = "#666"; ctx.textAlign = "center"; ctx.textBaseline = "top";
  const cStep = (cHi - cLo) > 70 ? 20 : 10;
  for (let c = Math.ceil(cLo / cStep) * cStep; c <= cHi; c += cStep)
    ctx.fillText(c + "°", xpos(c), padT + plotH + 5);
  ctx.textBaseline = "bottom";
  const fLo = toF(cLo), fHi = toF(cHi), fStep = (fHi - fLo) > 90 ? 40 : 20;
  for (let f = Math.ceil(fLo / fStep) * fStep; f <= fHi; f += fStep)
    ctx.fillText(f + "°", xpos(cFromF(f)), padT - 4);
  ctx.fillStyle = "#999"; ctx.textAlign = "left"; ctx.textBaseline = "top";
  ctx.fillText("air temperature — °C below, °F above", padL, padT + plotH + 14);

  // symmetric percentile lines: slider P -> lines at (100-P) and P.
  // Labels point inward and sit on two rows so they never overlap or clip.
  const P = state.pct;
  ctx.font = "11px system-ui, sans-serif"; ctx.textBaseline = "top";
  const rows = [padT + 1, padT + 15];
  for (const [frac, pc, side, row] of [[(100 - P) / 100, 100 - P, 1, 0], [P / 100, P, -1, 1]]) {
    const T = pctTemp(bins, total, frac), x = xpos(T);
    ctx.strokeStyle = "#222"; ctx.lineWidth = 1.25;
    ctx.beginPath(); ctx.moveTo(x, padT); ctx.lineTo(x, padT + plotH); ctx.stroke();
    ctx.fillStyle = "#222"; ctx.textAlign = side > 0 ? "left" : "right";
    ctx.fillText(`${ordinal(pc)} · ${fmtT(T)}`, x + side * 4, rows[row]);
  }
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
  drawHeat();
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
  $("stats").innerHTML = stat(fmtThtml(mean, 1), "mean temperature experienced")
    + stat((median == null ? "–" : fmtThtml(median)), "median")
    + stat(pct(above30), "of " + unit + " above " + fmtT(30))
    + stat(pct(below0), "below " + fmtT(0))
    + stat(totLabel, "in selection");
}
function stat(n, l) { return `<div class="stat"><div class="n">${n}</div><div class="l">${l}</div></div>`; }
function selPop() { let s = 0; for (let k = 0; k < N; k++) if (passes(k)) s += C.pop[k]; return s; }
function selArea() { let s = 0; for (let k = 0; k < N; k++) if (passes(k)) s += AREA[k]; return s; }

// ---------------------------------------------------------------------------
// Legend
// ---------------------------------------------------------------------------
const cfTick = (c) => `<span>${c | 0}°C<b>${toF(c).toFixed(0)}°F</b></span>`;   // stacked °C / °F

function tempLegend(bar, ticks) {
  const stops = [];
  for (let i = 0; i <= 10; i++) {
    const c = tempColor(TEMP_DOMAIN[0] + (i / 10) * (TEMP_DOMAIN[1] - TEMP_DOMAIN[0]));
    stops.push(rgb(c) + " " + i * 10 + "%");
  }
  bar.style.background = `linear-gradient(90deg, ${stops.join(",")})`;
  const a = TEMP_DOMAIN[0], b = TEMP_DOMAIN[1];
  ticks.innerHTML = cfTick(a) + cfTick((a + b) / 2) + cfTick(b).replace("°F", "°F+");
}

function updateLegend() {
  const bar = $("legend-bar"), ticks = $("legend-ticks"), label = $("legend-label");
  if (state.hoverBin != null) {
    const t = binTemp(state.hoverBin), lo = t - META.binw / 2, hi = t + META.binw / 2;
    label.textContent = `Where ${lo.toFixed(0)}–${hi.toFixed(0)}°C / ${toF(lo).toFixed(0)}–${toF(hi).toFixed(0)}°F is felt`;
    tempLegend(bar, ticks);
  } else if (state.mapview === "pop") {
    label.textContent = "Population per cell";
    bar.style.background = "linear-gradient(90deg, #fff, #181818)";
    ticks.innerHTML = "<span>low</span><span>high</span>";
  } else if (state.mapview === "range") {
    label.textContent = "Mean daily range (high − low)";
    bar.style.background = "linear-gradient(90deg, #fff, #141414)";
    ticks.innerHTML = `<span>0°</span><span>${fmtD(RANGE_MAX)}+</span>`;
  } else {
    label.textContent = VIEW_LABEL[state.mapview] || "Mean temperature";
    tempLegend(bar, ticks);
  }
}

// ---------------------------------------------------------------------------
// Controls + interaction
// ---------------------------------------------------------------------------
function buildControls() {
  seg("weight", (v) => { state.weight = v; update(); });
  $("mapview").onchange = (e) => { state.mapview = e.target.value; drawMap(); };
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

  const slider = $("pct");
  slider.value = state.pct;
  slider.oninput = () => { state.pct = +slider.value; updatePctLabel(); drawChart(); };
  updatePctLabel();

  $("clearbox").onclick = (e) => { e.preventDefault(); state.box = null; $("clearbox").style.display = "none"; update(); };
  wireMap();
  wireChart();
  wireHeat();
}

function updatePctLabel() {
  $("pct-label").textContent = `${ordinal(100 - state.pct)} – ${ordinal(state.pct)} percentile`;
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
        ro.textContent = `${fmtLat(C.lat[k])}, ${fmtLon(C.lon[k])} · ${(C.pop[k] / 1e6).toFixed(2)} M · ${metricAt(k)}`;
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

function metricAt(k) { // the currently-mapped statistic, for the hover readout
  const v = state.mapview;
  if (v === "range") return `daily range ${fmtD(C.drange[k], 1)}`;
  const label = { temp: "mean", mdhigh: "mean daily high", mdlow: "mean daily low",
                  tmax: "annual max", tmin: "annual min" }[v] || "mean";
  return `${label} ${fmtT(C[TEMPVIEW[v] || "tmean"][k], 1)}`;
}

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
    ro.style.left = Math.min(g.w - 160, x + 10) + "px";
    const clo = binTemp(b) - META.binw / 2, chi = binTemp(b) + META.binw / 2;
    ro.textContent = `${clo.toFixed(0)}–${chi.toFixed(0)}°C / ${toF(clo).toFixed(0)}–${toF(chi).toFixed(0)}°F · ${share.toFixed(1)}% · ${abs}`;
  });
  cv.addEventListener("mouseleave", clearHover);
  function clearHover() {
    if (state.hoverBin == null) return;
    state.hoverBin = null; ro.style.opacity = 0; drawMap(); drawChart();
  }
}

// ---------------------------------------------------------------------------
// Daily min/max heatmap: 2D histogram of (daily low, daily high) shaded by
// people-days. Grayscale on purpose — this is a count, not a temperature.
// ---------------------------------------------------------------------------
const binTemp2 = (b) => META.tmin + (b + 0.5) * META.binw2;
let heatState = null;

function heatBins() {
  const nb2 = META.nbins2, bins = new Float64Array(nb2 * nb2);
  let max = 0, total = 0;
  for (let k = 0; k < N; k++) {
    if (!passes(k)) continue;
    const w = wt(k);
    if (w <= 0) continue;
    const kv = DAILY.cell[k];
    for (let p = 0; p < kv.length; p += 2) {
      const v = (bins[kv[p]] += w * kv[p + 1]);
      if (v > max) max = v;
      total += w * kv[p + 1];
    }
  }
  return { bins, max, total };
}

function drawHeat() {
  const cv = $("heat");
  const nb2 = META.nbins2;
  const { bins, max } = heatBins();
  const padL = 42, padR = 32, padT = 24, padB = 30;
  const side = (cv.clientWidth || cv.parentElement.clientWidth) - padL - padR;
  const { ctx, w, h } = fit(cv, side + padT + padB);
  ctx.clearRect(0, 0, w, h);

  // occupied temperature window (shared by both axes so the diagonal is 45°)
  let gLo = nb2, gHi = 0;
  for (let i = 0; i < bins.length; i++) if (bins[i] > 0) {
    const mn = (i / nb2) | 0, mx = i % nb2;
    if (mn < gLo) gLo = mn; if (mx > gHi) gHi = mx; if (mn > gHi) gHi = mn; if (mx < gLo) gLo = mx;
  }
  if (gLo > gHi) { heatState = null; updateHeatLegend(); return; }
  gLo = Math.max(0, gLo - 1); gHi = Math.min(nb2 - 1, gHi + 1);
  const span = gHi - gLo + 1, cell = side / span;
  const xOf = (c) => padL + ((c - META.tmin) / META.binw2 - gLo) * cell;
  const yOf = (c) => padT + side - ((c - META.tmin) / META.binw2 - gLo) * cell;
  heatState = { bins, max, gLo, gHi, cell, padL, padT, side, w };
  const lmax = Math.log(1 + max);

  // cells
  for (let mn = gLo; mn <= gHi; mn++) {
    for (let mx = mn; mx <= gHi; mx++) {
      const v = bins[mn * nb2 + mx];
      if (v <= 0) continue;
      const f = Math.pow(Math.log(1 + v) / lmax, 0.6);
      ctx.fillStyle = grayShade(f);
      ctx.fillRect(padL + (mn - gLo) * cell, padT + side - (mx - gLo + 1) * cell, cell + 0.5, cell + 0.5);
    }
  }
  // diagonal (no daily swing) + plot border
  ctx.strokeStyle = "rgba(0,0,0,0.28)"; ctx.setLineDash([4, 3]);
  ctx.beginPath(); ctx.moveTo(xOf(binTemp2(gLo)), yOf(binTemp2(gLo))); ctx.lineTo(xOf(binTemp2(gHi)), yOf(binTemp2(gHi))); ctx.stroke();
  ctx.setLineDash([]);
  ctx.strokeStyle = "#ddd"; ctx.strokeRect(padL, padT, side, side);

  // dual axes: °C on bottom + left (with gridlines), °F on top + right
  ctx.font = "11px system-ui, sans-serif";
  const cLo = binTemp2(gLo), cHi = binTemp2(gHi);
  for (let c = Math.ceil(cLo / 10) * 10; c <= cHi; c += 10) {
    const gx = xOf(c), gy = yOf(c);
    ctx.strokeStyle = c === 0 ? "rgba(0,0,0,0.18)" : "#f0f0f0";
    ctx.beginPath(); ctx.moveTo(gx, padT); ctx.lineTo(gx, padT + side); ctx.stroke();
    ctx.beginPath(); ctx.moveTo(padL, gy); ctx.lineTo(padL + side, gy); ctx.stroke();
    ctx.fillStyle = "#666";
    ctx.textAlign = "center"; ctx.textBaseline = "top"; ctx.fillText(c + "°", gx, padT + side + 5);
    ctx.textAlign = "right"; ctx.textBaseline = "middle"; ctx.fillText(c + "°", padL - 5, gy);
  }
  ctx.fillStyle = "#aaa";
  const fLo = toF(cLo), fHi = toF(cHi);
  for (let f = Math.ceil(fLo / 20) * 20; f <= fHi; f += 20) {
    const c = cFromF(f), gx = xOf(c), gy = yOf(c);
    ctx.textAlign = "center"; ctx.textBaseline = "bottom"; ctx.fillText(f + "°", gx, padT - 4);
    ctx.textAlign = "left"; ctx.textBaseline = "middle"; ctx.fillText(f + "°", padL + side + 5, gy);
  }
  ctx.fillStyle = "#999"; ctx.textAlign = "center"; ctx.textBaseline = "top";
  ctx.fillText("daily low — °C bottom, °F top", padL + side / 2, padT + side + 17);
  ctx.save(); ctx.translate(11, padT + side / 2); ctx.rotate(-Math.PI / 2);
  ctx.textBaseline = "middle"; ctx.fillText("daily high — °C left, °F right", 0, 0); ctx.restore();
  updateHeatLegend();
}

function updateHeatLegend() {
  $("heat-legend-bar").style.background = "linear-gradient(90deg, #fff, #141414)";
  $("heat-legend-label").textContent = state.weight === "pop" ? "People-days" : "km²-days";
}

function wireHeat() {
  const cv = $("heat"), ro = $("heat-readout");
  cv.addEventListener("mousemove", (e) => {
    if (!heatState) return;
    const g = heatState, r = cv.getBoundingClientRect();
    const x = e.clientX - r.left, y = e.clientY - r.top;
    const mn = g.gLo + Math.floor((x - g.padL) / g.cell);
    const mx = g.gLo + Math.floor((g.padT + g.side - y) / g.cell);
    if (mn < g.gLo || mn > g.gHi || mx < g.gLo || mx > g.gHi || mx < mn) { ro.style.opacity = 0; return; }
    const v = g.bins[mn * META.nbins2 + mx];
    const lo = binTemp2(mn), hi = binTemp2(mx);
    ro.style.opacity = 1;
    ro.style.left = Math.min(g.w - 170, x + 12) + "px";
    ro.style.top = Math.max(4, y - 34) + "px";
    ro.textContent = `low ${lo.toFixed(0)}°C/${toF(lo).toFixed(0)}°F, high ${hi.toFixed(0)}°C/${toF(hi).toFixed(0)}°F · swing ~${fmtD(hi - lo)} · ${fmtBig(v, state.weight === "pop" ? "people-days" : "km²-days")}`;
  });
  cv.addEventListener("mouseleave", () => { ro.style.opacity = 0; });
}

load();
