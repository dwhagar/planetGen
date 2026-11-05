// html/static/galaxystrip.js
//
// MAP.375: the linked side strip for height on the Galaxy Map. An edge-on strip of the same stars the map draws,
// with the distance from the galactic center across and the height above or below the galactic plane up and down,
// so what the 3D view cannot show (how high a star lies, and how thin the disk is) is read at a glance. It is an
// option (Menu, "Side strip"), drawn on a 2D canvas under the map, and a drawer that slides up from the bottom
// edge on a phone. A star picked in the strip is selected on the map and the map's own selection rings in the strip:
// the two views share one selection (the map's `highlightPosition` and `showStarInfo`).
//
// The strip shows the stars of one WEDGE of the galaxy, a pie slice of `WEDGE_DEG` degrees from the galactic center
// (Boss, 5:48 PM EST on 2026-11-04: the same stars as the plane view's arc wedge, MAP.371), centered on the
// direction of the map's focus unless the page sets the center (`setWedge`, which MAP.371 will call so the two
// views agree); a selection outside the wedge turns it to face the selection, so what is selected is always in the strip. Heights are drawn HEIGHT_EXAGGERATION times taller than true against the distance across (the
// disk is about 3 percent as thick as it is wide); the readout says the exact height and how many scale
// heights of the thin and the thick disk it is. The arithmetic below is free of the page, so
// tests/js/galaxystrip.test.mjs holds it; `createSideStrip` is the page's part.
//
//     azimuth = atan2(y, x)       R = hypot(x, y)       height = z       (pc from the galactic center)
//     x_px = pad + (R - R_min) * pxPerPc        y_px = middle - z * pxPerPc * HEIGHT_EXAGGERATION

export const HEIGHT_EXAGGERATION = 2.2;
export const WEDGE_DEG = 60;
// The strip draws at most this many stars (the brightest first); 2,000 on a phone.
export const STRIP_STARS = 6000;
export const STRIP_STARS_PHONE = 2000;
// The pointer's reach on a star, px: the map's own (a mouse 10, a finger 22).
export const STRIP_MOUSE_REACH_PX = 10;
export const STRIP_TOUCH_REACH_PX = 22;
// The room kept round the plot, px.
export const STRIP_PAD_PX = 28;
// The least the strip's width covers across, pc, so a single star does not blow up to the whole panel.
export const MIN_ACROSS_PC = 200;

export function azimuthDeg(x, y) {
  const a = (Math.atan2(y, x) * 180) / Math.PI;
  return a < 0 ? a + 360 : a;
}

// The angle from `a` to `b` in degrees, -180 to 180.
export function angleDifferenceDeg(a, b) {
  let d = (b - a) % 360;
  if (d > 180) d -= 360;
  if (d <= -180) d += 360;
  return d;
}

// Whether a point of the galaxy lies in the wedge `{ centerDeg, widthDeg }` (always, with no wedge).
export function inWedge(x, y, wedge) {
  if (!wedge) return true;
  if (Math.hypot(x, y) < 1e-9) return true;
  return Math.abs(angleDifferenceDeg(wedge.centerDeg, azimuthDeg(x, y))) <= wedge.widthDeg / 2 + 1e-9;
}

// The wedge of `widthDeg` centered on the direction of the point (x, y) (the focus), or on `fallbackDeg` when it is at the center.
export function wedgeAbout(x, y, widthDeg = WEDGE_DEG, fallbackDeg = 0) {
  const center = Math.hypot(x, y) < 1e-6 ? fallbackDeg : azimuthDeg(x, y);
  return { centerDeg: center, widthDeg: widthDeg };
}

// The stars of the wedge, as strip entries `{ star, r, h }`, the brightest first and at most `limit`.
export function stripEntries(stars, wedge, limit = STRIP_STARS) {
  const found = [];
  for (let i = 0; i < stars.length; i++) {
    const s = stars[i];
    if (!inWedge(s.x, s.y, wedge)) continue;
    found.push({ star: s, r: Math.hypot(s.x, s.y), h: s.z });
  }
  found.sort((p, q) => (q.star.luminosity_sol || 0) - (p.star.luminosity_sol || 0));
  return found.length > limit ? found.slice(0, limit) : found;
}

// The frame of a strip `width` by `height` px: the distance across it covers (the entries' range, with the
// `include` points too, never under MIN_ACROSS_PC) and the scales both ways.
export function stripFrame(entries, width, height, include) {
  let lo = Infinity;
  let hi = -Infinity;
  const take = (r) => {
    if (r < lo) lo = r;
    if (r > hi) hi = r;
  };
  for (let i = 0; i < entries.length; i++) take(entries[i].r);
  if (include) include.forEach((point) => take(Math.hypot(point.x, point.y)));
  if (!Number.isFinite(lo)) {
    lo = 0;
    hi = MIN_ACROSS_PC;
  }
  if (hi - lo < MIN_ACROSS_PC) {
    const mid = (hi + lo) / 2;
    lo = mid - MIN_ACROSS_PC / 2;
    hi = mid + MIN_ACROSS_PC / 2;
  }
  const plotW = Math.max(1, width - 2 * STRIP_PAD_PX);
  const pxPerPc = plotW / (hi - lo);
  const middle = height / 2;
  const frame = {
    rMin: lo, rMax: hi, pxPerPc: pxPerPc, width: width, height: height, middle: middle,
    x: (r) => STRIP_PAD_PX + (r - lo) * pxPerPc,
    y: (h) => middle - h * pxPerPc * HEIGHT_EXAGGERATION,
    r: (px) => lo + (px - STRIP_PAD_PX) / pxPerPc,
    h: (py) => (middle - py) / (pxPerPc * HEIGHT_EXAGGERATION),
  };
  return frame;
}

// Whether the point (r, h) is drawn inside the frame's plot.
export function inFrame(frame, r, h) {
  const px = frame.x(r);
  const py = frame.y(h);
  return px >= 0 && px <= frame.width && py >= 0 && py <= frame.height;
}

// The entry nearest the pixel (px, py) within `reach` px, the brightest of those level with it, or null.
export function pickInStrip(entries, frame, px, py, reach) {
  let best = null;
  let bestD = Infinity;
  for (let i = 0; i < entries.length; i++) {
    const e = entries[i];
    const d = Math.hypot(frame.x(e.r) - px, frame.y(e.h) - py);
    if (d > reach) continue;
    if (d < bestD - 0.5 || (best && Math.abs(d - bestD) <= 0.5 && (e.star.luminosity_sol || 0) > (best.star.luminosity_sol || 0))) {
      best = e;
      bestD = d;
    }
  }
  return best;
}

// How high a point lies: the exact height and its size in scale heights of the thin and of the thick disk (null without a shape).
export function heightReadout(heightPc, shape) {
  const thin = shape && shape.disk_scale_height_pc > 0 ? Math.abs(heightPc) / shape.disk_scale_height_pc : null;
  const thick = shape && shape.thick_disk_scale_height_pc > 0 ? Math.abs(heightPc) / shape.thick_disk_scale_height_pc : null;
  return { heightPc: heightPc, thin: thin, thick: thick };
}

function oneDecimal(v) {
  return (Math.round(v * 10) / 10).toString();
}

// The readout as a line: "Height +134 pc above the plane (1.0 thin-disk and 0.2 thick-disk scale heights), 8.2 kpc from the center".
export function readoutText(r, heightPc, shape) {
  const read = heightReadout(heightPc, shape);
  const side = heightPc > 0 ? " above the plane" : heightPc < 0 ? " below the plane" : " on the plane";
  const sign = heightPc > 0 ? "+" : heightPc < 0 ? "-" : "";
  let text = "Height " + sign + Math.round(Math.abs(heightPc)).toLocaleString("en-US") + " pc" + side;
  if (read.thin !== null && read.thick !== null) {
    text += " (" + oneDecimal(read.thin) + " thin-disk and " + oneDecimal(read.thick) + " thick-disk scale heights)";
  }
  const across = r >= 1000 ? oneDecimal(r / 1000) + " kpc" : Math.round(r) + " pc";
  return text + ", " + across + " from the center";
}

// A round step (1, 2 or 5 times a power of ten) giving about `count` ticks across `span`.
export function niceStep(span, count = 6) {
  const raw = span / Math.max(1, count);
  const pow = Math.pow(10, Math.floor(Math.log10(raw)));
  const share = raw / pow;
  return (share < 1.5 ? 1 : share < 3.5 ? 2 : share < 7.5 ? 5 : 10) * pow;
}

// The ticks of a distance range: [{ r, label }] at round steps ("6 kpc" or "300 pc").
export function distanceTicks(rMin, rMax, count = 6) {
  const step = niceStep(rMax - rMin, count);
  const ticks = [];
  for (let r = Math.ceil(rMin / step) * step; r <= rMax + 1e-9; r += step) {
    ticks.push({ r: r, label: step >= 500 ? oneDecimal(r / 1000).replace(/\.0$/, "") + " kpc" : Math.round(r) + " pc" });
  }
  return ticks;
}

// ---- The page's part ----

const STAR_RGB = "238,240,255";

// The strip on the page. `host` (the element the strip fills), `getStars()` (the stars the map draws), `getShape()` (the
// density shape or null), `getFocus()` ({ x, y, z } the view's center), `onPick(star)` (a star picked in the strip), and
// `isPhone()`. Returns { setOn, isOn, setSelected, setWedge, redraw, update, entries, pickAt }.
export function createSideStrip(options) {
  const host = options.host;
  const doc = host.ownerDocument;
  const canvas = doc.createElement("canvas");
  canvas.className = "galaxy-strip-canvas";
  canvas.setAttribute("role", "img");
  canvas.setAttribute("aria-label", "Edge-on strip of the stars of the wedge: distance from the galactic center across, height up and down. Click a star to select it on the map.");
  const caption = doc.createElement("p");
  caption.className = "galaxy-strip-caption hint";
  caption.setAttribute("aria-live", "polite");
  const readout = doc.createElement("p");
  readout.className = "galaxy-strip-readout";
  readout.setAttribute("aria-live", "polite");
  host.appendChild(caption);
  host.appendChild(canvas);
  host.appendChild(readout);

  let on = false;
  let selected = null;       // { x, y, z } the map's selection
  let wedgeOverride = null;
  let cached = { signature: "", entries: [], frame: null, wedge: null };
  let lastSignature = "";

  // The wedge now: the page's (`setWedge`), else the one a selection outside the focus's wedge turned it to (kept while that
  // selection stands, so what is selected is always in the strip), else the one facing the focus.
  let turnedTo = null;
  function facingWedge() {
    const focus = options.getFocus ? options.getFocus() : { x: 0, y: 0, z: 0 };
    return wedgeAbout(focus.x, focus.y, WEDGE_DEG, cached.wedge ? cached.wedge.centerDeg : 0);
  }
  function wedgeNow() {
    return wedgeOverride || turnedTo || facingWedge();
  }

  function size() {
    const rect = host.getBoundingClientRect();
    const width = Math.max(120, Math.round(canvas.clientWidth || rect.width || 600));
    const height = Math.max(60, Math.round(canvas.clientHeight || 160));
    return { width: width, height: height };
  }

  function signatureOf(wedge, dims) {
    const focus = options.getFocus ? options.getFocus() : { x: 0, y: 0, z: 0 };
    return [on, options.version ? options.version() : 0, Math.round(wedge.centerDeg * 4), dims.width, dims.height,
      selected ? [selected.x, selected.y, selected.z].map(Math.round).join(",") : "-", Math.round(focus.x / 20), Math.round(focus.y / 20)].join("|");
  }

  function draw() {
    const dims = size();
    const wedge = wedgeNow();
    const phone = options.isPhone ? options.isPhone() : false;
    const ratio = Math.min(2, (doc.defaultView && doc.defaultView.devicePixelRatio) || 1);
    if (canvas.width !== Math.round(dims.width * ratio) || canvas.height !== Math.round(dims.height * ratio)) {
      canvas.width = Math.round(dims.width * ratio);
      canvas.height = Math.round(dims.height * ratio);
    }
    const stars = options.getStars ? options.getStars() : [];
    const entries = stripEntries(stars, wedge, phone ? STRIP_STARS_PHONE : STRIP_STARS);
    const frame = stripFrame(entries, dims.width, dims.height, selected ? [selected] : null);
    cached = { signature: signatureOf(wedge, dims), entries: entries, frame: frame, wedge: wedge };
    const ctx = canvas.getContext("2d");
    if (!ctx) return;
    ctx.setTransform(ratio, 0, 0, ratio, 0, 0);
    ctx.clearRect(0, 0, dims.width, dims.height);
    ctx.fillStyle = "rgba(8,12,28,0.92)";
    ctx.fillRect(0, 0, dims.width, dims.height);
    // The distance ticks, and the plane through the middle.
    ctx.font = "10px sans-serif";
    ctx.textBaseline = "bottom";
    ctx.strokeStyle = "rgba(140,160,220,0.16)";
    ctx.fillStyle = "rgba(170,185,230,0.8)";
    distanceTicks(frame.rMin, frame.rMax).forEach((tick) => {
      const px = frame.x(tick.r);
      ctx.beginPath();
      ctx.moveTo(px, 0);
      ctx.lineTo(px, dims.height - 12);
      ctx.stroke();
      ctx.fillText(tick.label, Math.max(2, px - 10), dims.height - 1);
    });
    ctx.strokeStyle = "rgba(160,180,240,0.5)";
    ctx.beginPath();
    ctx.moveTo(0, frame.middle);
    ctx.lineTo(dims.width, frame.middle);
    ctx.stroke();
    // The disk's thickness: one scale height of the thin disk (dashed) and of the thick disk (dotted) either side.
    const shape = options.getShape ? options.getShape() : null;
    if (shape && ctx.setLineDash) {
      [[shape.disk_scale_height_pc, [6, 4]], [shape.thick_disk_scale_height_pc, [2, 4]]].forEach((row) => {
        if (!(row[0] > 0)) return;
        ctx.setLineDash(row[1]);
        ctx.strokeStyle = "rgba(210,190,140,0.45)";
        [1, -1].forEach((side) => {
          const py = frame.y(side * row[0]);
          if (py < 0 || py > dims.height) return;
          ctx.beginPath();
          ctx.moveTo(0, py);
          ctx.lineTo(dims.width, py);
          ctx.stroke();
        });
      });
      ctx.setLineDash([]);
    }
    // The stars, the brightest last so they lie on top.
    let drawn = 0;
    for (let i = entries.length - 1; i >= 0; i--) {
      const e = entries[i];
      const px = frame.x(e.r);
      const py = frame.y(e.h);
      if (py < 0 || py > dims.height) continue;
      drawn++;
      const lum = e.star.luminosity_sol || 0;
      const radius = lum > 100 ? 2.2 : lum > 1 ? 1.6 : 1.1;
      ctx.fillStyle = "rgba(" + STAR_RGB + "," + (lum > 1 ? 0.95 : 0.6) + ")";
      ctx.beginPath();
      ctx.arc(px, py, radius, 0, 2 * Math.PI);
      ctx.fill();
    }
    // The map's selection, ringed.
    let inside = false;
    if (selected) {
      const r = Math.hypot(selected.x, selected.y);
      inside = inWedge(selected.x, selected.y, wedge) && inFrame(frame, r, selected.z);
      if (inside) {
        ctx.strokeStyle = options.accent || "#f0b44c";
        ctx.lineWidth = 2;
        ctx.beginPath();
        ctx.arc(frame.x(r), frame.y(selected.z), 8, 0, 2 * Math.PI);
        ctx.stroke();
        ctx.lineWidth = 1;
      }
    }
    caption.textContent = "Side strip: the " + Math.round(wedge.widthDeg) + " degree wedge toward " + Math.round(wedge.centerDeg)
      + " degrees from the galactic center, " + drawn.toLocaleString("en-US") + " stars, heights drawn " + HEIGHT_EXAGGERATION + " times true.";
    readout.textContent = selected
      ? readoutText(Math.hypot(selected.x, selected.y), selected.z, shape) + (inside ? "" : " (outside this wedge)")
      : "Click a star in the strip to select it on the map.";
    canvas.dataset.stars = String(drawn);
    canvas.dataset.wedge = String(Math.round(wedge.centerDeg));
    canvas.dataset.selected = inside ? "ringed" : "none";
    lastSignature = cached.signature;
  }

  function pointerPx(event) {
    const rect = canvas.getBoundingClientRect();
    return { x: event.clientX - rect.left, y: event.clientY - rect.top };
  }

  // The star under (px, py) of the strip, or null.
  function pickAt(px, py, touch) {
    if (!cached.frame) draw();
    const hit = pickInStrip(cached.entries, cached.frame, px, py, touch ? STRIP_TOUCH_REACH_PX : STRIP_MOUSE_REACH_PX);
    return hit ? hit.star : null;
  }

  canvas.addEventListener("click", (event) => {
    const p = pointerPx(event);
    const star = pickAt(p.x, p.y, event.pointerType === "touch" || (event.sourceCapabilities && event.sourceCapabilities.firesTouchEvents));
    if (star && options.onPick) options.onPick(star);
  });

  const api = {
    setOn(value) {
      on = !!value;
      host.hidden = !on;
      host.dataset.state = on ? "open" : "closed";
      if (on) draw();
    },
    isOn() { return on; },
    // The map's selection moved: ring it in the strip (null clears).
    setSelected(point) {
      selected = point ? { x: point.x, y: point.y, z: point.z } : null;
      if (!selected || inWedge(selected.x, selected.y, facingWedge())) {
        turnedTo = null;
      } else if (!turnedTo || !inWedge(selected.x, selected.y, turnedTo)) {
        turnedTo = { centerDeg: azimuthDeg(selected.x, selected.y), widthDeg: WEDGE_DEG };
      }
      if (on) draw();
    },
    // The wedge to show (MAP.371's, so the two views agree); null returns to following the focus.
    setWedge(wedge) {
      wedgeOverride = wedge || null;
      if (on) draw();
    },
    redraw() { if (on) draw(); },
    // Called each frame: redraws only when what the strip shows has changed.
    update() {
      if (!on) return;
      const dims = size();
      if (signatureOf(wedgeNow(), dims) !== lastSignature) draw();
    },
    entries() { return cached.entries; },
    // What the strip draws, for the tests: the wedge, and each star's place on the strip's canvas in CSS pixels.
    snapshot() {
      if (!cached.frame) return null;
      const frame = cached.frame;
      return {
        on: on, wedge: cached.wedge, rMin: frame.rMin, rMax: frame.rMax,
        stars: cached.entries.map((e) => ({ id: e.star.id, name: e.star.name || null, luminosity: e.star.luminosity_sol, r: e.r, h: e.h,
          x: frame.x(e.r), y: frame.y(e.h), visible: frame.y(e.h) >= 0 && frame.y(e.h) <= frame.height })),
      };
    },
    pickAt: pickAt,
    canvas: canvas,
  };
  host.hidden = true;
  return api;
}
