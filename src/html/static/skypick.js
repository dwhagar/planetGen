// html/static/skypick.js
//
// VIEW.37 (S19 of the System Sky plan, reports/system-sky-from-the-system-2026-11-04.md): picking a star of the sky a system sees, and what is said of
// it. The stars are drawn by the graphics card (skysphere.js), so a pick is worked out here, from the same arithmetic: each star's direction is its
// offset plus its motion since the list's epoch, then pushed back by the light time (`retardedOffset`), and that direction is turned into the camera's
// axes and put on the screen by the camera's projection. The star nearest the pointer within a few pixels is the pick, the brighter one when two are as near.
//
// A star carries its place in the list the server built (`sky.index`); `meta.charted` (docs/design/interface/system-sky-wire-format.md section 3.3) names
// the charted ones by it, with the printed star ID the page links to. A scatter, cluster or compact-object star has no ID to give until GEN.1254 and
// GEN.1255 land, so it is told but not linked.
//
// No DOM and no three.js here: tests/js/skypick.test.mjs holds it to this description.

const VERSION_QUERY = new URL(import.meta.url).search;
const { retardedOffset, KMS_TO_PC_PER_YEAR } = await import(`./skysphere.js${VERSION_QUERY}`);
const { apparentMagnitude } = await import(`./starmagnitude.js${VERSION_QUERY}`);

// How near (CSS pixels) the pointer must be to a star to pick it.
export const REACH_PX = 10;
// The first pass keeps a star whose direction without the light time is within this many pixels of the pointer more than REACH_PX (the light time moves
// a star by far less than this for all but the fastest).
export const PREFILTER_SLACK_PX = 30;
export const PC_TO_LY = 3.261563777;

// The spectral class of a temperature, by the usual cuts; null for none.
export function spectralClass(temperatureK) {
  if (!Number.isFinite(temperatureK)) return null;
  const cuts = [[30000, "O"], [10000, "B"], [7500, "A"], [6000, "F"], [5200, "G"], [3700, "K"]];
  for (const [limit, letter] of cuts) if (temperatureK >= limit) return letter;
  return "M";
}

const KIND_LABELS = ["Charted star", "Star not yet charted", "Compact object", "Star of a globular cluster"];

// "seen as it was N years ago", in words.
export function agoText(years) {
  if (!(years >= 0)) return "";
  if (years < 1) return "seen as it was less than a year ago";
  const whole = Math.round(years);
  return "seen as it was " + whole.toLocaleString("en-US") + (whole === 1 ? " year ago" : " years ago");
}

// The star of `sky` (skypack.js `unpackSky`, or skysphere.js `keepBrightest`) nearest the pointer, or null. `years` is the sky's own time after the list's
// epoch (`skyYears`); `view` is `{ matrix: [9 numbers, row-major, the galaxy's axes into the camera's], p00, p11 (the projection's x and y scales),
// width, height }` in CSS pixels; `pointer` is `{ x, y }` from the canvas's top left. `options`: `limit` (the limiting magnitude: a star dimmer is not
// drawn, so it is not picked), `reachPx`.
// Returns `{ place (in the sky's arrays), index, kind, distancePc, tauYears, magnitude, temperatureK, screenPx }`.
export function pickSkyStar(sky, years, view, pointer, options) {
  const o = options || {};
  const reach = Number.isFinite(o.reachPx) ? o.reachPx : REACH_PX;
  const limit = Number.isFinite(o.limit) ? o.limit : Infinity;
  const m = view.matrix;
  const toScreen = (x, y, z) => {
    const vx = m[0] * x + m[1] * y + m[2] * z;
    const vy = m[3] * x + m[4] * y + m[5] * z;
    const vz = m[6] * x + m[7] * y + m[8] * z;
    if (!(vz < 0)) return null;
    return [((view.p00 * vx / -vz) + 1) / 2 * view.width, (1 - (view.p11 * vy / -vz)) / 2 * view.height];
  };
  const rel = sky.velocityKms;
  let best = null;
  for (let i = 0; i < sky.count; i++) {
    const bx = sky.offsetPc[3 * i] + rel[3 * i] * KMS_TO_PC_PER_YEAR * years;
    const by = sky.offsetPc[3 * i + 1] + rel[3 * i + 1] * KMS_TO_PC_PER_YEAR * years;
    const bz = sky.offsetPc[3 * i + 2] + rel[3 * i + 2] * KMS_TO_PC_PER_YEAR * years;
    const near = toScreen(bx, by, bz);
    if (!near || Math.hypot(near[0] - pointer.x, near[1] - pointer.y) > reach + PREFILTER_SLACK_PX) continue;
    const seen = retardedOffset([sky.offsetPc[3 * i], sky.offsetPc[3 * i + 1], sky.offsetPc[3 * i + 2]],
      [rel[3 * i], rel[3 * i + 1], rel[3 * i + 2]], sky.observerVelocityKms, years);
    const magnitude = apparentMagnitude(sky.magnitude[i], Math.max(seen.distance, 1e-6));
    if (magnitude > limit) continue;
    const spot = toScreen(seen.apparent[0], seen.apparent[1], seen.apparent[2]);
    if (!spot) continue;
    const away = Math.hypot(spot[0] - pointer.x, spot[1] - pointer.y);
    if (away > reach) continue;
    if (!best || away < best.screenPx - 0.5 || (Math.abs(away - best.screenPx) <= 0.5 && magnitude < best.magnitude)) {
      best = {
        place: i, index: sky.index[i], kind: sky.kind[i], distancePc: seen.distance, tauYears: seen.tau, magnitude,
        temperatureK: sky.temperatureK[i], screenPx: away,
      };
    }
  }
  return best;
}

// Where the brightest stars of `sky` (at most `maxCount`, by apparent magnitude) are on the screen now, `[{ place, index, x, y, magnitude }]`: the stars that
// are in front of the camera, not dimmer than `limit` and inside the canvas. For the browser tests, and for anything that wants a star to point at.
export function skyStarsOnScreen(sky, years, view, options) {
  const o = options || {};
  const limit = Number.isFinite(o.limit) ? o.limit : Infinity;
  const m = view.matrix;
  const found = [];
  for (let i = 0; i < sky.count; i++) {
    const seen = retardedOffset([sky.offsetPc[3 * i], sky.offsetPc[3 * i + 1], sky.offsetPc[3 * i + 2]],
      [sky.velocityKms[3 * i], sky.velocityKms[3 * i + 1], sky.velocityKms[3 * i + 2]], sky.observerVelocityKms, years);
    const magnitude = apparentMagnitude(sky.magnitude[i], Math.max(seen.distance, 1e-6));
    if (magnitude > limit) continue;
    const [x, y, z] = seen.apparent;
    const vz = m[6] * x + m[7] * y + m[8] * z;
    if (!(vz < 0)) continue;
    const sx = ((view.p00 * (m[0] * x + m[1] * y + m[2] * z) / -vz) + 1) / 2 * view.width;
    const sy = (1 - (view.p11 * (m[3] * x + m[4] * y + m[5] * z) / -vz)) / 2 * view.height;
    if (sx < 0 || sy < 0 || sx > view.width || sy > view.height) continue;
    found.push({ place: i, index: sky.index[i], x: sx, y: sy, magnitude });
  }
  found.sort((a, b) => a.magnitude - b.magnitude);
  return found.slice(0, o.maxCount || 50);
}

// What is said of a pick: `charted` is `meta.charted` (or undefined), `starUrl` the page's template with `{id}` where a star's printed ID goes (or
// empty). Returns `{ title, lines: [text], fields: [[label, text]], href }`; `href` is null for a star with no page to go to.
export function describeSkyStar(pick, charted, starUrl) {
  const entry = charted && charted[String(pick.index)];
  const name = entry && entry[1];
  const letter = spectralClass(pick.temperatureK);
  const title = name || KIND_LABELS[pick.kind] || "Star";
  const lines = [];
  const bits = [];
  if (letter) bits.push("class " + letter);
  if (Number.isFinite(pick.temperatureK)) bits.push("about " + (Math.round(pick.temperatureK / 100) * 100).toLocaleString("en-US") + " K");
  if (name) bits.push(KIND_LABELS[pick.kind].toLowerCase());
  if (bits.length) lines.push(bits.join(", "));
  const ly = pick.distancePc * PC_TO_LY;
  lines.push(ly.toLocaleString("en-US", { maximumFractionDigits: ly < 10 ? 2 : 0 }) + " light-years away (" + pick.distancePc.toLocaleString("en-US", { maximumFractionDigits: pick.distancePc < 10 ? 2 : 0 }) + " pc)");
  lines.push("apparent magnitude " + pick.magnitude.toFixed(1));
  lines.push(agoText(pick.tauYears));
  const fields = [];
  if (bits.length) fields.push(["Class", bits.join(", ")]);
  fields.push(["Distance", lines[bits.length ? 1 : 0]]);
  fields.push(["Brightness", lines[bits.length ? 2 : 1]]);
  fields.push(["Light", agoText(pick.tauYears)]);
  const href = entry && entry[0] && starUrl ? starUrl.replace("{id}", encodeURIComponent(entry[0])) : null;
  return { title, lines, fields, href };
}

// The one line the hover tooltip shows.
export function hoverText(description) {
  return description.title + ": " + description.lines.join("; ");
}
