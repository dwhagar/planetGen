// html/static/spin.js
//
// The slow turn of a rendered body on its axis (MAP.328): one function for every view that draws a body as a sphere
// (the System Map's spheres, the 3D system view, the rendered example on an object page, the stars of the Sector Map),
// so they agree on how fast a body turns, which way, and when it stands still. The galactic views (the Galaxy Map,
// regions, sector blocks as points and glows) draw no turning body and do not use it.
//
// A body's look record (planetgen/web/maps/bodylook.py) carries what is stored for it: `spinH` (its rotation period,
// hours), `axis` (the unit spin axis in the system's frame, right-handed, so a retrograde body's axis points against
// its orbit's normal) and `tiltDeg` (the axial tilt, 0 to 180 degrees; past 90 it turns backwards). A body with no
// `spinH` has no stored period and is left still: nothing is made up for it.
//
// The rate is the stored period run through ONE shared law, for every body, so the order of speeds stays true:
// the seconds per turn on screen are `FASTEST_SECONDS_PER_TURN * (SLOWEST / FASTEST) ** t`, `t` being where the
// period falls on a log scale from `FASTEST_PERIOD_HOURS` (t = 0) to `SLOWEST_PERIOD_HOURS` (t = 1), clamped to the
// ends. A day of 5 hours shows a turn in about 31 s, an Earth day in about 70 s, the Sun's 25 days in about 180 s,
// a 47-day rotation in about 250 s, Venus's 243 days in about 480 s, and nothing turns faster than one turn in 20 s
// (it would strobe) or slower than one in 10 minutes (it would look stuck). A straight proportion cannot do both: real
// periods run from hours to hundreds of days, 100 times wider than the range a screen can show a turn in (a deliberate
// tweak, recorded as stand-in SP1 at tweak level L4 in docs/design/research/gap-filling-plan.md section 6.9). Every
// body turns far faster on screen than it does: the pages that show a picture say so (`DISPLAY_NOTE`).
//
// A body that is tidally locked (`locked`: its rotation period is its orbital period) does not turn by that law. It
// keeps one face toward what it orbits, so in a view that moves it along its orbit it turns once per orbit, which
// is the same as the turn of nothing at all as the star sees it. `turn` therefore takes the direction to the body it
// orbits and lays the sub-parent point (longitude 0) on it, wherever the orbit has taken the body; the pictures that
// show a body from the star's side (the object page) hold it still with its day side toward the star.
//
// The Spin switch (the System Map's Menu) and the reader's reduced-motion setting turn every body still; the choice
// is kept in localStorage, which can be missing or throw, so every access is guarded and the default is on.
//
// No DOM at import and no three.js import (the caller's own THREE is passed to `turn`), so tests/js/spin.test.mjs
// holds the numbers to the description above.

export const FASTEST_SECONDS_PER_TURN = 20;
export const SLOWEST_SECONDS_PER_TURN = 600;
export const FASTEST_PERIOD_HOURS = 2;
export const SLOWEST_PERIOD_HOURS = 10000;

export const DISPLAY_NOTE = "Shown turning much faster than the real rotation.";
export const LOCKED_NOTE = "Tidally locked: one face stays toward what it orbits.";

const KEY = "planetgen.spin";

// Seconds one turn takes on screen for a body of `periodHours`, or null for a body with no usable period.
export function secondsPerTurn(periodHours) {
  if (!(periodHours > 0) || !Number.isFinite(periodHours)) return null;
  const span = Math.log(SLOWEST_PERIOD_HOURS / FASTEST_PERIOD_HOURS);
  const t = Math.min(1, Math.max(0, Math.log(periodHours / FASTEST_PERIOD_HOURS) / span));
  return FASTEST_SECONDS_PER_TURN * Math.pow(SLOWEST_SECONDS_PER_TURN / FASTEST_SECONDS_PER_TURN, t);
}

// Radians per second about the body's axis (positive: right-handed about `axisOf(look)`), 0 for a body left still.
export function rateOf(look) {
  const seconds = look ? secondsPerTurn(look.spinH) : null;
  return seconds === null ? 0 : (2 * Math.PI) / seconds;
}

// Where a body starts, so a scene's bodies do not all show the same face at once: a fixed angle from its seed.
export function phaseOf(look) {
  return look && look.seed !== undefined ? ((Number(look.seed) >>> 0) % 6283) / 1000 : 0;
}

// Whether the body is tidally locked: the look record says its rotation period is its orbital period.
export function isLocked(look) {
  return !!(look && look.locked);
}

// The angle (about the axis, in the body's own frame) that lays the body's longitude 0 (the sphere's local +X, the middle
// of its equirectangular texture) along `direction` ([x, y, z] in
// the frame the body's axis is given in), the way a locked body keeps its face toward its parent; `axis` is the body's
// axis (`axisOf`). Null when the direction runs along the axis (no face to turn).
export function faceAngle(axis, direction) {
  const a = unit(axis);
  const d = unit(direction);
  if (!a || !d) return null;
  // The direction in the frame where the axis is +Y (the shortest turn, the one `turn` uses).
  const dot = a[1];
  let v;
  if (dot <= -1 + 1e-12) {
    v = [-d[0], -d[1], d[2]]; // three.js turns up to straight down about the Z axis
  } else {
    // Rodrigues: the rotation taking `a` to +Y is about a x Y.
    const cross = [-a[2], 0, a[0]];
    const f = 1 / (1 + dot);
    const cd = [cross[1] * d[2] - cross[2] * d[1], cross[2] * d[0] - cross[0] * d[2], cross[0] * d[1] - cross[1] * d[0]];
    const cc = [cross[1] * cd[2] - cross[2] * cd[1], cross[2] * cd[0] - cross[0] * cd[2], cross[0] * cd[1] - cross[1] * cd[0]];
    v = [d[0] + cd[0] + cc[0] * f, d[1] + cd[1] + cc[1] * f, d[2] + cd[2] + cc[2] * f];
  }
  if (Math.hypot(v[0], v[2]) < 1e-9) return null;
  return Math.atan2(-v[2], v[0]);
}

// How far the body has turned after `seconds` of spin time (see `tick`).
export function angleAt(look, seconds) {
  return phaseOf(look) + rateOf(look) * seconds;
}

function unit(vector) {
  const length = Math.hypot(vector[0], vector[1], vector[2]);
  return length > 0 ? [vector[0] / length, vector[1] / length, vector[2] / length] : null;
}

// The stored spin axis in the system's frame (the 3D view draws in it), else straight up.
export function axisOf(look) {
  const axis = look && Array.isArray(look.axis) && look.axis.length === 3 ? unit(look.axis) : null;
  return axis || [0, 1, 0];
}

// The axis a flat picture of the body shows: up, leaned over by the axial tilt in the picture's plane, so a body
// tilted past 90 degrees turns the other way.
export function pictureAxisOf(look) {
  const tilt = look && Number.isFinite(look.tiltDeg) ? (look.tiltDeg * Math.PI) / 180 : 0;
  return [Math.sin(tilt), Math.cos(tilt), 0];
}

// Where on a body a point of its flat picture lies: `(x, up, z)` is the unit vector from the sphere's middle to the point
// (x to the right, up, z toward the viewer), the body's axis leaning over by its axial tilt (`pictureAxisOf`). Returns
// `[longitude, latitude]` in radians in the body's own frame, longitude 0 facing the viewer when the body has not turned
// (`angleAt` = 0) and growing eastward, the way a prograde body turns, so after the body has turned `angle` the picture
// shows its longitude `longitude - angle` there.
export function bodyPoint(look, x, up, z) {
  const tilt = look && Number.isFinite(look.tiltDeg) ? (look.tiltDeg * Math.PI) / 180 : 0;
  const c = Math.cos(tilt);
  const s = Math.sin(tilt);
  const bodyX = x * c - up * s;
  const bodyUp = x * s + up * c;
  return [Math.atan2(bodyX, z), Math.asin(Math.max(-1, Math.min(1, bodyUp)))];
}

// --- Whether bodies turn ----------------------------------------------------

const memory = { on: true };

export function switchedOn() {
  try {
    return globalThis.localStorage.getItem(KEY) !== "off";
  } catch (e) {
    return memory.on;
  }
}

// Turns the Spin switch on or off and keeps the choice (for this page only where storage is unavailable).
export function setSwitch(on) {
  memory.on = !!on;
  try {
    if (on) globalThis.localStorage.removeItem(KEY);
    else globalThis.localStorage.setItem(KEY, "off");
  } catch (e) {
    // This page only.
  }
  return !!on;
}

export function reducedMotion() {
  return typeof globalThis.matchMedia === "function" && !!globalThis.matchMedia("(prefers-reduced-motion: reduce)").matches;
}

// True when bodies turn: the switch is on and the reader has not asked for reduced motion.
export function active() {
  return switchedOn() && !reducedMotion();
}

// --- Spin time ----------------------------------------------------------------

// Spin time runs only while bodies turn, so switching Spin off stops them where they are and on takes them on from there.
const clock = { seconds: 0, lastMs: null };

// Advances spin time to `nowMs` (a frame's timestamp, milliseconds) and returns it in seconds. A long gap (a hidden tab)
// counts for at most a second, so bodies do not jump when the tab returns.
export function tick(nowMs) {
  if (clock.lastMs !== null && active()) clock.seconds += Math.min(1, Math.max(0, (nowMs - clock.lastMs) / 1000));
  clock.lastMs = nowMs;
  return clock.seconds;
}

export function seconds() {
  return clock.seconds;
}

// --- Drawing ------------------------------------------------------------------

const scratch = {};

// Turns `object` (a three.js sphere, poles on its local Y) about `axis` ([x, y, z], default `axisOf(look)`) to where
// `look` says it is after `spinSeconds`: its Y is laid on the axis, then it is turned about it. A locked body, when `toward`
// (the direction to the body it orbits, in the frame of the axis) is given, is turned to face it instead.
export function turn(THREE, object, look, spinSeconds, axis, toward) {
  if (!scratch.up || scratch.THREE !== THREE) {
    scratch.THREE = THREE;
    scratch.up = new THREE.Vector3(0, 1, 0);
    scratch.axis = new THREE.Vector3();
    scratch.lean = new THREE.Quaternion();
    scratch.spin = new THREE.Quaternion();
  }
  const a = axis || axisOf(look);
  scratch.axis.set(a[0], a[1], a[2]).normalize();
  scratch.lean.setFromUnitVectors(scratch.up, scratch.axis);
  const face = isLocked(look) && toward ? faceAngle(a, toward) : null;
  scratch.spin.setFromAxisAngle(scratch.up, face === null ? angleAt(look, spinSeconds) : face);
  object.quaternion.copy(scratch.lean).multiply(scratch.spin);
  return object;
}
