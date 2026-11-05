// html/static/skypack.js
//
// VIEW.29: reads the packed sky list the sky route sends and the cache keeps. planetgen/web/lib/skypack.py writes
// it and describes the layout; this file is its reader, so the format has one definition on each side and a test
// (tests/test_skypack.py with tests/js/skypack.test.mjs) holds the two together.
//
// A list is "PGK1", the star count (u32), the time the offsets hold at (f64 Unix seconds), the limiting magnitude
// (f32), the velocity step (f32 km/s a code), the observer's velocity (3 x f64 km/s), the JSON length (u32), the
// JSON and then the columns: offset x, y, z (f32 parsecs from the observer at that time, not seen back by the
// light time), velocity x, y, z (i16 codes, km/s relative to the observer), the dimmed absolute magnitude (i16,
// 1/512 magnitude), the temperature (u8 over log10 2..6, 0 for none), the kind (u8) and the index (u32).
// The graphics card places a star at offset + relativeVelocity * t - velocity * tau, with tau = distance / c and
// velocity = relativeVelocity + observerVelocity.

export const TEMPERATURE_LOG_RANGE = [2, 6];
export const MAGNITUDE_STEPS_PER_MAG = 512;
export const KIND_NAMES = ["charted", "bright_star", "phenomenon", "cluster"];
export const BYTES_PER_STAR = 12 + 6 + 2 + 1 + 1 + 4;

const MAGIC = "PGK1";
const HEADER_BYTES = 52;
const U8 = 255;

// One column of `count` values read with `Constructor` (a typed array) from `bytes` at `offset`; copies, so the
// offset need not be aligned.
function column(bytes, offset, count, Constructor) {
  const size = count * Constructor.BYTES_PER_ELEMENT;
  const copy = bytes.slice(offset, offset + size);
  if (copy.length !== size) throw new Error("a packed sky ends early");
  return new Constructor(copy.buffer);
}

// The three columns of an axis group joined, one star after the other (x, y, z, x, y, z, ...): what a points
// object's position attribute takes.
function interleave(planes, count, scale) {
  const out = new Float32Array(3 * count);
  for (let axis = 0; axis < 3; axis++) {
    const plane = planes[axis];
    for (let i = 0; i < count; i++) out[3 * i + axis] = plane[i] * scale;
  }
  return out;
}

// A packed sky (a Uint8Array) as `{ count, epochUnix, limit, velocityScale, observerVelocityKms, meta, offsetPc,
// velocityKms, magnitude, temperatureK, kind, index }`: the offsets and velocities interleaved Float32Arrays of
// 3 x count, `magnitude` the dimmed absolute magnitude, `temperatureK` NaN where none, `index` a Float64Array.
// Throws on bytes that are not one whole packed sky.
export function unpackSky(bytes) {
  if (bytes.length < HEADER_BYTES) throw new Error("not a packed sky");
  if (String.fromCharCode(bytes[0], bytes[1], bytes[2], bytes[3]) !== MAGIC) throw new Error("not a packed sky");
  const view = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
  const count = view.getUint32(4, true);
  const epochUnix = view.getFloat64(8, true);
  const limit = view.getFloat32(16, true);
  const velocityScale = view.getFloat32(20, true);
  const observerVelocityKms = [view.getFloat64(24, true), view.getFloat64(32, true), view.getFloat64(40, true)];
  const metaLength = view.getUint32(48, true);
  if (bytes.length !== HEADER_BYTES + metaLength + count * BYTES_PER_STAR) throw new Error("not a packed sky");
  let at = HEADER_BYTES;
  const meta = JSON.parse(new TextDecoder("utf-8").decode(bytes.subarray(at, at + metaLength)));
  at += metaLength;
  const planes = (Constructor) => {
    const out = [];
    for (let axis = 0; axis < 3; axis++) {
      out.push(column(bytes, at, count, Constructor));
      at += count * Constructor.BYTES_PER_ELEMENT;
    }
    return out;
  };
  const offsetPc = interleave(planes(Float32Array), count, 1);
  const velocityKms = interleave(planes(Int16Array), count, velocityScale);
  const codes = column(bytes, at, count, Int16Array);
  at += 2 * count;
  const magnitude = new Float32Array(count);
  for (let i = 0; i < count; i++) magnitude[i] = codes[i] / MAGNITUDE_STEPS_PER_MAG;
  const temperatureK = new Float32Array(count);
  for (let i = 0; i < count; i++) {
    const code = bytes[at + i];
    temperatureK[i] = code === 0 ? NaN
      : Math.pow(10, TEMPERATURE_LOG_RANGE[0] + (code - 1) / (U8 - 1) * (TEMPERATURE_LOG_RANGE[1] - TEMPERATURE_LOG_RANGE[0]));
  }
  at += count;
  const kind = bytes.slice(at, at + count);
  at += count;
  const index = Float64Array.from(column(bytes, at, count, Uint32Array));
  return { count, epochUnix, limit, velocityScale, observerVelocityKms, meta, offsetPc, velocityKms, magnitude, temperatureK, kind, index };
}
