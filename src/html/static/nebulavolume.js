// html/static/nebulavolume.js
//
// MAP.365: the nebulae of the Galaxy Map drawn from their shape, the choosing half (nebulavolumelayer.js draws them). The server
// samples the galaxy's gas and star-forming density for each node of the map (`GET /galaxy/nebula-volume?key=<node key>`,
// planetgen/galaxy/nebula_volume.py): points on an absolute lattice, each with a strength and the family of the nebula that owns it.
// A nebula is therefore a cloud of soft splats whose outline is the density's own: irregular, a dense core fading unevenly, coarse
// when its node is coarse and with more filaments as the nodes get finer. Two nodes draw the points that lie in them, so a nebula
// that crosses a border is drawn without a seam, and where nebulae overlap on the screen their colors add up under a soft cap
// (`cappedSum`), so a crowd of them never washes to white.
//
// No DOM and no three.js here: tests/js/nebulavolume.test.mjs holds it to this description.

// What the server packs (planetgen.galaxy.nebula_volume.pack).
export const MAGIC = "NVOL";
export const PACK_VERSION = 1;
const HEADER_BYTES = 14;

// A node's volume is asked for once a nebula in it is this many pixels across (radius), and at most this many are held at once, the
// biggest on the screen first; MAX_REQUESTS are in flight at a time.
export const MIN_RADIUS_PX = 14;
export const FULL_RADIUS_PX = 28;
export const MAX_LOADED = 24;
export const MAX_REQUESTS = 3;
// A splat is this many lattice steps wide (they overlap, so the sum is smooth) and no smaller than MIN_SPLAT_PX on the screen.
export const SPLAT_STEPS = 2.2;
export const MIN_SPLAT_PX = 2;
export const MAX_SPLAT_PX = 420;
// How much one point at full strength adds (a nebula is about twelve steps across), and the most the nebulae add up to on the map
// (the largest value any color channel of the layer reaches, 0 to 1): a soft cap, so the sum bends toward it and never passes it.
export const POINT_ALPHA = 0.045;
export const CAP = 0.5;

function ease(x) {
  const t = Math.min(1, Math.max(0, x));
  return t * t * (3 - 2 * t);
}

// The points of a packed volume: `{ level, spacingPc, count, ix, iy, iz, strength, family }` (the lattice steps are Int32Arrays; a point is
// at `i * spacingPc` parsecs from the galaxy's origin). Throws on bytes that are not a volume.
export function decodeVolume(buffer) {
  const bytes = buffer instanceof Uint8Array ? buffer : new Uint8Array(buffer);
  if (bytes.length < HEADER_BYTES || String.fromCharCode(bytes[0], bytes[1], bytes[2], bytes[3]) !== MAGIC) {
    throw new Error("not a nebula volume");
  }
  if (bytes[4] !== PACK_VERSION) throw new Error("nebula volume version " + bytes[4]);
  const view = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
  const count = view.getUint32(6, true);
  const spacingPc = view.getFloat32(10, true);
  if (bytes.length !== HEADER_BYTES + count * 14) throw new Error("nebula volume is the wrong length");
  const column = (offset) => {
    const out = new Int32Array(count);
    for (let i = 0; i < count; i++) out[i] = view.getInt32(offset + 4 * i, true);
    return out;
  };
  const start = HEADER_BYTES;
  return {
    level: bytes[5], spacingPc, count,
    ix: column(start), iy: column(start + 4 * count), iz: column(start + 8 * count),
    strength: bytes.slice(start + 12 * count, start + 13 * count),
    family: bytes.slice(start + 13 * count, start + 14 * count),
  };
}

// The positions (a Float32Array of x, y, z) of a volume's points relative to `frame` ([x, y, z] parsecs), so a float32 holds them.
export function positionsOf(volume, frame) {
  const out = new Float32Array(3 * volume.count);
  for (let i = 0; i < volume.count; i++) {
    out[3 * i] = volume.ix[i] * volume.spacingPc - frame[0];
    out[3 * i + 1] = volume.iy[i] * volume.spacingPc - frame[1];
    out[3 * i + 2] = volume.iz[i] * volume.spacingPc - frame[2];
  }
  return out;
}

// The middle of a volume's points, parsecs.
export function centerOf(volume) {
  const sum = [0, 0, 0];
  for (let i = 0; i < volume.count; i++) {
    sum[0] += volume.ix[i];
    sum[1] += volume.iy[i];
    sum[2] += volume.iz[i];
  }
  const n = Math.max(1, volume.count);
  return sum.map((v) => (v / n) * volume.spacingPc);
}

// The width in pixels of one splat at `distance` parsecs, for a camera with focal length `focalPx`.
export function splatPx(spacingPc, distance, focalPx) {
  const px = (SPLAT_STEPS * spacingPc * focalPx) / Math.max(distance, 1e-6);
  return Math.min(MAX_SPLAT_PX, Math.max(MIN_SPLAT_PX, px));
}

// The soft cap: `m` (the largest channel of a sum) bent so it approaches `cap` and never reaches it; near zero it is `m`.
export function softCap(m, cap) {
  const limit = cap === undefined ? CAP : cap;
  return limit * (1 - Math.exp(-Math.max(0, m) / limit));
}

// The color the nebulae over one pixel make: `terms` is a list of [r, g, b] (each already its color times its strength there),
// added, then bent so the brightest channel stays under `cap`, keeping the hue (every channel is scaled by the same amount).
export function cappedSum(terms, cap) {
  const sum = [0, 0, 0];
  (terms || []).forEach((term) => {
    sum[0] += term[0];
    sum[1] += term[1];
    sum[2] += term[2];
  });
  const m = Math.max(sum[0], sum[1], sum[2]);
  if (!(m > 0)) return [0, 0, 0];
  const scale = softCap(m, cap) / m;
  return sum.map((v) => v * scale);
}

// MAP.215: a star, a glow or a region seen through a nebula takes its hue. `nebulaSum` is the nebulae's raw sum at the pixel (the
// layer's target before the cap): the hue is the sum's direction (its brightest channel scaled to 1) and the strength grows with its size
// and levels off at TINT_MAX (a subtle tint, denser cores tint more); a little light is lost as well (TINT_DIM), so a dense cloud
// dims what is behind it. `tintOf` is the formula; NEBULA_TINT_GLSL the same in the shaders (`nebulaTinted`).
export const TINT_MAX = 0.6;
export const TINT_SCALE = 0.25;
export const TINT_DIM = 0.2;
export function tintOf(sum) {
  const m = Math.max(sum[0], sum[1], sum[2]);
  if (!(m > 0)) return { strength: 0, hue: [1, 1, 1] };
  return { strength: TINT_MAX * (1 - Math.exp(-m / TINT_SCALE)), hue: sum.map((v) => v / m) };
}
// The color `colour` ([r, g, b]) takes behind a pixel whose nebulae add up to `sum`.
export function tinted(colour, sum) {
  const { strength, hue } = tintOf(sum);
  return colour.map((c, i) => (c * (1 - strength) + c * hue[i] * strength) * (1 - TINT_DIM * strength));
}
export const NEBULA_TINT_GLSL = [
  "uniform sampler2D nebulaSum;",
  "uniform vec2 drawSize;",
  "uniform float nebulaTint;",
  "vec3 nebulaTinted(vec3 c) {",
  "  if (nebulaTint <= 0.0) return c;",
  "  vec3 s = texture2D(nebulaSum, gl_FragCoord.xy / drawSize).rgb;",
  "  float m = max(s.r, max(s.g, s.b));",
  "  if (m <= 0.0) return c;",
  `  float t = ${TINT_MAX.toFixed(3)} * (1.0 - exp(-m / ${TINT_SCALE.toFixed(3)})) * nebulaTint;`,
  "  vec3 hue = s / m;",
  `  return mix(c, c * hue, t) * (1.0 - ${TINT_DIM.toFixed(3)} * t);`,
  "}",
].join("\n");

// How strongly a node's volume shows: from the largest nebula in it on the screen, `radiusPx`, and the share the cut gives the node.
export function volumeAlpha(radiusPx, share) {
  return ease((radiusPx - MIN_RADIUS_PX) / (FULL_RADIUS_PX - MIN_RADIUS_PX)) * Math.min(1, Math.max(0, share));
}

// The nebulae a drawn node lists: its groups (`nebula_regions`, a center and a radius each) and its single clouds. Returns
// `{ key, block, reaches }` for a node that holds any, `reaches` a list of `{ x, y, z, radiusPc }` (the groups' and the clouds' spheres).
export function sourceOf(entry) {
  const node = entry && entry.node;
  if (!node) return null;
  const reaches = [];
  (node.nebula_regions || []).forEach((region) => {
    if (region && region.center && region.radius_pc > 0) {
      reaches.push({ x: region.center[0], y: region.center[1], z: region.center[2], radiusPc: region.radius_pc });
    }
  });
  (node.clouds || []).forEach((cloud) => {
    if (cloud && cloud.radius_pc > 0) reaches.push({ x: cloud.x, y: cloud.y, z: cloud.z, radiusPc: cloud.radius_pc });
  });
  return reaches.length ? { key: entry.key, block: entry.block, reaches } : null;
}

// MAP.142: a source for one nebula drawn alone (the Sector Map's): `GET /galaxy/nebula/<id>/volume` holds the whole of it, so it is its own key
// (`nebula:<id>`), has no node block to take a share from and shows at `weight` (0 to 1: the sector's fade, 0 for a hidden kind) of its strength.
export function singleSource(id, sphere, weight) {
  return {
    key: "nebula:" + id, single: true, block: null, weight: Math.min(1, Math.max(0, weight === undefined ? 1 : weight)),
    reaches: [{ x: sphere.x, y: sphere.y, z: sphere.z, radiusPc: sphere.radiusPc }],
  };
}

// The id of a single source's key (`nebula:<id>`), or null for a node key.
export function singleIdOf(key) {
  return typeof key === "string" && key.indexOf("nebula:") === 0 ? key.slice("nebula:".length) : null;
}

// The pixel radius of the biggest nebula of `source` as the camera at `camera` ([x, y, z] parsecs) with focal length `focalPx` sees it.
export function biggestPx(source, camera, focalPx) {
  let best = 0;
  source.reaches.forEach((sphere) => {
    const distance = Math.max(1e-6, Math.hypot(sphere.x - camera[0], sphere.y - camera[1], sphere.z - camera[2]));
    best = Math.max(best, (sphere.radiusPc * focalPx) / distance);
  });
  return best;
}

// Which sources to hold: those with a nebula at least MIN_RADIUS_PX across and in front of the camera, the biggest first, at most
// MAX_LOADED; `[{ source, radiusPx }]`.
export function wanted(sources, camera, focalPx, onScreen) {
  const found = [];
  sources.forEach((source) => {
    const radiusPx = biggestPx(source, camera, focalPx);
    if (radiusPx >= MIN_RADIUS_PX && (!onScreen || onScreen(source))) found.push({ source, radiusPx });
  });
  found.sort((a, b) => (b.radiusPx - a.radiusPx) || (a.source.key < b.source.key ? -1 : 1));
  return found.slice(0, MAX_LOADED);
}

// A fetcher of volumes that asks the server once for each node key and hands every later caller the same state, at most MAX_REQUESTS
// at a time. `fetchBytes(key) -> Promise<ArrayBuffer>`; `onChange()` is called whenever one arrives or fails. A failed key is tried
// again after RETRY_MS.
export const RETRY_MS = 15000;
export function createVolumeStore(fetchBytes, onChange, clock) {
  const now = clock || (() => Date.now());
  const states = new Map();
  const queue = [];
  let inFlight = 0;
  function pump() {
    while (inFlight < MAX_REQUESTS && queue.length) {
      const key = queue.shift();
      const state = states.get(key);
      if (!state || state.status !== "queued") continue;
      state.status = "loading";
      inFlight += 1;
      Promise.resolve().then(() => fetchBytes(key)).then((bytes) => {
        state.volume = decodeVolume(bytes);
        state.status = "ready";
      }).catch(() => {
        state.status = "failed";
        state.failedAt = now();
      }).then(() => {
        inFlight -= 1;
        if (onChange) onChange(key);
        pump();
      });
    }
  }
  return {
    // Asks for `key`'s volume if it is not held; returns its state `{ status: "queued"|"loading"|"ready"|"failed", volume? }`.
    request(key) {
      let state = states.get(key);
      if (state && state.status === "failed" && now() - state.failedAt >= RETRY_MS) state = null;
      if (!state) {
        state = { status: "queued", volume: null };
        states.set(key, state);
        queue.push(key);
        pump();
      }
      return state;
    },
    get: (key) => states.get(key) || null,
    // Forgets the volumes of nodes not in `keep` (a Set), so the memory follows the view.
    release(keep) {
      Array.from(states.keys()).forEach((key) => {
        const state = states.get(key);
        if (!keep.has(key) && state.status !== "loading") states.delete(key);
      });
    },
    size: () => states.size,
  };
}
