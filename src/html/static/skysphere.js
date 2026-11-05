// html/static/skysphere.js
//
// VIEW.33 (S15 of the System Sky plan, reports/system-sky-from-the-system-2026-11-04.md): the stars of the sky a system sees, drawn as one
// points object centred on the camera. The list is a packed sky (skypack.js, VIEW.29): offsets from the observer in parsecs at the list's
// own time, each star's velocity relative to the observer, its dimmed absolute magnitude and its temperature. The vertex shader does the
// rest every frame, so the orbit clock runs the sky:
//
//     base     = offset + relative velocity * t                    where the star is now, relative to the observer
//     apparent = base - velocity * tau,  tau = |apparent| / c      where it was when its light left it (two passes, three for a star
//                                                                  faster than FAST_KMS), velocity = relative + the observer's
//     direction = apparent / |apparent|,  distance = |apparent|
//     m        = M + 5 log10(distance / 10 pc)                     skylight.js: opacity, size and colour from m, the fade at the cut
//
// The points sit on a sphere round the camera (the object's position follows the camera and its depth is not tested), so the sky is at
// infinity: a move of the camera inside the system changes no direction. The offsets are floats on the card (the server works in doubles
// and sends offsets from the observer), which keeps a star at 1 kpc to a hundredth of an arcsecond.
//
// The arithmetic that is also in the shader is in JavaScript here (`retardedOffset`, `lightPasses`) so a node test and a browser test hold
// the two to one answer. `THREE` is passed in, so plain node can test the rest.

const VERSION_QUERY = new URL(import.meta.url).search;
const { SKY_LIGHT_GLSL, skyLimit, starColour } = await import(`./skylight.js${VERSION_QUERY}`);
const { apparentMagnitude } = await import(`./starmagnitude.js${VERSION_QUERY}`);

export const LIGHT_SPEED_KMS = 299792.458;
export const PARSEC_KM = 3.0856775814913673e13;
export const YEAR_SECONDS = 31557600;
export const SECONDS_PER_YEAR = YEAR_SECONDS;
// 1 km/s in parsecs a year, and the speed of light in parsecs a year.
export const KMS_TO_PC_PER_YEAR = YEAR_SECONDS / PARSEC_KM;
export const LIGHT_PC_PER_YEAR = (LIGHT_SPEED_KMS * YEAR_SECONDS) / PARSEC_KM;
// A star faster than this (relative to the observer) takes a third light-time pass, as the plan says (hypervelocity stars).
export const FAST_KMS = 1000;
// The most stars drawn: a phone draws fewer (the plan's 15,000), the rest of the list is kept for the cut.
export const MAX_STARS = 60000;
export const MAX_STARS_PHONE = 15000;
// How far from the camera the sky is drawn, as a share of the way between its near and far planes (the depth is not tested).
export const SKY_RADIUS_SHARE = 0.5;

// How many light-time passes a star moving `speedKms` (relative to the observer) takes.
export function lightPasses(speedKms) {
  return speedKms > FAST_KMS ? 3 : 2;
}

// The years of the sky's own time after the list's epoch, from the orbit clock: `clockYears` is the clock's years after the scene's epoch
// (`sceneEpochUnix`), the list holds at `sky.epochUnix`.
export function skyYears(sky, sceneEpochUnix, clockYears) {
  return clockYears + (sceneEpochUnix - sky.epochUnix) / YEAR_SECONDS;
}

// Where a star appears: its offset (pc) at the list's epoch, its velocity relative to the observer and the observer's (km/s), `years` after the
// epoch. The same arithmetic as the vertex shader, as `{ apparent: [x, y, z], distance, tau }` (tau in years).
export function retardedOffset(offset, relativeKms, observerKms, years) {
  const rel = relativeKms.map((v) => v * KMS_TO_PC_PER_YEAR);
  const star = relativeKms.map((v, i) => (v + observerKms[i]) * KMS_TO_PC_PER_YEAR);
  const base = offset.map((o, i) => o + rel[i] * years);
  const speed = Math.hypot(relativeKms[0], relativeKms[1], relativeKms[2]);
  let apparent = base;
  let tau = 0;
  for (let pass = 0; pass < lightPasses(speed); pass++) {
    apparent = base.map((b, i) => b - star[i] * tau);
    tau = Math.hypot(apparent[0], apparent[1], apparent[2]) / LIGHT_PC_PER_YEAR;
  }
  apparent = base.map((b, i) => b - star[i] * tau);
  return { apparent, distance: Math.hypot(apparent[0], apparent[1], apparent[2]), tau };
}

// The apparent magnitude of star `i` of `sky` where it stands at the list's epoch (no light time: the order the cap keeps stars in).
function epochMagnitude(sky, i) {
  const d = Math.hypot(sky.offsetPc[3 * i], sky.offsetPc[3 * i + 1], sky.offsetPc[3 * i + 2]);
  return apparentMagnitude(sky.magnitude[i], d);
}

// The `maxCount` brightest stars of `sky` (by apparent magnitude at the epoch), in that order, as a sky of the same shape; `sky` itself when it
// holds no more. A phone draws fewer than a desktop, so the dimmest go first, the way the limit would take them.
export function keepBrightest(sky, maxCount) {
  if (sky.count <= maxCount) return sky;
  const order = Array.from({ length: sky.count }, (_, i) => i);
  const m = Float32Array.from(order, (i) => epochMagnitude(sky, i));
  order.sort((a, b) => (m[a] - m[b]) || (a - b));
  const keep = order.slice(0, maxCount);
  const pick = (array, width) => {
    const out = new array.constructor(keep.length * width);
    keep.forEach((from, to) => { for (let k = 0; k < width; k++) out[to * width + k] = array[from * width + k]; });
    return out;
  };
  return {
    ...sky, count: keep.length, offsetPc: pick(sky.offsetPc, 3), velocityKms: pick(sky.velocityKms, 3), magnitude: pick(sky.magnitude, 1),
    temperatureK: pick(sky.temperatureK, 1), kind: pick(sky.kind, 1), index: pick(sky.index, 1),
  };
}

// The colour (linear, a hue) of every star of `sky`, 3 floats each: its blackbody, white for a star with no temperature.
export function starColours(sky) {
  const out = new Float32Array(3 * sky.count);
  for (let i = 0; i < sky.count; i++) {
    const rgb = starColour(sky.temperatureK[i]);
    out[3 * i] = rgb[0];
    out[3 * i + 1] = rgb[1];
    out[3 * i + 2] = rgb[2];
  }
  return out;
}

export const VERTEX_SHADER = [
  "attribute vec3 relativeKms;",
  "attribute float absMag;",
  "attribute vec3 hue;",
  "uniform float years;",
  "uniform vec3 observerKms;",
  "uniform float limit;",
  "uniform float pixelRatio;",
  "uniform float skyRadius;",
  "varying vec3 vColour;",
  "varying float vOpacity;",
  SKY_LIGHT_GLSL,
  "const float KMS = " + KMS_TO_PC_PER_YEAR.toExponential(9) + ";",
  "const float LIGHT = " + LIGHT_PC_PER_YEAR.toFixed(9) + ";",
  "void main() {",
  "  vec3 rel = relativeKms * KMS;",
  "  vec3 star = (relativeKms + observerKms) * KMS;",
  "  vec3 base = position + rel * years;",
  "  vec3 apparent = base;",
  "  float tau = 0.0;",
  "  for (int i = 0; i < 2; i++) {",
  "    apparent = base - star * tau;",
  "    tau = length(apparent) / LIGHT;",
  "  }",
  "  if (length(relativeKms) > " + FAST_KMS.toFixed(1) + ") {",
  "    apparent = base - star * tau;",
  "    tau = length(apparent) / LIGHT;",
  "  }",
  "  apparent = base - star * tau;",
  "  float distancePc = max(length(apparent), 1e-6);",
  "  float m = apparentMagnitude(absMag, distancePc);",
  "  float opacity = skyOpacity(m, limit);",
  "  vOpacity = opacity;",
  "  vColour = skyTint(hue, m);",
  "  vec3 direction = apparent / distancePc;",
  "  gl_Position = projectionMatrix * modelViewMatrix * vec4(direction * skyRadius, 1.0);",
  "  gl_PointSize = opacity > 0.0 ? skySizePx(m) * pixelRatio : 0.0;",
  "}",
].join("\n");

// A Gaussian: the point is SIZE_SIGMAS (4) sigma across, so sigma is a quarter of its width and exp(-8 r^2) is that for r in 0..0.5.
export const FRAGMENT_SHADER = [
  "varying vec3 vColour;",
  "varying float vOpacity;",
  "void main() {",
  "  vec2 c = gl_PointCoord - vec2(0.5);",
  "  float g = exp(-8.0 * dot(c, c)) * (1.0 - smoothstep(0.35, 0.5, length(c)));",
  "  gl_FragColor = vec4(vColour * (vOpacity * g), 1.0);",
  "}",
].join("\n");

// The sky of `sky` (skypack.js `unpackSky`) as a points object. `options`: `phone` (a lower cap and a lower limit), `limit` (the limiting
// magnitude, default skylight.skyLimit), `pixelRatio` (default 1), `maxStars` (default MAX_STARS or MAX_STARS_PHONE).
// Returns `{ object, count, setLimit(limit), limit(), update(camera, { years, pixelRatio }), dispose() }`:
// `object` goes into the scene (drawn first, depth not tested), `update` follows the camera and sets the sky's time (`years` after the list's
// epoch: `skyYears`). Rotate `object.quaternion` to turn the sky's galactic axes into the scene's.
export function createSkySphere(THREE, sky, options) {
  const o = options || {};
  const cap = o.maxStars || (o.phone ? MAX_STARS_PHONE : MAX_STARS);
  const shown = keepBrightest(sky, cap);
  const geometry = new THREE.BufferGeometry();
  geometry.setAttribute("position", new THREE.BufferAttribute(shown.offsetPc, 3));
  geometry.setAttribute("relativeKms", new THREE.BufferAttribute(shown.velocityKms, 3));
  geometry.setAttribute("absMag", new THREE.BufferAttribute(shown.magnitude, 1));
  geometry.setAttribute("hue", new THREE.BufferAttribute(starColours(shown), 3));
  let limit = Number.isFinite(o.limit) ? o.limit : skyLimit(!!o.phone);
  const uniforms = {
    years: { value: 0 },
    observerKms: { value: new THREE.Vector3(...shown.observerVelocityKms) },
    limit: { value: limit },
    pixelRatio: { value: o.pixelRatio || 1 },
    skyRadius: { value: 1 },
  };
  const material = new THREE.ShaderMaterial({
    uniforms, vertexShader: VERTEX_SHADER, fragmentShader: FRAGMENT_SHADER,
    transparent: true, depthTest: false, depthWrite: false, blending: THREE.AdditiveBlending,
  });
  const object = new THREE.Points(geometry, material);
  object.frustumCulled = false;
  object.renderOrder = -10;
  object.name = "sky-sphere";
  return {
    object, count: shown.count, kept: shown,
    setLimit(value) { limit = value; uniforms.limit.value = value; },
    limit() { return limit; },
    update(camera, state) {
      object.position.copy(camera.position);
      object.updateMatrixWorld();
      uniforms.skyRadius.value = Math.max(1e-6, camera.near + SKY_RADIUS_SHARE * (camera.far - camera.near));
      if (state && Number.isFinite(state.years)) uniforms.years.value = state.years;
      if (state && state.pixelRatio) uniforms.pixelRatio.value = state.pixelRatio;
    },
    dispose() {
      geometry.dispose();
      material.dispose();
    },
  };
}
