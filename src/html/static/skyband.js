// html/static/skyband.js
//
// VIEW.34 (S16 of the System Sky plan, reports/system-sky-from-the-system-2026-11-04.md): the glow band of the sky a system sees, the light of
// everything dimmer than the limiting magnitude, drawn as a sphere shader round the camera under the stars (skysphere.js, VIEW.33) and the
// sprites (skysprites.js, VIEW.36). The server's glow image (S9) is sampled by direction; the shader dithers the 8-bit output and adds the grain
// of the stars too faint to resolve, so the band does not look like smooth paint.
//
// The glow image (the contract the server builds to, docs/reference/html-interface.md):
//   * an equirectangular map in the sky's own axes (the packed sky's: galactic, x toward the centre, z toward the north pole), `width` x `height`
//     texels, row 0 at the north pole (latitude +90), column 0 at longitude -180 (the texel centres are half a texel in), longitude = atan2(y, x),
//     latitude = asin(z), the left edge on the anti-centre; 256 x 128 by default;
//   * three floats a texel (red, green, blue), linear light, in units of one mag/arcsec^2 of 22: value = 10^(-0.4 (SB - SB_REF)) for the
//     surface brightness SB (V mag/arcsec^2) of that colour channel's share of the unresolved light. The Milky Way at the Sun runs 21.6 to 23.0
//     (sky-view.md), so the band is values of about 0.4 to 1.5. Zero is no light.
//
// What is drawn, per pixel with the interpolated value v (the texels are bilinear, wrapping in longitude and clamped at the poles):
//     level = PEAK (1 - exp(-GAIN v))                          the soft knee: a faint band shows, a bright core does not blow out
//     grain = 1 + GRAIN (n - 0.5) 2                            n a seeded hash of the cell of the sky under the pixel, cells GRAIN_CELL_DEG wide,
//                                                              equal in area (the longitude cells widen toward the poles' cos)
//     out   = level grain + (h - 0.5) / 255                    h a hash of the screen pixel: the dither, one level either way
//
// Stand-ins (GEN.725 rows by the TODO thread): SB_REF 22 is the middle of the quoted range; PEAK, GAIN, GRAIN and GRAIN_CELL_DEG are chosen by eye
// (tier Unverified, tweak level L3), the grain from needed output (the plan: "a seeded noise in the shader"). `THREE` is passed in, so plain
// node can test the rest.

export const SB_REF = 22;
// Display: the brightest the band is drawn (a share of white), and how fast it gets there.
export const PEAK = 0.4;
export const GAIN = 0.5;
// The grain: a relative variation of the light of +-GRAIN, constant over cells of this size on the sky.
export const GRAIN = 0.2;
export const GRAIN_CELL_DEG = 0.12;
export const DEFAULT_WIDTH = 256;
export const DEFAULT_HEIGHT = 128;
// The sphere sits this far between the camera's near and far planes (the depth is not tested, so it only has to be inside them).
export const BAND_RADIUS_SHARE = 0.5;

// The surface brightness (V mag/arcsec^2) of an image value, and the value of a surface brightness.
export function valueOfSurfaceBrightness(sb) {
  return Math.pow(10, -0.4 * (sb - SB_REF));
}

export function surfaceBrightnessOfValue(value) {
  return value > 0 ? SB_REF - 2.5 * Math.log10(value) : Infinity;
}

// The share of white a value is drawn at.
export function displayLevel(value) {
  return value > 0 ? PEAK * (1 - Math.exp(-GAIN * value)) : 0;
}

// The texture coordinates (0..1, u from longitude -180 eastward, v from the north pole down) of a direction.
export function equirectUV(direction) {
  const length = Math.hypot(direction[0], direction[1], direction[2]) || 1;
  const lon = Math.atan2(direction[1] / length, direction[0] / length);
  const lat = Math.asin(Math.min(1, Math.max(-1, direction[2] / length)));
  return [(lon + Math.PI) / (2 * Math.PI), 0.5 - lat / Math.PI];
}

// The direction (unit vector) at the centre of the texel (column, row) of a map `width` x `height`: the inverse of `equirectUV` there.
export function texelDirection(column, row, width, height) {
  const lon = ((column + 0.5) / width) * 2 * Math.PI - Math.PI;
  const lat = (0.5 - (row + 0.5) / height) * Math.PI;
  return [Math.cos(lat) * Math.cos(lon), Math.cos(lat) * Math.sin(lon), Math.sin(lat)];
}

// The bilinear value of channel `channel` (0..2) of `image` ({ width, height, data: Float32Array of width*height*3 }) towards `direction`: the texture
// unit's answer, wrapping in longitude and clamping at the poles.
export function sample(image, direction, channel) {
  const [u, v] = equirectUV(direction);
  const x = u * image.width - 0.5;
  const y = v * image.height - 0.5;
  const x0 = Math.floor(x);
  const y0 = Math.floor(y);
  const fx = x - x0;
  const fy = y - y0;
  const at = (column, row) => {
    const c = ((column % image.width) + image.width) % image.width;
    const r = Math.min(image.height - 1, Math.max(0, row));
    return image.data[3 * (r * image.width + c) + channel];
  };
  return (at(x0, y0) * (1 - fx) + at(x0 + 1, y0) * fx) * (1 - fy) + (at(x0, y0 + 1) * (1 - fx) + at(x0 + 1, y0 + 1) * fx) * fy;
}

// A glow image of a constant value in each channel, or, with `fill(direction) -> [r, g, b]`, of that function of each texel's direction.
export function makeImage(width, height, fill) {
  const data = new Float32Array(width * height * 3);
  for (let row = 0; row < height; row++) {
    for (let column = 0; column < width; column++) {
      const rgb = typeof fill === "function" ? fill(texelDirection(column, row, width, height)) : fill;
      data.set(rgb, 3 * (row * width + column));
    }
  }
  return { width, height, data };
}

// The image as half floats padded to four channels (what the texture holds), via three.js's `DataUtils`.
export function halfFloatData(THREE, image) {
  const out = new Uint16Array(image.width * image.height * 4);
  const one = THREE.DataUtils.toHalfFloat(1);
  for (let i = 0; i < image.width * image.height; i++) {
    for (let k = 0; k < 3; k++) out[4 * i + k] = THREE.DataUtils.toHalfFloat(image.data[3 * i + k]);
    out[4 * i + 3] = one;
  }
  return out;
}

export const VERTEX_SHADER = [
  "varying vec3 vDirection;",
  "void main() {",
  "  vDirection = position;",
  "  gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);",
  "}",
].join("\n");

export const FRAGMENT_SHADER = [
  "uniform sampler2D glow;",
  "uniform float seed;",
  "uniform float grainStrength;",
  "uniform float grainCells;",
  "varying vec3 vDirection;",
  "const float PI = 3.14159265358979;",
  "const float PEAK = " + PEAK.toFixed(6) + ";",
  "const float GAIN = " + GAIN.toFixed(6) + ";",
  "float hash(vec2 p) {",
  "  p = fract(p * vec2(123.34, 456.21));",
  "  p += dot(p, p + 45.32);",
  "  return fract(p.x * p.y);",
  "}",
  "void main() {",
  "  vec3 d = normalize(vDirection);",
  "  float lon = atan(d.y, d.x);",
  "  float lat = asin(clamp(d.z, -1.0, 1.0));",
  "  vec3 value = texture2D(glow, vec2((lon + PI) / (2.0 * PI), 0.5 - lat / PI)).rgb;",
  "  vec3 level = PEAK * (1.0 - exp(-GAIN * value));",
  // Cells of nearly equal area: rows of equal latitude height, each cut into as many cells in longitude as its circle is long.
  "  float rows = grainCells * 0.5;",
  "  float row = floor((lat + PI * 0.5) / PI * rows);",
  "  float rowLat = (row + 0.5) / rows * PI - PI * 0.5;",
  "  float perRow = max(1.0, floor(grainCells * cos(rowLat)));",
  "  float col = floor((lon + PI) / (2.0 * PI) * perRow);",
  "  float n = hash(vec2(col + seed * 17.0, row + seed * 31.0));",
  "  level *= 1.0 + grainStrength * (n - 0.5) * 2.0;",
  "  level += (hash(gl_FragCoord.xy + seed) - 0.5) / 255.0;",
  "  gl_FragColor = vec4(max(level, 0.0), 1.0);",
  "}",
].join("\n");

// The band of `image` as a mesh. `options`: `seed` (the grain's, default 1), `grain` (default GRAIN, 0 for none). Returns
// `{ object, setImage(image), setSeed(seed), update(camera), dispose() }`: `object` goes into the scene (drawn first, depth not tested); `update` keeps
// it on the camera. Rotate `object.quaternion` to turn the sky's galactic axes into the scene's, as for skysphere.js.
export function createSkyBand(THREE, image, options) {
  const o = options || {};
  const make = (img) => {
    const texture = new THREE.DataTexture(halfFloatData(THREE, img), img.width, img.height, THREE.RGBAFormat, THREE.HalfFloatType);
    texture.wrapS = THREE.RepeatWrapping;
    texture.wrapT = THREE.ClampToEdgeWrapping;
    texture.magFilter = THREE.LinearFilter;
    texture.minFilter = THREE.LinearFilter;
    texture.generateMipmaps = false;
    texture.needsUpdate = true;
    return texture;
  };
  let texture = make(image);
  const cells = Math.round(360 / GRAIN_CELL_DEG);
  const uniforms = {
    glow: { value: texture },
    seed: { value: Number.isFinite(o.seed) ? o.seed : 1 },
    grainStrength: { value: Number.isFinite(o.grain) ? o.grain : GRAIN },
    grainCells: { value: cells },
  };
  const material = new THREE.ShaderMaterial({
    uniforms, vertexShader: VERTEX_SHADER, fragmentShader: FRAGMENT_SHADER,
    side: THREE.BackSide, transparent: false, depthTest: false, depthWrite: false,
  });
  const geometry = new THREE.SphereGeometry(1, 48, 24);
  const object = new THREE.Mesh(geometry, material);
  object.frustumCulled = false;
  object.renderOrder = -12;
  object.name = "sky-band";
  return {
    object,
    setImage(next) {
      const old = texture;
      texture = make(next);
      uniforms.glow.value = texture;
      old.dispose();
    },
    setSeed(seed) { uniforms.seed.value = seed; },
    update(camera) {
      object.position.copy(camera.position);
      const radius = Math.max(1e-6, camera.near + BAND_RADIUS_SHARE * (camera.far - camera.near));
      object.scale.setScalar(radius);
      object.updateMatrixWorld();
    },
    dispose() {
      geometry.dispose();
      material.dispose();
      texture.dispose();
    },
  };
}
