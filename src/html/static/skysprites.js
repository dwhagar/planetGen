// html/static/skysprites.js
//
// VIEW.36 (S18 of the System Sky plan, reports/system-sky-from-the-system-2026-11-04.md): the nebulae and supernova remnants of the sky a system
// sees, drawn as sprites from VIEW.28's data (planetgen/physics/nebula_sprites.py), next to the stars of skysphere.js. A luminous black hole, neutron star
// or quasar is a point in the packed list itself (kind "phenomenon") and is drawn by skysphere.js like any star at its luminosity.
//
// A sprite is a record of the sky route's JSON (`nebulae` of the sky response): `{ id, family, direction: [x, y, z], angularRadiusRad, inside,
// extinctionMag, surfaceBrightness, visible, seen }`: the unit direction of its centre where it was when its light left it, the angular radius
// (`asin(R / d)`), the dimming through its middle, its mean V surface brightness (null for none) and the two verdicts of `sprite_data`.
//
// What is drawn, by what makes the light (nebula_sprites.py):
//   * a nebula with a surface brightness under the naked-eye limit (the record's `visible`): a soft disc in its family's colour. The ball is a uniform
//     one, so what a line of sight crosses is a chord, sqrt(1 - r^2) at a share r of the radius: the disc follows it, scaled so its mean over the
//     disc is the surface brightness. Its strength is SB_RAMP mag/arcsec^2 of smoothstep to GLOW_MAX under the limit (the sky does not blind);
//   * a supernova remnant that has been seen: the same law through a shell SHELL_THICKNESS of its radius thick, so a ring;
//   * a cloud that gives no light of its own but dims what lies behind it (a dark cloud, a diffuse one with extinction): black, with the opacity
//     1 - 10^(-0.4 A sqrt(1 - r^2)) of its chord, drawn after the stars so it dims them. It dims every star in its direction, near or far, which
//     is right for a star behind it and an approximation for one in front (the per-star extinction is VIEW.8's dust);
//   * a remnant not yet seen, one that is too faint, one that holds the observer (`inside`: it covers the sky, the haze is the dust of VIEW.8) and one
//     wider than MAX_ANGLE (a quad on a tangent plane cannot cover it) are not drawn.
//
// Stand-ins: SB_RAMP and GLOW_MAX are chosen by eye (tier Unverified, L3), like skylight.js's; the limit (23) and the shell thickness (1/12) are
// nebula_sprites.py's registered values. `THREE` is passed in, so plain node can test the rest.

const VERSION_QUERY = new URL(import.meta.url).search;
const { NEBULA_LOOKS } = await import(`./nebulalook.js${VERSION_QUERY}`);
const { colourOf, familyColour } = await import(`./galaxynebulae.js${VERSION_QUERY}`);

// The mean surface brightness, V mag/arcsec^2, under which a sprite is visible (nebula_sprites.NAKED_EYE_LIMIT_MAG_ARCSEC2).
export const SB_LIMIT = 23;
// A sprite goes from nothing at the limit to GLOW_MAX over this many magnitudes brighter.
export const SB_RAMP = 5;
// The brightest a sprite's mean is drawn (a share of its colour).
export const GLOW_MAX = 0.55;
// A shell's thickness over its radius (nebula_sprites.SHELL_THICKNESS_OVER_RADIUS).
export const SHELL_THICKNESS = 1 / 12;
// A sprite wider than this (radians) is not drawn.
export const MAX_ANGLE = (80 * Math.PI) / 180;
// Under this a cloud's dimming is not worth a sprite (magnitudes).
export const MIN_EXTINCTION = 0.02;
export const REMNANT_COLOUR = "#f0a050";
export const REMNANT_FAMILIES = ["supernova-remnant", "plerion"];

function ease(x) {
  const t = Math.min(1, Math.max(0, x));
  return t * t * (3 - 2 * t);
}

// The strength of a sprite of mean surface brightness `sb` (mag/arcsec^2): GLOW_MAX at SB_RAMP under the limit and brighter, nothing at the limit.
export function glowStrength(sb) {
  return Number.isFinite(sb) ? GLOW_MAX * ease((SB_LIMIT - sb) / SB_RAMP) : 0;
}

// The mean over the unit disc of the chord through a shell from `inner` to 1 (of the ball when inner is 0): (2/3)(1 - inner^3).
export function chordMean(inner) {
  return (2 / 3) * (1 - inner * inner * inner);
}

// What a sprite looks like at a share `r` (0 centre, 1 edge) of its radius, as a share of its own mean: the chord through the shell, over its mean.
export function profile(r, inner) {
  if (r >= 1) return 0;
  const through = Math.sqrt(1 - r * r) - (r < inner ? Math.sqrt(inner * inner - r * r) : 0);
  return through / chordMean(inner);
}

// The opacity of a cloud of dimming `extinctionMag` (through its middle) at a share `r` of its radius: 1 - 10^(-0.4 A sqrt(1 - r^2)).
export function darkOpacity(extinctionMag, r) {
  if (r >= 1) return 0;
  return 1 - Math.pow(10, -0.4 * extinctionMag * Math.sqrt(1 - r * r));
}

// How a record is drawn: `{ kind: "glow" | "dark" | "none", strength, inner, colour: [r, g, b], extinctionMag }`.
export function spriteLook(record) {
  const none = { kind: "none", strength: 0, inner: 0, colour: [0, 0, 0], extinctionMag: 0 };
  if (!record || record.inside || record.seen === false || !(record.angularRadiusRad > 0) || record.angularRadiusRad > MAX_ANGLE) return none;
  const remnant = REMNANT_FAMILIES.includes(record.family);
  if (Number.isFinite(record.surfaceBrightness)) {
    const strength = record.visible ? glowStrength(record.surfaceBrightness) : 0;
    if (!(strength > 0)) return none;
    const colour = remnant ? colourOf(REMNANT_COLOUR) : familyColour(record.family, NEBULA_LOOKS);
    return { kind: "glow", strength, inner: remnant ? 1 - SHELL_THICKNESS : 0, colour, extinctionMag: 0 };
  }
  if (!(record.extinctionMag >= MIN_EXTINCTION)) return none;
  return { kind: "dark", strength: 1, inner: 0, colour: [0, 0, 0], extinctionMag: record.extinctionMag };
}

// The two unit vectors across a sprite's tangent plane at `direction`.
export function tangentBasis(direction) {
  const d = direction;
  const length = Math.hypot(d[0], d[1], d[2]) || 1;
  const n = [d[0] / length, d[1] / length, d[2] / length];
  const helper = Math.abs(n[2]) < 0.9 ? [0, 0, 1] : [1, 0, 0];
  const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
  const unit = (v) => { const l = Math.hypot(v[0], v[1], v[2]); return [v[0] / l, v[1] / l, v[2] / l]; };
  const e1 = unit(cross(n, helper));
  return { n, e1, e2: cross(n, e1) };
}

const CORNERS = [[-1, -1], [1, -1], [1, 1], [-1, 1]];

// The geometry arrays of `looks` ([{ record, look }]) as quads: `position` (3 a vertex, a point of the tangent plane at unit distance), `offset` (the
// vertex's offset in that plane, 2), `theta` (1), `colour` (3), `shape` (strength, inner, extinction: 3), `index`.
export function quads(looks) {
  const n = looks.length;
  const position = new Float32Array(12 * n);
  const offset = new Float32Array(8 * n);
  const theta = new Float32Array(4 * n);
  const colour = new Float32Array(12 * n);
  const shape = new Float32Array(12 * n);
  const index = new Uint32Array(6 * n);
  looks.forEach(({ record, look }, i) => {
    const t = record.angularRadiusRad;
    const k = Math.tan(t);
    const { n: dir, e1, e2 } = tangentBasis(record.direction);
    CORNERS.forEach(([cx, cy], c) => {
      const v = 4 * i + c;
      for (let a = 0; a < 3; a++) position[3 * v + a] = dir[a] + k * (cx * e1[a] + cy * e2[a]);
      offset[2 * v] = k * cx;
      offset[2 * v + 1] = k * cy;
      theta[v] = t;
      for (let a = 0; a < 3; a++) colour[3 * v + a] = look.colour[a];
      shape[3 * v] = look.strength;
      shape[3 * v + 1] = look.inner;
      shape[3 * v + 2] = look.extinctionMag;
    });
    index.set([4 * i, 4 * i + 1, 4 * i + 2, 4 * i, 4 * i + 2, 4 * i + 3], 6 * i);
  });
  return { position, offset, theta, colour, shape, index };
}

const VERTEX = [
  "attribute vec2 offset;",
  "attribute float theta;",
  "attribute vec3 colour;",
  "attribute vec3 shape;",
  "uniform float spriteRadius;",
  "varying vec2 vOffset;",
  "varying float vTheta;",
  "varying vec3 vColour;",
  "varying vec3 vShape;",
  "void main() {",
  "  vOffset = offset;",
  "  vTheta = theta;",
  "  vColour = colour;",
  "  vShape = shape;",
  "  gl_Position = projectionMatrix * modelViewMatrix * vec4(position * spriteRadius, 1.0);",
  "}",
].join("\n");

// The share of the radius at the pixel: the angle from the sprite's centre (the atan of the offset in the tangent plane at unit distance) over its radius.
const SHARE = [
  "varying vec2 vOffset;",
  "varying float vTheta;",
  "varying vec3 vColour;",
  "varying vec3 vShape;",
  "float share() { return atan(length(vOffset)) / vTheta; }",
].join("\n");

export const GLOW_FRAGMENT = [
  SHARE,
  "void main() {",
  "  float r = share();",
  "  if (r >= 1.0) discard;",
  "  float inner = vShape.y;",
  "  float through = sqrt(1.0 - r * r) - (r < inner ? sqrt(inner * inner - r * r) : 0.0);",
  "  float mean = (2.0 / 3.0) * (1.0 - inner * inner * inner);",
  "  float edge = 1.0 - smoothstep(0.96, 1.0, r);",
  "  float strength = min(1.0, vShape.x * through / mean) * edge;",
  "  gl_FragColor = vec4(vColour * strength, 1.0);",
  "}",
].join("\n");

export const DARK_FRAGMENT = [
  SHARE,
  "void main() {",
  "  float r = share();",
  "  if (r >= 1.0) discard;",
  "  float alpha = 1.0 - pow(10.0, -0.4 * vShape.z * sqrt(1.0 - r * r));",
  "  gl_FragColor = vec4(0.0, 0.0, 0.0, alpha * (1.0 - smoothstep(0.96, 1.0, r)));",
  "}",
].join("\n");

// The sprites of `records` as meshes: `{ objects, glow, dark, update(camera), dispose() }`. `objects` (the glow mesh, drawn before the stars, and the dark
// one, after them) go into the scene; `update` follows the camera, as skysphere.js's does. Rotate each object's quaternion for the scene's axes.
export function createSkySprites(THREE, records, options) {
  const looks = (records || []).map((record) => ({ record, look: spriteLook(record) })).filter((entry) => entry.look.kind !== "none");
  const uniforms = { spriteRadius: { value: 1 } };
  const make = (kind, fragmentShader, blending, renderOrder) => {
    const own = looks.filter((entry) => entry.look.kind === kind);
    if (!own.length) return null;
    const data = quads(own);
    const geometry = new THREE.BufferGeometry();
    geometry.setAttribute("position", new THREE.BufferAttribute(data.position, 3));
    geometry.setAttribute("offset", new THREE.BufferAttribute(data.offset, 2));
    geometry.setAttribute("theta", new THREE.BufferAttribute(data.theta, 1));
    geometry.setAttribute("colour", new THREE.BufferAttribute(data.colour, 3));
    geometry.setAttribute("shape", new THREE.BufferAttribute(data.shape, 3));
    geometry.setIndex(new THREE.BufferAttribute(data.index, 1));
    const material = new THREE.ShaderMaterial({
      uniforms, vertexShader: VERTEX, fragmentShader, transparent: true, depthTest: false, depthWrite: false, blending, side: THREE.DoubleSide,
    });
    const mesh = new THREE.Mesh(geometry, material);
    mesh.frustumCulled = false;
    mesh.renderOrder = renderOrder;
    mesh.userData.count = own.length;
    return mesh;
  };
  const glow = make("glow", GLOW_FRAGMENT, THREE.AdditiveBlending, -11);
  const dark = make("dark", DARK_FRAGMENT, THREE.NormalBlending, -9);
  glow && (glow.name = "sky-nebulae");
  dark && (dark.name = "sky-dust");
  const objects = [glow, dark].filter(Boolean);
  return {
    objects, glow, dark, count: { glow: glow ? glow.userData.count : 0, dark: dark ? dark.userData.count : 0 },
    // The sprites sit on a sphere nearer than the stars', so the corners of the widest quad stay inside the far plane.
    update(camera) {
      const radius = Math.max(camera.near * 2, camera.far / 12);
      uniforms.spriteRadius.value = radius;
      objects.forEach((object) => { object.position.copy(camera.position); object.updateMatrixWorld(); });
    },
    dispose() {
      objects.forEach((object) => { object.geometry.dispose(); object.material.dispose(); });
    },
  };
}
