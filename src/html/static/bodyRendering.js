// html/static/bodyRendering.js
//
// Shared three.js building blocks for rendering a celestial body as a real
// textured, glowing sphere -- used by both static/systemmap.js (planets/
// moons/stars, one marker at a time via a scissored shared canvas) and
// static/sectormap.js (every star/nebula/asteroid field/black hole/
// neutron star in a sector, all live in the same free-flying 3D scene).
// Factored out here rather than duplicated in both files once sectormap.js
// needed the exact same star-granulation texture and fresnel glow shader
// systemmap.js had already built for its own star spheres.
//
// Every export is a plain function/string constant -- no shared mutable
// state, no THREE.js import of its own (the caller's own `THREE` module
// instance is passed in explicitly) so this stays usable from either
// file's own `import * as THREE from "./vendor/three.module.min.js"`
// without risking two separate THREE instances ever existing in one page.

// A standard fresnel rim-glow: brightest where the surface normal points
// away from the camera (the limb), near-zero head-on -- rendered on a
// slightly larger, back-face-only, additively-blended sphere around the
// body itself, the common cheap "atmosphere/corona" trick (no real
// scattering simulation, just a glow that reads as one). glowPower/
// glowStrength are per-material uniforms (not baked into the shader
// source) so one shared ShaderMaterial can be reused for a subtle planet
// atmosphere rim in one place and a bright star/nebula corona in another,
// just by setting different uniform values -- see makeGlowMaterial below.
export var GLOW_VERTEX_SHADER = [
  "varying vec3 vNormal;",
  "varying vec3 vViewDir;",
  "void main() {",
  "  vNormal = normalize(normalMatrix * normal);",
  "  vec4 viewPosition = modelViewMatrix * vec4(position, 1.0);",
  "  vViewDir = normalize(-viewPosition.xyz);",
  "  gl_Position = projectionMatrix * viewPosition;",
  "}",
].join("\n");

export var GLOW_FRAGMENT_SHADER = [
  "uniform vec3 glowColor;",
  "uniform float glowPower;",
  "uniform float glowStrength;",
  "varying vec3 vNormal;",
  "varying vec3 vViewDir;",
  "void main() {",
  "  float rim = 1.0 - max(dot(vNormal, vViewDir), 0.0);",
  "  float intensity = pow(rim, glowPower) * glowStrength;",
  "  gl_FragColor = vec4(glowColor, intensity);",
  "}",
].join("\n");

// Builds one glow ShaderMaterial -- `THREE` is the caller's own module
// instance (see this file's own header comment for why it's passed in
// rather than imported here). `glowPower` (lower = glow spreads in
// further from the limb, not just a thin rim) and `glowStrength` (overall
// brightness) both default to a subtle "planet atmosphere" reading;
// callers wanting a brighter/wider "this is a light source" corona (a
// star, a neutron star, an accreting black hole) pass their own values.
export function makeGlowMaterial(THREE, colorHex, glowPower, glowStrength) {
  return new THREE.ShaderMaterial({
    uniforms: {
      glowColor: { value: new THREE.Color(colorHex) },
      glowPower: { value: glowPower != null ? glowPower : 2.5 },
      glowStrength: { value: glowStrength != null ? glowStrength : 1.0 },
    },
    vertexShader: GLOW_VERTEX_SHADER,
    fragmentShader: GLOW_FRAGMENT_SHADER,
    side: THREE.BackSide,
    blending: THREE.AdditiveBlending,
    transparent: true,
    depthWrite: false,
  });
}

// A star's photosphere reads as mottled granulation, not a clean flat
// disc -- varies in BOTH UV directions (unlike a gas giant's own height-
// only banding, correct for a horizontally-banded planet but not a star):
// several layered sine "octaves" at different frequencies/phases/axes,
// cross-modulated against each other, give an irregular blotchy pattern
// (bright granulation cells and darker starspot-like patches) rather than
// the regular grid a plain 2D sine product would produce. `THREE` is the
// caller's own module instance, same reasoning as makeGlowMaterial above.
export function makeStarSurfaceTexture(THREE, baseColorHex) {
  var size = 128;
  var canvasEl = document.createElement("canvas");
  canvasEl.width = canvasEl.height = size;
  var ctx = canvasEl.getContext("2d");
  var base = new THREE.Color(baseColorHex);
  var imageData = ctx.createImageData(size, size);
  for (var y = 0; y < size; y++) {
    var v = y / size;
    for (var x = 0; x < size; x++) {
      var u = x / size;
      var n =
        0.5 +
        0.22 * Math.sin(u * Math.PI * 18 + Math.sin(v * 11) * 2.0) +
        0.18 * Math.sin(v * Math.PI * 14 + Math.cos(u * 9) * 2.4) +
        0.12 * Math.sin((u + v) * Math.PI * 23) +
        0.1 * Math.sin((u - v) * Math.PI * 27 + 1.3);
      n = Math.max(0, Math.min(1, n));
      // Kept bright overall (0.72-1.12x) -- a self-luminous surface, so
      // this is texture/granulation, never allowed to read as "shadow".
      var shade = 0.72 + 0.4 * n;
      var idx = (y * size + x) * 4;
      imageData.data[idx] = Math.min(255, Math.round(base.r * 255 * shade));
      imageData.data[idx + 1] = Math.min(255, Math.round(base.g * 255 * shade));
      imageData.data[idx + 2] = Math.min(255, Math.round(base.b * 255 * shade));
      imageData.data[idx + 3] = 255;
    }
  }
  ctx.putImageData(imageData, 0, 0);
  var texture = new THREE.CanvasTexture(canvasEl);
  texture.colorSpace = THREE.SRGBColorSpace;
  return texture;
}
