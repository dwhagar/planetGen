// html/static/superbubble.js
//
// MAP.236: how a superbubble (GEN.254) is drawn on the Galaxy Map. A bubble was three faint great circles, so its bound read as two
// crossing ellipses and the inside had no edge. It is now the sphere itself, seen from inside its far wall: a faint, even fill
// (more transparent than the old circles, 0.14) and a crisp rim a couple of pixels wide where the sphere meets the line of sight,
// brighter than the fill, so the true silhouette is the bound. Back faces only, so a camera outside sees one layer of fill and a
// camera inside the bubble (the Sun is in one) sees the shell round it, never two layers.
//
// Pure data and shader text, free of three.js, so a test can hold the numbers; galaxymap3d.js builds the material.

// The look in each theme: the colour (CSS hex), the fill's alpha and the rim's alpha. The rim is always the brighter of the two, and
// the fill is always under the 0.14 the three circles had.
export const DARK_LOOK = { color: "#6fb3d8", fillAlpha: 0.06, rimAlpha: 0.7 };
export const LIGHT_LOOK = { color: "#1f6a99", fillAlpha: 0.07, rimAlpha: 0.75 };

// How wide the rim is, in the shader's own measure (about two pixels on a bubble that fills the map).
export const RIM_PX = 4.0;

export function bubbleLook(light) {
  return light ? LIGHT_LOOK : DARK_LOOK;
}

// The colour as [r, g, b] in 0..1, as written to the screen (the shader writes it as it is).
export function colorComponents(hex) {
  const value = parseInt(hex.replace("#", ""), 16);
  return [((value >> 16) & 255) / 255, ((value >> 8) & 255) / 255, (value & 255) / 255];
}

export const VERTEX_SHADER = `
varying vec3 vNormal;
varying vec3 vView;
void main() {
  vec4 viewPosition = modelViewMatrix * vec4(position, 1.0);
  vNormal = normalize(normalMatrix * normal);
  vView = -viewPosition.xyz;
  gl_Position = projectionMatrix * viewPosition;
}
`;

// `facing` is 1 where the sphere's wall faces the eye and 0 where it meets the line of sight (the silhouette); the rim is the
// pixels within `rimPx` of that edge, measured by how fast `facing` changes across the screen.
export const FRAGMENT_SHADER = `
uniform vec3 color;
uniform float fillAlpha;
uniform float rimAlpha;
uniform float rimPx;
varying vec3 vNormal;
varying vec3 vView;
void main() {
  float facing = abs(dot(normalize(vNormal), normalize(vView)));
  float width = max(fwidth(facing), 1e-4) * rimPx;
  float rim = 1.0 - smoothstep(0.0, width, facing);
  gl_FragColor = vec4(color, mix(fillAlpha, rimAlpha, rim));
}
`;
