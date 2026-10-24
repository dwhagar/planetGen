// MAP.153, MAP.152 and MAP.261: stars fade in and out as the camera nears or leaves the node that lists them.
//
// The map draws the sector-block nodes of its cut (static/galaxynodecut.js: big blocks far, small ones near) and hands
// over from a node to its children over a band of projected size, so nothing pops. A node is shown for what its
// ancestors have left of it (`weight`, galaxynodecut.js), and across a hand-over the weights of a node and its children
// add up to 1: a star listed in both stays whole, one listed in one comes in or goes out smoothly.
//
// A star's `fade` is its node, `{ block: {m, shell, band, wedge} }`, or `null` for a star in no node (a landmark),
// which is whole.

export function smoothstep(x) {
  const t = Math.min(1, Math.max(0, x));
  return t * t * (3 - 2 * t);
}

// A star in no node (a landmark) is whole.
export const POINT_FADE = null;

// A star's opacity from the node that lists it (`fade`, see above): `weigh(block)` is the cut's weight calculator
// (galaxynodecut.weightCalculator) for the camera.
export function starShare(fade, weigh) {
  return fade ? weigh(fade.block) : 1;
}

// MAP.155: a cloud sprite fades in between CLOUD_RAMP_PX[0] and
// CLOUD_RAMP_PX[1] pixels of radius instead of switching on at 2 px.
export const CLOUD_RAMP_PX = [1, 4];

export function cloudSizeShare(radiusPx) {
  return smoothstep((radiusPx - CLOUD_RAMP_PX[0]) / (CLOUD_RAMP_PX[1] - CLOUD_RAMP_PX[0]));
}

// The sum of `starShare` over `fades` for a weight calculator.
export function opacitySum(fades, weigh) {
  let total = 0;
  for (const fade of fades) {
    total += starShare(fade, weigh);
  }
  return total;
}
