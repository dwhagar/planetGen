// html/static/galaxyglow.js
//
// MAP.266: the unresolved glow of the Galaxy Map, the choosing half (galaxyglowlayer.js draws it). A node lists only its
// brightest stars; the rest of its light, which the camera cannot tell apart, is one soft blob at the centre of that
// light, drawn for the share the cut gives the node (galaxynodecut.weightCalculator) and so handed over to its
// children's own stars and blobs at the rate the weights trade places: a node and its children always add up to one, so
// nothing pops and nothing doubles (docs/design-drafts single-space-map-zoom sections 2.2 and 2.4).
//
// A node's unlisted light is what its aggregate stores (`aggregate.lights`, the pre-placed stars nothing was built
// from, MAP.259; `aggregate.stats.luminosity_sol`, its generated sectors) less the stars the node itself lists. The
// centre of the first part is exact (the stored flux-weighted centre less the listed stars'); the second has no stored
// centre and sits at the middle of the block's box. A node with no unlisted light has no blob.
//
// The look is a tone curve over the blob's apparent magnitude, judged GLOW_MAG magnitudes deeper than a point of light
// because the eye reads a glow by its surface brightness (design 2.4), ending at GLOW_ALPHA_MAX so the glow never
// outshines a star, and fading as its block grows past a few hundred pixels (the stars are resolved by then, MAP.324). No DOM and no three.js here: tests/js/galaxyglow.test.mjs holds it to this description.

const VERSION_QUERY = new URL(import.meta.url).search;
const { geometry } = await import(`./galaxynodecut.js${VERSION_QUERY}`);
const { absoluteMagnitude, apparentMagnitude } = await import(`./starmagnitude.js${VERSION_QUERY}`);

// How many magnitudes deeper than a point a glow is judged.
export const GLOW_MAG = 5;
// The apparent magnitude (after GLOW_MAG) at which a glow begins to show, and the span over which it comes to full strength.
export const GLOW_LIMIT_MAG = 8;
export const GLOW_RAMP_MAG = 6;
// The strongest a glow is drawn: soft, under a star's own core. Glows are drawn with a maximum blend, not an additive
// one, so overlapping glows never add up past this (MAP.324).
export const GLOW_ALPHA_MAX = 0.4;
// The glow fades as the stars it stands for resolve (MAP.324): a blob whose block is, on the screen, this many pixels in
// radius (the first number) begins to fade and is gone at the second, the point at which the camera is among the block's
// own stars and a haze over them only hides them. Radii in CSS pixels, before GLOW_RADIUS_PX caps what is drawn.
export const GLOW_RESOLVE_PX = [100, 420];
// A blob's radius in pixels is its block's own radius on screen, kept between these (a blob under a pixel is a point
// the stars already are; one over the cap would wash the canvas).
export const GLOW_RADIUS_PX = [3, 220];
// The temperature of light with no stored temperature (a generated sector's), Kelvin.
export const GLOW_DEFAULT_TEMPERATURE_K = 5200;

function listedLight(stars) {
  let light = 0;
  const moment = [0, 0, 0];
  (stars || []).forEach((star) => {
    const lum = star.luminosity_sol > 0 ? star.luminosity_sol : 0;
    light += lum;
    moment[0] += lum * star.x;
    moment[1] += lum * star.y;
    moment[2] += lum * star.z;
  });
  return { light, moment };
}

// A node's blob, `{ key, block, x, y, z, luminosity, temperatureK, radiusPc }`, or `null` for a node with no light left
// over. `block` is the parsed key (`{m, shell, band, wedge}`), `edgePc` the sector edge.
export function blobOf(key, block, node, edgePc) {
  if (!node) return null;
  const aggregate = node.aggregate || {};
  const lights = aggregate.lights || null;
  const listedStars = listedLight(node.stars);
  const listedGenerated = listedLight(node.generated);
  const geo = geometry(block, edgePc);
  const middle = [0, 1, 2].map((a) => (geo.lo[a] + geo.hi[a]) / 2);

  let preLight = 0;
  let preMoment = [0, 0, 0];
  let temperature = null;
  if (lights && lights.star_luminosity_sol > 0 && lights.star_centre) {
    preLight = Math.max(0, lights.star_luminosity_sol - listedStars.light);
    preMoment = [0, 1, 2].map((a) => Math.max(0, lights.star_luminosity_sol * lights.star_centre[a] - listedStars.moment[a]));
    // A centre worked out from sums can drift when the listed stars are most of the light; keep it inside the block.
    if (preLight > 0) {
      preMoment = preMoment.map((value, a) => Math.min(Math.max(value / preLight, geo.lo[a]), geo.hi[a]) * preLight);
    }
    temperature = lights.star_temperature_k || null;
  }
  const stats = aggregate.stats || null;
  const generatedTotal = stats && stats.luminosity_sol > 0 ? stats.luminosity_sol : 0;
  const generatedLight = Math.max(0, generatedTotal - listedGenerated.light);

  const total = preLight + generatedLight;
  if (!(total > 0)) return null;
  const centre = [0, 1, 2].map((a) => (preMoment[a] + generatedLight * middle[a]) / total);
  const radiusPc = Math.sqrt(geo.size * geo.size) / 2;
  return {
    key, block, x: centre[0], y: centre[1], z: centre[2], luminosity: total,
    temperatureK: temperature || GLOW_DEFAULT_TEMPERATURE_K, radiusPc,
  };
}

// The blobs of the held nodes `entries` (`[{ key, block, node }]`), in their order.
export function blobsOf(entries, edgePc) {
  const found = [];
  entries.forEach((entry) => {
    const blob = blobOf(entry.key, entry.block, entry.node, edgePc);
    if (blob) found.push(blob);
  });
  return found;
}

function ease(x) {
  const t = Math.min(1, Math.max(0, x));
  return t * t * (3 - 2 * t);
}

// How strongly a blob shows from the camera at `camera` ([x, y, z] parsecs) with focal length `focalPx`, as
// `{ alpha, radiusPx }` before the node's share (the share multiplies `alpha`). A blob the camera is inside or at is
// at full radius cap.
export function glowLook(blob, camera, focalPx) {
  const dx = blob.x - camera[0];
  const dy = blob.y - camera[1];
  const dz = blob.z - camera[2];
  const distance = Math.max(Math.sqrt(dx * dx + dy * dy + dz * dz), 1e-6);
  const magnitude = apparentMagnitude(absoluteMagnitude(blob.luminosity), distance) - GLOW_MAG;
  const blockPx = (blob.radiusPc * focalPx) / distance;
  const unresolved = 1 - ease((blockPx - GLOW_RESOLVE_PX[0]) / (GLOW_RESOLVE_PX[1] - GLOW_RESOLVE_PX[0]));
  const alpha = GLOW_ALPHA_MAX * ease((GLOW_LIMIT_MAG - magnitude) / GLOW_RAMP_MAG) * unresolved;
  const radiusPx = Math.min(GLOW_RADIUS_PX[1], Math.max(GLOW_RADIUS_PX[0], blockPx));
  return { alpha, radiusPx, distance };
}

// The light of a blob a node shows, for the check that a hand-over keeps the total: its luminosity times its share.
export function sharedLight(blob, share) {
  return blob.luminosity * share;
}
