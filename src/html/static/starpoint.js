// html/static/starpoint.js
//
// MAP.208: how big a star, or a black hole, neutron star or quasar, is drawn
// on the Galaxy Map and in the sector scene it opens into. Far from the
// camera it is a point of light POINT_PX across whose size never changes
// with the zoom, the brightness alone telling one from another; only near
// sector level does it grow, smoothly, into the disc and glow
// static/galaxymap3d.js and static/sectorscene.js give it.
//
// ONE constant sets where that happens: NEAR_SECTOR_PC, the distance from the
// camera at which a point is fully grown. Both scenes use it (the shaders
// below repeat `pointGrow` in GLSL, built from these constants), so the
// hand-off from the Galaxy Map to a sector shows no jump in size (MAP.203
// matched the two once already). The growth runs over the last POINT_RAMP_STEPS
// zoom steps of the map's block ladder (each ZOOM_STEP times the last), in
// the logarithm of the distance:
//
//     grow = 1 - smoothstep(log(d / NEAR_SECTOR_PC) / log(POINT_FAR_PC / NEAR_SECTOR_PC))
//     size = POINT_PX + (fullSize - POINT_PX) * grow
//
// The system map keeps its own scale (its bodies are drawn at their drawn
// size, and it has no background stars of its own).

// The camera distance at which a point is fully grown: a sector is about
// 4 pc across and the camera fits it from about 8 to 12 pc, so every light
// in an opened sector is at its full size with room to zoom out a little.
export const NEAR_SECTOR_PC = 16;
// The block ladder's step (243, 81, 27, 9, 3, 1) and how many steps the growth spans.
export const ZOOM_STEP = 3;
export const POINT_RAMP_STEPS = 2;
// Beyond this distance a point is a point.
export const POINT_FAR_PC = NEAR_SECTOR_PC * Math.pow(ZOOM_STEP, POINT_RAMP_STEPS);
// A point of light's width and the opacity of its core, faint star to bright
// star, when it is far (a share of the star's own brightness grows it back
// to the full look as it nears).
export const POINT_PX = 1.5;
export const POINT_BRIGHT = [0.45, 1];
// Boss (2026-10-23): a star glows at every distance, its brightness and color its own. A far point keeps a halo
// POINT_HALO_PX wider than itself and POINT_HALO_GLOW strong, faint star to bright star (by luminosity share),
// growing into the full halo near a sector; its core is whitened only POINT_WHITEN (a grown core GROWN_WHITEN), so
// the star's blackbody color shows in the point as well as in its halo.
export const POINT_HALO_PX = [3, 12];
export const POINT_HALO_GLOW = [0.3, 0.75];
export const POINT_WHITEN = 0.25;
export const GROWN_WHITEN = 0.6;

function lerpShare(range, share) {
  return range[0] + (range[1] - range[0]) * Math.min(1, Math.max(0, share));
}

// The halo's width beyond the point, in pixels, of a far light of luminosity share `share` (0..1).
export function pointHaloPx(share) {
  return lerpShare(POINT_HALO_PX, share);
}

// The halo's strength of a far light of luminosity share `share` (0..1).
export function pointHaloGlow(share) {
  return lerpShare(POINT_HALO_GLOW, share);
}

function ease01(x) {
  const t = Math.min(1, Math.max(0, x));
  return t * t * (3 - 2 * t);
}

// 0 for a point (at POINT_FAR_PC or beyond), 1 fully grown (at NEAR_SECTOR_PC or nearer).
export function pointGrow(distancePc) {
  if (!(distancePc > 0)) return 1;
  const t = Math.log(distancePc / NEAR_SECTOR_PC) / Math.log(POINT_FAR_PC / NEAR_SECTOR_PC);
  return 1 - ease01(t);
}

// The width, in pixels, of a light whose full width is `fullPx`, seen from
// `distancePc`, `haloPx` (`pointHaloPx`) the far halo around its point. `arrive` (0..1, MAP.155) holds a star
// arriving with the zoom at a point.
export function pointSizePx(fullPx, distancePc, arrive = 1, haloPx = 0) {
  const far = POINT_PX + haloPx;
  return far + Math.max(fullPx - far, 0) * pointGrow(distancePc) * Math.min(1, Math.max(0, arrive));
}

// The core's width the same way, and the core's share of the light's width (the fragment shader's
// `vCore`): 1 for a point with no halo, which is all core.
export function pointCoreShare(fullCorePx, fullPx, distancePc, arrive = 1, haloPx = 0) {
  const grow = pointGrow(distancePc) * Math.min(1, Math.max(0, arrive));
  const core = POINT_PX + Math.max(fullCorePx - POINT_PX, 0) * grow;
  return core / pointSizePx(fullPx, distancePc, arrive, haloPx);
}

// The opacity of a far star's point, from its luminosity share (0..1): the
// only thing that tells one star from another at a distance (apart from color).
export function pointBrightness(luminosityShare) {
  return POINT_BRIGHT[0] + (POINT_BRIGHT[1] - POINT_BRIGHT[0]) * Math.min(1, Math.max(0, luminosityShare));
}

// A phenomenon is drawn as bright as a star of its mass, never brighter
// than the brightest star drawn (`brightestSol`, or none when there is no star):
// the luminosity to draw it at, solar luminosities. The server works out
// `luminositySol` (a tile point's `map_luminosity_sol`, planetgen.physics.stellar_evolution.phenomenon_map_luminosity_sol:
// a star of its mass, a black hole a little dimmer, no mass a one solar mass star).
export function phenomenonLuminositySol(luminositySol, brightestSol) {
  const own = luminositySol > 0 ? luminositySol : 1;
  return brightestSol > 0 ? Math.min(own, brightestSol) : own;
}

// Whether "Charted only" (MAP.111) dims this tile entry: a star outside a charted
// sector, or a scattered phenomenon whose sector is not built yet. A generated star
// and a phenomenon of a built sector are in a filled sector. The server marks a
// scattered point `charted` once its sector is built (`galaxy_scattered_points_in_box`).
export function isUncharted(entry) {
  if (entry.phenomenon) return Boolean(entry.scattered) && !entry.charted;
  return !(entry.generated || entry.charted);
}

// How much of its light an uncharted object keeps under "Charted only" (MAP.252), 0 to 1, by the order objects
// are cut in (MAP.251, docs/design/interface/zoom-star-visibility.md section 10): a star of a satellite galaxy (`tier` 2) is
// dimmed hardest, then a halo star (`tier` 1), then two ramps at the same time, a star from small and dim to large
// and bright and a phenomenon from small to large by mass (`rampShare`, 0 to 1: `luminosityShare` of a star,
// `massRampShare` of a phenomenon), the most massive and the brightest keeping the most. A charted object keeps all.
export const UNCHARTED_KEEP_SATELLITE = 0.03;
export const UNCHARTED_KEEP_HALO = 0.08;
export const UNCHARTED_KEEP_RAMP = [0.12, 0.5];
export function unchartedKeep(entry, rampShare) {
  if (!isUncharted(entry)) return 1;
  if (entry.tier === 2) return UNCHARTED_KEEP_SATELLITE;
  if (entry.tier === 1) return UNCHARTED_KEEP_HALO;
  const share = Math.min(1, Math.max(0, rampShare));
  return UNCHARTED_KEEP_RAMP[0] + (UNCHARTED_KEEP_RAMP[1] - UNCHARTED_KEEP_RAMP[0]) * share;
}

// A phenomenon's place on its ramp: its mass in powers of ten over `MASS_RAMP_DECADES` (a solar mass is the
// smallest, the same range as the server's `scattered_point_size`); a phenomenon with no mass sits in the middle.
export const MASS_RAMP_DECADES = 8;
export function massRampShare(massSolar) {
  if (!(massSolar > 0)) return 0.5;
  return Math.min(1, Math.max(0, Math.log10(Math.max(massSolar, 1)) / MASS_RAMP_DECADES));
}

// The distance rule (MAP.252): the farther an uncharted object is from the camera than its focus, the less it keeps:
// at the focus's distance all of `keep`, from four times that far `UNCHARTED_FAR_KEEP` of it (log in between).
export const UNCHARTED_FAR_KEEP = 0.35;
export const UNCHARTED_FAR_RATIO = 4;
export function unchartedFarShare(distancePc, focusPc) {
  if (!(focusPc > 0) || !(distancePc > focusPc)) return 0;
  return Math.min(1, Math.log(distancePc / focusPc) / Math.log(UNCHARTED_FAR_RATIO));
}

// The same `unchartedFarShare` in GLSL: `float unchartedFar(float distancePc, float focusPc)`, the share of the
// distance rule to apply (0 at the focus, 1 four times as far).
export const UNCHARTED_FAR_GLSL = [
  "float unchartedFar(float distancePc, float focusPc) {",
  "  if (focusPc <= 0.0 || distancePc <= focusPc) return 0.0;",
  "  return min(1.0, log(distancePc / focusPc) / " + Math.log(UNCHARTED_FAR_RATIO).toFixed(6) + ");",
  "}",
].join("\n");

// The same `pointGrow` in GLSL: `float pointGrow(float distancePc)`.
export const STAR_POINT_GLSL = [
  "float pointGrow(float distancePc) {",
  "  float t = clamp(log(max(distancePc, 1e-6) / " + NEAR_SECTOR_PC.toFixed(4) + ") / "
    + Math.log(POINT_FAR_PC / NEAR_SECTOR_PC).toFixed(6) + ", 0.0, 1.0);",
  "  return 1.0 - t * t * (3.0 - 2.0 * t);",
  "}",
].join("\n");
