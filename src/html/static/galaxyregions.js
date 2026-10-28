// html/static/galaxyregions.js
//
// MAP.316: the charted regions of the Galaxy Map, the choosing half (galaxyglowlayer.js draws them, galaxyregionlist.js lists
// them). A region is a 9-block that holds a generated sector (`region_stats`, the API's `charted-regions`); at galactic scale
// it is the thing drawn for known space, a soft sphere at the region's extent (the bounding sphere of its largest
// connected clump of generated sectors, at least the 12 pc framing sphere) coloured by the share of the region charted, with
// its name on the info card. Every region above 0 percent is drawn (the names on the map are regionlabels.js's, MAP.314). As the camera comes in and the region's own stars resolve, the sphere
// fades (REGION_RESOLVE_PX), handing the picture to the stars, as the unresolved glow does (MAP.324). No DOM and no three.js here:
// tests/js/galaxyregions.test.mjs holds it to this description.

// The strongest a region is drawn (alpha) at a share charted of 0 and of 1: a lone sector is a faint cool sphere, a filled
// region a warm one, none outshining a star's own core (the glow layer draws with a maximum blend, so overlaps never add).
export const REGION_ALPHA = [0.2, 0.5];
// A region whose sphere is, on the screen, this many pixels in radius (the first number) begins to fade and is gone at the
// second: the camera is among its stars, which a haze would only hide. CSS pixels.
export const REGION_RESOLVE_PX = [160, 600];
// A region's radius in pixels is its own radius on screen, kept between these: a far region stays a mark that can be seen
// and aimed at, a near one never washes the canvas.
export const REGION_RADIUS_PX = [5, 280];
// The colour stops of a region by share charted (its square root, so a standard fill of about 14 percent is already
// plainly warmer than a lone sector): cool blue for a lone sector, teal, then gold and amber for a well filled one. Linear RGB
// 0 to 1.
const COLOR_STOPS = [
  [0, [0.36, 0.52, 1.0]],
  [0.35, [0.3, 0.85, 0.9]],
  [0.7, [0.96, 0.84, 0.38]],
  [1, [1.0, 0.55, 0.24]],
];

function ease(x) {
  const t = Math.min(1, Math.max(0, x));
  return t * t * (3 - 2 * t);
}

// The colour [r, g, b] (0 to 1) of a region that is `share` (0 to 1) charted.
export function shareColor(share) {
  const t = Math.sqrt(Math.min(1, Math.max(0, share || 0)));
  for (let i = 1; i < COLOR_STOPS.length; i++) {
    const [t1, c1] = COLOR_STOPS[i];
    if (t <= t1) {
      const [t0, c0] = COLOR_STOPS[i - 1];
      const f = (t - t0) / (t1 - t0);
      return c0.map((v, k) => v + (c1[k] - v) * f);
    }
  }
  return COLOR_STOPS[COLOR_STOPS.length - 1][1].slice();
}

// The blobs of the regions the API listed (`regions`: the `charted-regions` payload's list), in their order:
// `{ key, block, x, y, z, radiusPc, share, color, name, parent, sectorCount, capacity }`. A region with no extent is left out.
export function regionBlobs(regions) {
  const found = [];
  (regions || []).forEach((region) => {
    if (!region || !region.centre || !(region.radius > 0)) return;
    found.push({
      key: region.key,
      block: { m: region.m, shell: region.shell, band: region.band, wedge: region.wedge },
      x: region.centre[0], y: region.centre[1], z: region.centre[2], radiusPc: region.radius,
      share: region.share, color: shareColor(region.share), name: region.name || null, parent: region.parent || null,
      sectorCount: region.sector_count, capacity: region.capacity,
    });
  });
  return found;
}

// How strongly a region shows from the camera at `camera` ([x, y, z] parsecs) with focal length `focalPx`, as
// `{ alpha, radiusPx, distance }`. The same shape as the glow's `glowLook`, so the one layer draws both.
export function regionLook(blob, camera, focalPx) {
  const dx = blob.x - camera[0];
  const dy = blob.y - camera[1];
  const dz = blob.z - camera[2];
  const distance = Math.max(Math.sqrt(dx * dx + dy * dy + dz * dz), 1e-6);
  const spherePx = (blob.radiusPc * focalPx) / distance;
  const resolved = ease((spherePx - REGION_RESOLVE_PX[0]) / (REGION_RESOLVE_PX[1] - REGION_RESOLVE_PX[0]));
  const strength = REGION_ALPHA[0] + (REGION_ALPHA[1] - REGION_ALPHA[0]) * Math.sqrt(Math.min(1, Math.max(0, blob.share || 0)));
  return {
    alpha: strength * (1 - resolved),
    radiusPx: Math.min(REGION_RADIUS_PX[1], Math.max(REGION_RADIUS_PX[0], spherePx)),
    distance,
    spherePx,
  };
}

// The regions in the order the list shows them: by share charted (most first) or by name (the unnamed last), the cell
// breaking ties so the order is the same on every machine. `by` is "share" (the default) or "name". A copy.
export function sortRegions(regions, by) {
  const list = (regions || []).slice();
  const byKey = (a, b) => (a.key < b.key ? -1 : a.key > b.key ? 1 : 0);
  if (by === "name") {
    list.sort((a, b) => {
      if (!a.name !== !b.name) return a.name ? -1 : 1;
      const d = a.name && b.name ? a.name.localeCompare(b.name) : 0;
      return d || byKey(a, b);
    });
  } else {
    list.sort((a, b) => (b.share - a.share) || byKey(a, b));
  }
  return list;
}

// The percent text of a share, "14 percent charted" (a share under one percent as "under 1 percent"): the list's and
// the tooltip's.
export function chartedText(share) {
  const percent = 100 * (share || 0);
  if (percent > 0 && percent < 1) return "under 1 percent charted";
  return Math.round(percent) + " percent charted";
}
