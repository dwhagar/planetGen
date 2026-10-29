// html/static/galaxynebulae.js
//
// MAP.176: the nebula regions of the Galaxy Map, the choosing half (galaxyglowlayer.js draws them). A node of 9 to 729 sectors lists
// the groups of the nebulae inside it, its child blocks that hold nebulae (`nebula_regions`, MAP.173: count, volume, centre, reach,
// family shares, dominant class), the largest first. A region is a soft sprite at the group's centre, as wide as the group, in the
// colour of the family with the most volume and with a thin ring to tell it from a single nebula. The groups are a fixed grid (the
// block ladder), so a parent's groups are the union of its children's and the picture is the same for everyone at the same zoom:
// a node draws its groups for the share the cut gives it (galaxynodecut.weightCalculator), a block hands over to its children, and
// each group is gone into its own groups as the camera comes in. A region shows from 2 px across, fades out when it fills the
// screen (the node has long since split it), and at most REGION_BUDGET are drawn, the strongest on the screen. No DOM and no
// three.js here: tests/js/galaxynebulae.test.mjs holds it to this description.

// The most regions drawn at once; the weakest of more are not drawn, and the last tenth of the budget fades out so none pops.
export const REGION_BUDGET = 600;
export const BUDGET_FADE_SHARE = 0.1;
// The families a region is made of, in the order the server gives their shares (planetgen.galaxy.nebula_blocks.FAMILIES).
export const FAMILIES = ["diffuse", "emission", "reflection", "planetary", "dark", "supernova-remnant"];
export const FAMILY_LABELS = {
  diffuse: "Diffuse nebulae", emission: "Emission nebulae", reflection: "Reflection nebulae", planetary: "Planetary nebulae",
  dark: "Dark nebulae", "supernova-remnant": "Supernova remnants",
};
// The colour of a supernova remnant (nebulalook.js has no look for one: a remnant is drawn from its own texture).
const REMNANT_COLOUR = { dark: "#f0a050", light: "#c2570c" };
// A region is drawn from this radius in pixels (a nebula is, at 2 px across), its strength follows its count, and it fades
// out between the two radii of RESOLVE_PX: the camera is among its nebulae, which a haze would only hide.
export const MIN_RADIUS_PX = 1;
export const RESOLVE_PX = [220, 520];
export const RADIUS_PX = [2, 320];
export const ALPHA = [0.16, 0.42];
// The count at which a region is as strong as it gets.
export const STRONG_COUNT = 100;

function ease(x) {
  const t = Math.min(1, Math.max(0, x));
  return t * t * (3 - 2 * t);
}

// "#rrggbb" as [r, g, b] in 0..1.
export function colourOf(hex) {
  return [1, 3, 5].map((i) => parseInt(hex.slice(i, i + 2), 16) / 255);
}

// The colour [r, g, b] (0 to 1) of a family in `looks` (nebulalook.js `nebulaLooks`; `light` whether the page is light).
export function familyColour(family, looks, light) {
  if (family === "supernova-remnant") return colourOf(light ? REMNANT_COLOUR.light : REMNANT_COLOUR.dark);
  const look = looks[family] || looks.default;
  return colourOf(look[0]);
}

// The family with the most volume in a region's `shares` ({family: share}), ties by the order of FAMILIES; null with none.
export function dominantFamily(shares) {
  let best = null;
  let bestShare = 0;
  FAMILIES.forEach((family) => {
    const share = shares && shares[family];
    if (share > bestShare) {
      best = family;
      bestShare = share;
    }
  });
  return best;
}

// The blobs of the regions of the nodes in `entries` (`[{key, block, node}]`: the nodes being drawn): one for each group of
// each node's `nebula_regions` whose family is shown (`isShown(family) -> boolean`). A blob is
// `{ key, block, x, y, z, radiusPc, count, volume, family, nebulaClass, color, ring }`: `block` the node it was listed by (the share
// the cut gives that node is the region's), `ring` true for a group of more than one nebula. In the order of the nodes' lists.
export function regionBlobs(entries, looks, light, isShown) {
  const found = [];
  const seen = new Set();
  (entries || []).forEach((entry) => {
    ((entry.node && entry.node.nebula_regions) || []).forEach((region) => {
      if (!region || !region.centre || !(region.radius_pc > 0) || !(region.count > 0) || seen.has(region.key)) return;
      const family = dominantFamily(region.shares);
      if (!family || (isShown && !isShown(family))) return;
      seen.add(region.key);
      found.push({
        key: region.key, block: entry.block, x: region.centre[0], y: region.centre[1], z: region.centre[2],
        radiusPc: region.radius_pc, count: region.count, volume: region.volume_pc3, family, nebulaClass: region.class || null,
        color: familyColour(family, looks, light), ring: region.count > 1,
      });
    });
  });
  return found;
}

// How strongly a region shows from the camera at `camera` ([x, y, z] parsecs) with focal length `focalPx`, as
// `{ alpha, radiusPx, distance, spherePx }`; the same shape as the glow's and the charted regions', so the one layer draws all three.
export function regionLook(blob, camera, focalPx) {
  const dx = blob.x - camera[0];
  const dy = blob.y - camera[1];
  const dz = blob.z - camera[2];
  const distance = Math.max(Math.sqrt(dx * dx + dy * dy + dz * dz), 1e-6);
  const spherePx = (blob.radiusPc * focalPx) / distance;
  const seen = ease((spherePx - MIN_RADIUS_PX / 2) / MIN_RADIUS_PX);
  const resolved = ease((spherePx - RESOLVE_PX[0]) / (RESOLVE_PX[1] - RESOLVE_PX[0]));
  const strong = Math.min(1, Math.log10(Math.max(1, blob.count)) / Math.log10(STRONG_COUNT));
  return {
    alpha: (ALPHA[0] + (ALPHA[1] - ALPHA[0]) * strong) * seen * (1 - resolved),
    radiusPx: Math.min(RADIUS_PX[1], Math.max(RADIUS_PX[0], spherePx)),
    distance,
    spherePx,
  };
}

// Holds `alphas` (an array of strengths, changed in place and returned) to `budget` drawn: the weakest of more are zeroed, the
// last BUDGET_FADE_SHARE of the budget scaled down towards the cut so a region at the edge of the budget does not pop.
// Ties go to the earlier index, so the same input always draws the same regions.
export function cutToBudget(alphas, budget) {
  const limit = budget === undefined ? REGION_BUDGET : budget;
  let positive = 0;
  for (let i = 0; i < alphas.length; i++) if (alphas[i] > 0) positive += 1;
  if (positive <= limit) return alphas;
  const order = [];
  for (let i = 0; i < alphas.length; i++) if (alphas[i] > 0) order.push(i);
  order.sort((p, q) => (alphas[q] - alphas[p]) || (p - q));
  const fadeFrom = Math.floor(limit * (1 - BUDGET_FADE_SHARE));
  for (let rank = 0; rank < order.length; rank++) {
    const i = order[rank];
    if (rank >= limit) alphas[i] = 0;
    else if (rank >= fadeFrom) alphas[i] *= (limit - rank) / (limit - fadeFrom + 1);
  }
  return alphas;
}

// The words of a region for its hover and card: "12 emission nebulae, mostly Emission nebulae" style text.
export function regionText(blob) {
  const noun = blob.count === 1 ? "nebula" : "nebulae";
  return blob.count + " " + noun + ", mostly " + FAMILY_LABELS[blob.family].toLowerCase();
}
