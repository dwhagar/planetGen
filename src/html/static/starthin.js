// html/static/starthin.js
//
// MAP.366: the Galaxy Map's tiles keep the same number of stars whatever their size on screen, so where the cut leaves a big
// tile next to a small one at the same distance, the small one draws the same stars in a patch of the screen a tenth of
// the size and the picture shows square chunks of thin and thick. The server lists the brightest NODE_LIST_CAP stars of
// a tile's box (bright and generated lists each), so a tile that reaches the cap is cut by the tile's size and not by the
// view. This is the view's side of it: a tile that reached the cap shows its brightest stars for the share of its list that
// keeps the density on screen near that of the bigger tiles in the view (the top REFERENCE_QUANTILE of the capped tiles
// by size on screen), never below DENSITY_FLOOR of the list, fading the dimmer stars out by their rank in the tile
// instead of dropping them, so the picture changes smoothly as the camera moves. A tile under the cap lists every star it has and
// is left as it is, and so is a phenomenon (a black hole, a pulsar, a quasar).
//
// The constants are interface tuning, not physics (docs/design/interface/zoom-star-visibility.md, MAP.366).

// The most bright stars, and the most generated stars, a block node lists (planetgen.db.nodes.NODE_STARS and
// NODE_GENERATED_STARS; a test holds them equal). A sector node (m = 1) lists everything it has.
export const NODE_LIST_CAP = 96;
// The share of its list a tile that reached the cap keeps at the least, however small it is on screen.
export const DENSITY_FLOOR = 0.3;
// The tile size on screen (pixels, the longest edge of the block) that sets the density: this share of the capped tiles
// are smaller, and a tile at or over it keeps its whole list.
export const REFERENCE_QUANTILE = 0.6;
// The reference never falls under this many pixels (a view of tiny tiles only), and a tile counts as no bigger than PX_CAP
// (one over the screen shows no more of its list for being bigger).
export const REFERENCE_MIN_PX = 150;
export const PX_CAP = 1200;
// The fewest stars of a capped tile shown whatever its size, and the width of the fade from kept to gone, in stars.
export const KEEP_MIN = 4;
export const SOFT_SHARE = 0.3;

function ease01(x) {
  const t = Math.min(1, Math.max(0, x));
  return t * t * (3 - 2 * t);
}

// Ranks the stars of a node by luminosity (0 the brightest): sets `nodeRank` on each star that is not a phenomenon and returns
// `{ listed, capped }`, the number ranked and whether the node reached the cap of either of its lists (`generated` marks
// the star of a generated system). A phenomenon is never thinned and gets no rank.
export function rankStars(stars, block) {
  const ranked = stars.filter((star) => !star.phenomenon);
  ranked.sort((p, q) => (q.luminosity_sol || 0) - (p.luminosity_sol || 0));
  let bright = 0;
  let generated = 0;
  ranked.forEach((star, rank) => {
    star.nodeRank = rank;
    if (star.generated) generated += 1;
    else bright += 1;
  });
  const capped = block.m > 1 && (bright >= NODE_LIST_CAP || generated >= NODE_LIST_CAP);
  return { listed: ranked.length, capped };
}

// The size on screen the density is held to: the REFERENCE_QUANTILE of the capped tiles' sizes (`sizes`, pixels), at
// least REFERENCE_MIN_PX; `null` with no capped tile (nothing is thinned).
export function referencePx(sizes) {
  const sorted = sizes.map((px) => Math.min(PX_CAP, px)).sort((p, q) => p - q);
  if (!sorted.length) return null;
  const at = sorted[Math.min(sorted.length - 1, Math.floor(REFERENCE_QUANTILE * sorted.length))];
  return Math.max(REFERENCE_MIN_PX, at);
}

// The share of its list a capped tile of `px` pixels keeps against the reference: the square of the ratio of sizes (the same
// density on screen), between DENSITY_FLOOR and 1.
export function keptShare(px, reference) {
  if (reference === null) return 1;
  const ratio = Math.min(PX_CAP, px) / reference;
  return Math.min(1, Math.max(DENSITY_FLOOR, ratio * ratio));
}

// The opacity (0 to 1) of the star of `rank` in a tile that lists `listed` stars and keeps `kept` of them: whole up to the
// kept count, then fading out over SOFT_SHARE of it.
export function rankShare(rank, listed, kept) {
  if (kept >= 1) return 1;
  const keep = Math.max(KEEP_MIN, listed * kept);
  const soft = Math.max(2, SOFT_SHARE * keep);
  return ease01((keep - rank) / soft);
}
