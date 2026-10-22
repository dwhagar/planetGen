// html/static/tileview.js
//
// MAP.249: what the Galaxy Map decides about a view from its tiles alone. A tile is read once and never changes
// (a replaced tile is a new object), so whether the drawn stars are the ones a view needs follows from which tiles
// it holds and in which revision, with no walk over the tens of thousands of stars in them. galaxymap3d.js asks
// here before it builds anything: a view whose tiles are the ones already drawn costs a pass over the tile keys.

const revisions = new WeakMap();
let lastRevision = 0;

// A number that names this tile object; a tile that arrives again (refetched, read from storage anew) is another.
export function tileRevision(tile) {
  let revision = revisions.get(tile);
  if (revision === undefined) {
    lastRevision += 1;
    revision = lastRevision;
    revisions.set(tile, revision);
  }
  return revision;
}

// The tiles of a view, in order: `keys` (the finest, `detail` ones first), `getTile(key)` -> the tile or undefined.
// Returns `signature` (the held tiles by key, detail mark and revision, equal for equal views) and `missing`, the
// keys no tile is held for yet. The order of `keys` does not matter to the signature.
export function drawnSignature(keys, detail, getTile) {
  const missing = [];
  const held = [];
  for (let i = 0; i < keys.length; i += 1) {
    const key = keys[i];
    const tile = getTile(key);
    if (tile === undefined) {
      missing.push(key);
    } else {
      held.push(key + (detail && detail.has(key) ? "d" : "") + "#" + tileRevision(tile));
    }
  }
  return { signature: held.sort().join(";"), missing };
}

// The entries (`{id, star, ...}`, the stars drawn now) that lie in one of the tiles named by `missing`, to carry over
// until those tiles arrive (MAP.48). One pass over the entries: each is placed in a tile of each level the missing
// tiles come at, by arithmetic, and kept when that tile is missing. Before, every missing tile walked all the
// entries (about 250 tiles by 45,000 stars, 11 million box tests, on each of the several calls a zoom makes).
// `root` is the root tile's edge (the tiles tile the cube centred on the origin); `taken(id)` says an entry is
// already drawn from a present tile. Returns the entries in the order given, each once.
export function entriesInTiles(entries, missing, root, taken) {
  if (!missing.length || !entries.length) {
    return [];
  }
  const parsed = missing.map((key) => key.split("/").map(Number));
  const top = parsed.reduce((most, parts) => Math.max(most, parts[0]), 0);
  const side = Math.pow(2, top);
  const code = (level, i, j, k) => ((level * side + i) * side + j) * side + k;
  const wanted = new Set(parsed.map((parts) => code(parts[0], parts[1], parts[2], parts[3])));
  const levels = Array.from(new Set(parsed.map((parts) => parts[0])));
  const origin = -root / 2;
  const edges = levels.map((level) => root / Math.pow(2, level));
  const carried = [];
  for (let n = 0; n < entries.length; n += 1) {
    const entry = entries[n];
    if (taken(entry.id)) {
      continue;
    }
    const star = entry.star;
    for (let l = 0; l < levels.length; l += 1) {
      const edge = edges[l];
      const i = Math.floor((star.x - origin) / edge);
      const j = Math.floor((star.y - origin) / edge);
      const k = Math.floor((star.z - origin) / edge);
      if (i >= 0 && j >= 0 && k >= 0 && i < side && j < side && k < side && wanted.has(code(levels[l], i, j, k))) {
        carried.push(entry);
        break;
      }
    }
  }
  return carried;
}
