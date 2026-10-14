// html/static/tilepack.js
//
// MAP.159: reads the packed Galaxy Map tiles `/galaxy/tiles` sends and the
// tile store keeps. planetgen/web/lib/tilepack.py writes them and describes
// the layout; this file is its reader, so the format has one definition on
// each side and a test (tests/test_tilepack.py with tests/js/tilepack.test.mjs)
// holds the two together.
//
// A response is "PGTB", the header's JSON length (u32), the header (stamp,
// generation, ..., `tiles`: [[key, byte length], ...]) and the tiles' bytes
// in that order. A tile is "PGT1", the JSON length, its JSON (clouds, points,
// incomplete) and two star columns (`stars`, `generated`); a star's id is a
// whole number, a printed object ID or anything else JSON holds, whatever
// the list has (tilepack.py describes the three).

export const LUMINOSITY_LOG_RANGE = [-8, 8];
export const TEMPERATURE_LOG_RANGE = [2, 6];
export const RADIUS_LOG_RANGE = [-3, 4];

const MAGIC = "PGT1";
const RESPONSE_MAGIC = "PGTB";
const U16 = 65535;
const U8 = 255;

function magicAt(bytes, at) {
  return String.fromCharCode(bytes[at], bytes[at + 1], bytes[at + 2], bytes[at + 3]);
}

// The value a log code stands for (null for code 0), the twin of
// tilepack.py's `_decode_log`.
function decodeLog(code, range, top) {
  if (code === 0) return null;
  return Math.pow(10, range[0] + (code - 1) / (top - 1) * (range[1] - range[0]));
}

// One column of `count` values read with `Constructor` (a typed array) from
// `bytes` at `offset`; copies, so the offset need not be aligned.
function column(bytes, offset, count, Constructor) {
  const size = count * Constructor.BYTES_PER_ELEMENT;
  const copy = bytes.slice(offset, offset + size);
  if (copy.length !== size) throw new Error("a packed tile ends early");
  return new Constructor(copy.buffer);
}

function readStars(bytes, view, start) {
  let at = start;
  const count = view.getUint32(at, true);
  at += 4;
  if (!count) return { stars: [], end: at };
  const axes = [];
  for (let i = 0; i < 3; i++) {
    axes.push([view.getFloat64(at, true), view.getFloat64(at + 8, true)]);
    at += 16;
  }
  const positions = [];
  for (let i = 0; i < 3; i++) {
    positions.push(column(bytes, at, count, Uint16Array));
    at += 2 * count;
  }
  const luminosity = column(bytes, at, count, Uint16Array);
  at += 2 * count;
  const temperature = bytes.subarray(at, at + count);
  at += count;
  const radius = bytes.subarray(at, at + count);
  at += count;
  const classes = bytes.subarray(at, at + count);
  at += count;
  const charted = bytes.subarray(at, at + count);
  at += count;
  const kind = bytes[at];
  at += 1;
  let ids;
  if (kind === 0) {
    const low = column(bytes, at, count, Uint32Array);
    at += 4 * count;
    const high = column(bytes, at, count, Uint32Array);
    at += 4 * count;
    ids = (i) => low[i] + high[i] * 4294967296;
  } else if (kind === 1) {
    const width = bytes[at];
    at += 1;
    const text = new TextDecoder("ascii").decode(bytes.subarray(at, at + width * count));
    at += width * count;
    ids = (i) => text.substr(i * width, width);
  } else if (kind === 2) {
    const length = view.getUint32(at, true);
    at += 4;
    const list = JSON.parse(new TextDecoder("utf-8").decode(bytes.subarray(at, at + length)));
    at += length;
    ids = (i) => list[i];
  } else {
    throw new Error("a packed tile has an unknown id kind");
  }
  if (at > bytes.length) throw new Error("a packed tile ends early");

  const stars = new Array(count);
  for (let i = 0; i < count; i++) {
    const star = {
      id: ids(i),
      x: axes[0][0] + positions[0][i] * axes[0][1],
      y: axes[1][0] + positions[1][i] * axes[1][1],
      z: axes[2][0] + positions[2][i] * axes[2][1],
      luminosity_sol: decodeLog(luminosity[i], LUMINOSITY_LOG_RANGE, U16) || 0,
    };
    const t = decodeLog(temperature[i], TEMPERATURE_LOG_RANGE, U8);
    if (t !== null) star.temperature_k = t;
    const r = decodeLog(radius[i], RADIUS_LOG_RANGE, U8);
    if (r !== null) star.radius_sol = r;
    if (classes[i]) star.star_type = String.fromCharCode(classes[i]);
    star.charted = charted[i] === 1;
    stars[i] = star;
  }
  return { stars, end: at };
}

// A packed tile (a Uint8Array) as the tile object the page reads: `stars`
// and `generated` (a star's `charted` is true when its system is known), and
// the JSON sections as sent. Throws on bytes that are not one whole tile.
export function unpackTile(bytes) {
  if (bytes.length < 8 || magicAt(bytes, 0) !== MAGIC) throw new Error("not a packed tile");
  const view = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
  const metaLength = view.getUint32(4, true);
  const metaEnd = 8 + metaLength;
  if (metaEnd > bytes.length) throw new Error("a packed tile ends early");
  const tile = JSON.parse(new TextDecoder("utf-8").decode(bytes.subarray(8, metaEnd)));
  let at = metaEnd;
  for (const section of ["stars", "generated"]) {
    const read = readStars(bytes, view, at);
    tile[section] = read.stars;
    at = read.end;
  }
  if (at !== bytes.length) throw new Error("a packed tile has bytes left over");
  return tile;
}

// A `/galaxy/tiles` body (an ArrayBuffer or Uint8Array) as
// `{ header, packed: Map key -> Uint8Array, tiles: { key: tile } }`: the
// header without its `tiles` list, each tile's own bytes (what the tile
// store keeps) and the tiles read.
export function unpackResponse(body) {
  const bytes = body instanceof Uint8Array ? body : new Uint8Array(body);
  if (bytes.length < 8 || magicAt(bytes, 0) !== RESPONSE_MAGIC) throw new Error("not a packed tile response");
  const view = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
  const headerEnd = 8 + view.getUint32(4, true);
  if (headerEnd > bytes.length) throw new Error("a packed tile response ends early");
  const header = JSON.parse(new TextDecoder("utf-8").decode(bytes.subarray(8, headerEnd)));
  const listed = header.tiles || [];
  delete header.tiles;
  const packed = new Map();
  const tiles = {};
  let at = headerEnd;
  for (const [key, length] of listed) {
    if (at + length > bytes.length) throw new Error("a packed tile response ends early");
    const part = bytes.slice(at, at + length);
    packed.set(key, part);
    tiles[key] = unpackTile(part);
    at += length;
  }
  return { header, packed, tiles };
}
