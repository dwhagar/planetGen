// tests/js/tileview.test.mjs -- static/tileview.js (MAP.249): a view is known from its tiles alone. The guard here is
// that the pass never reads a tile's stars, so its cost follows the number of tiles (a view holds about 250) and not
// the number of stars (up to 45,000), and that it still tells a changed view from an unchanged one.

import { test } from "node:test";
import assert from "node:assert/strict";

import { importPage } from "./fakedom.mjs";

const V = await importPage("tileview.js");

// A tile whose star arrays count every read, so a pass over the stars shows.
function tile(starCount, reads) {
  const guard = (list) => new Proxy(list, { get(target, name) { reads.count += 1; return target[name]; } });
  return { stars: guard(new Array(starCount).fill(null).map((_, i) => ({ id: i }))), generated: guard([]), points: guard([]) };
}

function view(tileCount, starsPerTile) {
  const reads = { count: 0 };
  const tiles = new Map();
  const keys = [];
  for (let i = 0; i < tileCount; i += 1) {
    const key = `7/${i}/0/0`;
    keys.push(key);
    tiles.set(key, tile(starsPerTile, reads));
  }
  return { keys, tiles, reads, getTile: (key) => tiles.get(key) };
}

test("equal views have equal signatures, in any order, and a replaced or missing tile changes it", () => {
  const v = view(6, 3);
  const first = V.drawnSignature(v.keys, new Set(), v.getTile);
  assert.deepEqual(first.missing, []);
  assert.equal(V.drawnSignature([...v.keys].reverse(), new Set(), v.getTile).signature, first.signature);
  assert.notEqual(V.drawnSignature(v.keys, new Set([v.keys[0]]), v.getTile).signature, first.signature, "detail is part of it");
  const replaced = new Map(v.tiles);
  replaced.set(v.keys[2], { stars: [], generated: [], points: [] });
  assert.notEqual(V.drawnSignature(v.keys, new Set(), (key) => replaced.get(key)).signature, first.signature);
  const gone = V.drawnSignature(v.keys, new Set(), (key) => (key === v.keys[1] ? undefined : v.tiles.get(key)));
  assert.deepEqual(gone.missing, [v.keys[1]]);
  assert.notEqual(gone.signature, first.signature);
});

test("a view of 250 tiles and 45,000 stars is signed without reading a star, in a few milliseconds", () => {
  const v = view(250, 180);
  V.drawnSignature(v.keys, new Set(), v.getTile);   // the revisions are numbered the first time
  const started = performance.now();
  const rounds = 20;
  for (let i = 0; i < rounds; i += 1) V.drawnSignature(v.keys, new Set(), v.getTile);
  const each = (performance.now() - started) / rounds;
  assert.equal(v.reads.count, 0, "the signature read a tile's stars");
  assert.ok(each < 5, `a pass over 250 tiles took ${each.toFixed(2)} ms; the old walk over 45,000 stars took about 160 ms`);
});

// --- The stars carried into tiles that have not arrived (MAP.249) -------------------------------------------------

const ROOT = 524288;

function tileBoxOf(key) {
  const [level, i, j, k] = key.split("/").map(Number);
  const edge = ROOT / 2 ** level;
  const lo = [-ROOT / 2 + i * edge, -ROOT / 2 + j * edge, -ROOT / 2 + k * edge];
  return [lo, lo.map((v) => v + edge)];
}

function randomEntries(count, spread, seed) {
  let state = seed;
  const next = () => { state = (state * 1664525 + 1013904223) % 4294967296; return state / 4294967296; };
  return Array.from({ length: count }, (_, id) => ({
    id: `b${id}@x`, star: { x: (next() - 0.5) * spread, y: (next() - 0.5) * spread, z: (next() - 0.5) * spread },
  }));
}

test("the stars carried into missing tiles are the ones the old walk over every tile found, once each", () => {
  const entries = randomEntries(4000, 6000, 7);
  // Missing tiles of three levels, overlapping, around the middle of the galaxy; a tile holding the middle at each.
  const missing = [];
  for (const level of [6, 8, 10]) {
    const edge = ROOT / 2 ** level;
    const base = Math.floor((0 + ROOT / 2) / edge);
    for (let i = -2; i <= 2; i += 1) for (let j = -1; j <= 1; j += 1) missing.push(`${level}/${base + i}/${base + j}/${base}`);
  }
  const taken = new Set(["b3@x", "b10@x"]);
  const found = V.entriesInTiles(entries, missing, ROOT, (id) => taken.has(id));
  const boxes = missing.map(tileBoxOf);
  const inside = (star, [lo, hi]) => [0, 1, 2].every((a) => [star.x, star.y, star.z][a] >= lo[a] && [star.x, star.y, star.z][a] < hi[a]);
  const wanted = entries.filter((entry) => !taken.has(entry.id) && boxes.some((box) => inside(entry.star, box)));
  assert.ok(wanted.length > 20, "the case has stars to carry");
  assert.deepEqual(found.map((entry) => entry.id), wanted.map((entry) => entry.id));
  assert.deepEqual(V.entriesInTiles(entries, [], ROOT, () => false), []);
});

test("carrying stars into 250 missing tiles from 45,000 drawn stars takes a pass, not a walk per tile", () => {
  const entries = randomEntries(45000, 3000, 11);
  const missing = [];
  const level = 9;
  const base = Math.floor(ROOT / 2 / (ROOT / 2 ** level));
  for (let i = -3; i < 3; i += 1) for (let j = -3; j < 3; j += 1) for (let k = -3; k < 3; k += 1) if (missing.length < 250) missing.push(`${level}/${base + i}/${base + j}/${base + k}`);
  V.entriesInTiles(entries, missing, ROOT, () => false);   // warm the engine
  const started = performance.now();
  const found = V.entriesInTiles(entries, missing, ROOT, () => false);
  const took = performance.now() - started;
  assert.ok(found.length > 1000);
  assert.ok(took < 60, `${took.toFixed(1)} ms; the walk per tile took about 665 ms on this case`);
});
