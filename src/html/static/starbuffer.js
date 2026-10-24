// html/static/starbuffer.js
//
// MAP.267: the Galaxy Map's one star buffer. Two pieces of arithmetic, both free of the page so a test can hold them.
//
// 1. Relative-to-eye pairs (design docs/design-drafts single-space-map-zoom section 2.5). A star's position is stored as
//    a high and a low float32 (`splitPair`: the high part is the position rounded to float32, the low part what that
//    rounding left, itself a float32). The shader subtracts the camera's own high and low parts, high from high and low
//    from low, and adds the two: `(high - camHigh) + (low - camLow)`. Each subtraction is exact enough that the star
//    keeps under a quarter pixel down to about 1e-6 pc (0.2 AU) from the camera at any distance from the origin, where
//    one float32 would be a parsec off 260 kpc out. `relativeOffset` is that sum as the GPU does it (every operation
//    rounded to float32, `Math.fround`), for the test that measures the error.
//
// 2. A dense slot store (`SlotStore`) for the node-by-node fill. Each node's stars are one chunk at the end of the
//    buffer; a node that leaves is removed by moving the stars at the tail into the holes (swap-remove), so the buffer
//    is always dense from 0 to `n` and a node landing is a write at the end of it, never a rebuild.

export const SINGLE = Math.fround;

// A value as a pair of float32s, `[high, low]`: high is the value rounded to float32, low the rest (also rounded).
export function splitPair(value) {
  const high = SINGLE(value);
  return [high, SINGLE(value - high)];
}

// The offset from the camera to a star on one axis, as the shader works it out in float32:
// `(high - camHigh) + (low - camLow)`.
export function relativeOffset(high, low, camHigh, camLow) {
  return SINGLE(SINGLE(high - camHigh) + SINGLE(low - camLow));
}

// The same offset when the star is one float32 and the camera is taken from it whole (what a single-float buffer does).
export function singleOffset(value, camera) {
  return SINGLE(SINGLE(value) - SINGLE(camera));
}

// Merges sorted slot numbers into `[start, count]` runs of adjacent slots.
export function runsOf(slots) {
  const runs = [];
  slots.forEach((slot) => {
    const last = runs[runs.length - 1];
    if (last && slot === last[0] + last[1]) last[1] += 1;
    else runs.push([slot, 1]);
  });
  return runs;
}

// A dense run of slots 0..n-1 owned by chunks. `add(key, count)` takes the next `count` slots; `remove(key)` frees a
// chunk's slots by moving the tail into them and returns the moves `[{ from, to }]` to apply to every parallel array,
// in order (each move copies slot `from` over slot `to`), after which the store is dense again at the new `n`.
export class SlotStore {
  constructor() {
    this.n = 0;
    this.owners = [];
    this.chunks = new Map();   // key -> { slots: Set<number>, tag }
  }

  has(key) {
    return this.chunks.has(key);
  }

  tagOf(key) {
    const chunk = this.chunks.get(key);
    return chunk ? chunk.tag : undefined;
  }

  keys() {
    return Array.from(this.chunks.keys());
  }

  // Takes `count` slots at the end for `key` (not already held), tagged `tag` (a revision to compare later); returns
  // the first slot.
  add(key, count, tag) {
    if (this.chunks.has(key)) throw new Error(`chunk ${key} is already held`);
    const start = this.n;
    const slots = new Set();
    for (let i = 0; i < count; i += 1) {
      slots.add(start + i);
      this.owners.push(key);
    }
    this.chunks.set(key, { slots, tag });
    this.n += count;
    return start;
  }

  remove(key) {
    const chunk = this.chunks.get(key);
    if (!chunk) return [];
    const moves = [];
    const descending = Array.from(chunk.slots).sort((a, b) => b - a);
    descending.forEach((slot) => {
      const last = this.n - 1;
      if (slot !== last) {
        const owner = this.owners[last];
        const held = this.chunks.get(owner).slots;
        held.delete(last);
        held.add(slot);
        this.owners[slot] = owner;
        moves.push({ from: last, to: slot });
      }
      this.owners.pop();
      this.n -= 1;
    });
    this.chunks.delete(key);
    return moves;
  }

  clear() {
    this.n = 0;
    this.owners = [];
    this.chunks = new Map();
  }

  // The slots of a chunk, ascending.
  slotsOf(key) {
    const chunk = this.chunks.get(key);
    return chunk ? Array.from(chunk.slots).sort((a, b) => a - b) : [];
  }
}
