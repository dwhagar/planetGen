// tests/js/picker.test.mjs -- static/picker.js, the one selection every
// map and page shares (NAV.13): select, step out, step in, step sideways,
// the change event and the trail for a breadcrumb.

import { test } from "node:test";
import assert from "node:assert/strict";

import { STATIC_URL } from "./fakedom.mjs";

const { createSelection } = await import(new URL("picker.js", STATIC_URL).href);

// sector > system > (star, planet) ; planet > (moon a, moon b)
const TREE = {
  galaxy: { parent: null, kids: ["sector"] },
  sector: { parent: "galaxy", kids: ["system"] },
  system: { parent: "sector", kids: ["star", "planet"] },
  star: { parent: "system", kids: [] },
  planet: { parent: "system", kids: ["moon-a", "moon-b"] },
  "moon-a": { parent: "planet", kids: [] },
  "moon-b": { parent: "planet", kids: [] },
};

function selection() {
  return createSelection({
    parentOf: (ref) => TREE[ref.id].parent && { id: TREE[ref.id].parent },
    childrenOf: (ref) => TREE[ref.id].kids.map((id) => ({ id })),
    siblingsOf: (ref) => {
      const parent = TREE[ref.id].parent;
      return parent ? TREE[parent].kids.map((id) => ({ id })) : [ref];
    },
    same: (a, b) => a.id === b.id,
  });
}

test("up walks moon to planet to system to sector to the galaxy, then stops", () => {
  const s = selection();
  s.select({ id: "moon-a" });
  const seen = [];
  for (let n = 0; n < 6; n++) {
    const next = s.up();
    if (!next) break;
    seen.push(next.id);
  }
  assert.deepEqual(seen, ["planet", "system", "sector", "galaxy"]);
  assert.equal(s.current().id, "galaxy");
  assert.equal(s.up(), null);
});

test("into takes only a child of the selection", () => {
  const s = selection();
  s.select({ id: "system" });
  assert.equal(s.into({ id: "moon-a" }), null);
  assert.equal(s.current().id, "system");
  assert.equal(s.into({ id: "planet" }).id, "planet");
  assert.equal(s.into({ id: "moon-b" }).id, "moon-b");
});

test("sideways wraps among the siblings, of any kind", () => {
  const s = selection();
  s.select({ id: "star" });
  assert.equal(s.sideways(1).id, "planet");
  assert.equal(s.sideways(1).id, "star");
  assert.equal(s.sideways(-1).id, "planet");
  s.select({ id: "galaxy" });
  assert.equal(s.sideways(1), null);
});

test("one change event per move, saying how; quiet moves and stops are silent", () => {
  const s = selection();
  const events = [];
  const stop = s.onChange((e) => events.push([e.how, e.previous && e.previous.id, e.ref.id]));
  s.select({ id: "planet" });
  s.up();
  s.into({ id: "planet" });
  s.sideways(1);
  s.select({ id: "planet" }, { quiet: true });
  s.select({ id: "planet" });
  assert.deepEqual(events, [
    ["select", null, "planet"],
    ["up", "planet", "system"],
    ["into", "system", "planet"],
    ["sideways", "planet", "star"],
  ]);
  stop();
  s.up();
  assert.equal(events.length, 4);
});

test("trail is the chain from the top down to the selection", () => {
  const s = selection();
  s.select({ id: "moon-b" });
  assert.deepEqual(s.trail().map((r) => r.id), ["galaxy", "sector", "system", "planet", "moon-b"]);
});
