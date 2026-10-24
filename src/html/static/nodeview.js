// html/static/nodeview.js
//
// MAP.249: what the Galaxy Map decides about a view from its nodes alone. A node is read once and never changes
// (a replaced node is a new object), so whether the drawn stars are the ones a view needs follows from which nodes
// it holds and in which revision, with no walk over the tens of thousands of stars in them. galaxymap3d.js asks
// here before it builds anything: a view whose nodes are the ones already drawn costs a pass over the node keys.

const revisions = new WeakMap();
let lastRevision = 0;

// A number that names this node object; a node that arrives again (refetched, read from storage anew) is another.
export function nodeRevision(node) {
  let revision = revisions.get(node);
  if (revision === undefined) {
    lastRevision += 1;
    revision = lastRevision;
    revisions.set(node, revision);
  }
  return revision;
}

// The nodes of a view: `keys`, `getNode(key)` -> the node or undefined. Returns `signature` (the held nodes by key and
// revision, equal for equal views) and `missing`, the keys no node is held for yet. The order of `keys` does not
// matter to the signature.
export function drawnSignature(keys, getNode) {
  const missing = [];
  const held = [];
  for (let i = 0; i < keys.length; i += 1) {
    const key = keys[i];
    const node = getNode(key);
    if (node === undefined) {
      missing.push(key);
    } else {
      held.push(key + "#" + nodeRevision(node));
    }
  }
  return { signature: held.sort().join(";"), missing };
}
