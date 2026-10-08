// html/static/picker.js
//
// The one selection every map and page shares (NAV.13): the current
// selection is a reference to an object (a stage, a sector, a system, a
// star, a planet, a moon, ...), and moving it is one of four calls:
//   select(ref)   make `ref` the selection;
//   up()          step out to its parent (moon to planet to system to
//                 sector to the galaxy);
//   into(ref)     step in to one of its children;
//   sideways(dir) move to the next (1) or previous (-1) sibling, wrapping.
// Keys, clicks and the breadcrumb all go through these, and every change
// fires one event to the listeners added with `onChange`, so every display
// follows the same selection.
//
// What the objects ARE, and how they relate, is the resolver's business
// (the page or map that creates the selection):
//   parentOf(ref)    the parent reference, or null at the top;
//   childrenOf(ref)  the children (optional; `into` then checks the ref is
//                    one of them);
//   siblingsOf(ref)  the references at the same level, in order,
//                    including `ref` itself (optional);
//   same(a, b)       whether two references are the same object (default:
//                    ===).
// A reference is any plain object; the module only compares and hands
// them back. No DOM here: this file loads under node.

export function createSelection(resolver) {
  const same = resolver.same || function (a, b) { return a === b; };
  const listeners = [];
  let current = null;

  // Makes `ref` the selection and tells the listeners {ref, previous,
  // how}. `how` is "select", "up", "into" or "sideways"; `quiet` (a
  // display saying what it already shows) records the change without it.
  function move(ref, how, quiet) {
    if (ref == null) return null;
    const previous = current;
    if (previous != null && same(previous, ref)) {
      current = ref;
      return ref;
    }
    current = ref;
    if (!quiet) {
      listeners.slice().forEach(function (fn) { fn({ ref: ref, previous: previous, how: how }); });
    }
    return ref;
  }

  function siblings(ref) {
    return resolver.siblingsOf ? resolver.siblingsOf(ref) || [] : [];
  }

  return {
    current: function () { return current; },

    select: function (ref, options) {
      return move(ref, "select", options && options.quiet);
    },

    // The parent, or null (and no change) at the top.
    up: function () {
      if (current == null) return null;
      return move(resolver.parentOf(current), "up");
    },

    // `ref`, if it is one of the current selection's children; null (and
    // no change) otherwise.
    into: function (ref) {
      if (current == null || ref == null) return null;
      if (resolver.childrenOf) {
        const kids = resolver.childrenOf(current) || [];
        if (!kids.some(function (kid) { return same(kid, ref); })) return null;
      }
      return move(ref, "into");
    },

    // The next (dir > 0) or previous sibling, wrapping; null with fewer
    // than two.
    sideways: function (dir) {
      if (current == null) return null;
      const all = siblings(current);
      if (all.length < 2) return null;
      let at = all.findIndex(function (item) { return same(item, current); });
      if (at < 0) at = 0;
      const step = dir < 0 ? -1 : 1;
      return move(all[(at + step + all.length) % all.length], "sideways");
    },

    // The chain from the top down to the selection: a breadcrumb's steps.
    trail: function () {
      const chain = [];
      let ref = current;
      while (ref != null && chain.length < 64) {
        chain.unshift(ref);
        ref = resolver.parentOf(ref);
      }
      return chain;
    },

    // Calls `fn({ref, previous, how})` on every change; returns a function
    // that stops it.
    onChange: function (fn) {
      listeners.push(fn);
      return function () {
        const at = listeners.indexOf(fn);
        if (at >= 0) listeners.splice(at, 1);
      };
    },
  };
}
