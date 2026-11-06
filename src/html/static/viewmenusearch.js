// html/static/viewmenusearch.js
//
// UX.279: the search over the View menu's entries (design: design-drafts/view-menu-many-options-ux274-2026-11-05.md, section 4). A pure
// function, no DOM: it takes the entries `registry.entries()` (maplayers.js, UX.277) hands out and the text a visitor typed, and says
// which entries match and which headings and groups the menu should open to show them. The box that uses it is UX.281; nothing shows
// it yet. tests/js/viewmenusearch.test.mjs holds it.

// The words of a text, lower case, split at anything that is not a letter or a digit: "World classes" -> ["world", "classes"].
export function wordsOf(text) {
  return String(text == null ? "" : text).toLowerCase().split(/[^\p{L}\p{N}]+/u).filter(Boolean);
}

// True when every word the visitor typed starts a word of the entry's label, its group's name or one of its keywords ("bla" finds
// "Black holes"; "world cl" finds the world classes; a word found nowhere finds nothing).
export function entryMatches(entry, query) {
  const wanted = wordsOf(query);
  if (!wanted.length) return true;
  const haystack = wordsOf([entry.label, entry.group].concat(entry.keywords || []).join(" "));
  return wanted.every(function (word) {
    return haystack.some(function (have) { return have.startsWith(word); });
  });
}

// The search: {query, active, matches, ids, headings, groups, count, text}.
//   active    false for a blank query (the menu shows everything, as the visitor left it);
//   matches   the matching entries in the menu's order;
//   ids       their ids;
//   headings  the headings that have a match, in the order they first appear, for the menu to open;
//   groups    [{heading, group}] the groups that have a match, for the menu to open (a group with no name is left out);
//   count     how many entries match;
//   text      what a screen reader is told: "5 matches", "1 match" or "No matches".
export function searchEntries(entries, query) {
  const active = wordsOf(query).length > 0;
  const matches = active ? entries.filter(function (entry) { return entryMatches(entry, query); }) : entries.slice();
  const headings = [];
  const groups = [];
  if (active) {
    matches.forEach(function (entry) {
      if (headings.indexOf(entry.heading) < 0) headings.push(entry.heading);
      if (entry.group && !groups.some(function (g) { return g.heading === entry.heading && g.group === entry.group; })) {
        groups.push({ heading: entry.heading, group: entry.group });
      }
    });
  }
  const count = matches.length;
  return {
    query: String(query == null ? "" : query),
    active: active,
    matches: matches,
    ids: matches.map(function (entry) { return entry.id; }),
    headings: headings,
    groups: groups,
    count: count,
    text: !active ? "" : count === 0 ? "No matches" : count === 1 ? "1 match" : count + " matches",
  };
}
