// html/static/poiprefs.js
//
// MAP.306: the visitor's choice for each point-of-interest layer of the Galaxy Map's Menu (shown, pulled forward or hidden), remembered
// in this browser so the next map opens as the last one was left. The address still wins: a map opened with a layer named in it
// (`hide`, `show`, `front`) is as the address says, and what is remembered is only used when the address names none of them.
// No DOM and no three.js; `storage` is the browser's localStorage (or anything with getItem and setItem), which can be missing or throw,
// so every access is guarded. tests/js/poiprefs.test.mjs holds it.

export const STORAGE_KEY = "planetgen.map.poi";

const WORDS = ["hidden", "shown", "forward"];

// The remembered choices, `{layerId: "hidden" | "shown" | "forward"}`, or `{}` when there are none.
export function loadPrefs(storage) {
  try {
    const parsed = JSON.parse(storage.getItem(STORAGE_KEY) || "{}");
    const found = {};
    Object.keys(parsed || {}).forEach(function (id) { if (WORDS.indexOf(parsed[id]) >= 0) found[id] = parsed[id]; });
    return found;
  } catch (error) {
    return {};
  }
}

// Remembers the state of each of `ids` in `layers` (a layer registry).
export function savePrefs(storage, layers, ids) {
  const now = {};
  ids.forEach(function (id) { now[id] = layers.stateOf(id); });
  try {
    storage.setItem(STORAGE_KEY, JSON.stringify(now));
  } catch (error) {
    // Private windows and full stores refuse; the choice then lasts until the page is closed.
  }
  return now;
}

// Opens `ids` as they were remembered, unless the address already chose for one of them (`defaults`: each id's state word with
// nothing chosen). Returns true when the remembered choices were used.
export function applyPrefs(storage, layers, ids, defaults) {
  if (ids.some(function (id) { return layers.stateOf(id) !== defaults[id]; })) return false;
  const saved = loadPrefs(storage);
  let used = false;
  ids.forEach(function (id) {
    const word = saved[id];
    if (!word || word === layers.stateOf(id)) return;
    used = true;
    if (word === "forward") layers.setPulled(id, true);
    else layers.setOn(id, word === "shown");
  });
  return used;
}
