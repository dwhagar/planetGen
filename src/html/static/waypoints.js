// html/static/waypoints.js
//
// Waypoints (NAV.49, docs/design/orbits/course-routing.md section 10): objects picked on the maps as the stops of a course. A star
// system or a phenomenon (anything with a NAV endpoint) is added from its info panel with the Waypoint button and removed the same
// way; the Waypoints menu lists them in the order they were picked, with Remove and Clear all. They are kept per browser, in
// localStorage, one list per database ("planetgen.waypoints.<db name>"), through the saved-list store (savedlist.js) with its local
// backend, so they stay until cleared and a reload or another tab sees the same list. The writes answer with a promise; the reads
// (list, find) are at once. Every storage read and write is wrapped by the store: with storage blocked there are just no waypoints, and
// adding one says so.
//
// An entry is {name, kind, value, url?, created}: `value` is the NAV endpoint ("system:12", "nebula:3"), the entry's key.
//
// Importing this module wires whatever the page has (no inline script, see the Content-Security-Policy in web/__init__.py):
// - [data-waypoints-menu]: the Galaxy Map's Waypoints menu (a <details>), its [data-waypoints-panel] filled with the list and its
//   [data-waypoints-count] with the number of waypoints.
// The database comes from the first element's data-bookmark-db, as the bookmarks' does.

const VERSION_QUERY = new URL(import.meta.url).search;
const { createStore, localBackend } = await import(`./savedlist.js${VERSION_QUERY}`);

export const MAX_WAYPOINTS = 20;
const FLAG = "\u2691";
const KEY_PREFIX = "planetgen.waypoints.";
// `<kind>:<printed object ID>` (objectref.js; API.23): hex parts, up to three, joined by dashes.
const ENDPOINT_RE = /^[a-z_]+:[0-9A-Fa-f]+(?:-[0-9A-Fa-f]+){0,2}$/;

// --- The store -------------------------------------------------------------------

function database() {
  const el = document.querySelector("[data-bookmark-db]");
  return el ? el.getAttribute("data-bookmark-db") || "" : "";
}

// True for something a course can go through: an entry with a NAV endpoint (a map view or a sector, which the bookmarks keep, has none).
export function isWaypoint(entry) {
  return !!entry && typeof entry === "object" && typeof entry.value === "string" && ENDPOINT_RE.test(entry.value)
    && entry.kind !== "stage" && entry.kind !== "sector" && typeof entry.name === "string";
}

function keyOf(entry) {
  return isWaypoint(entry) ? entry.value : null;
}

function build(entry) {
  const saved = { name: entry.name.trim() || entry.value, kind: entry.value.split(":")[0], value: entry.value, created: new Date().toISOString() };
  if (entry.url) saved.url = entry.url;
  return saved;
}

export const store = createStore({
  backend: localBackend({ prefix: KEY_PREFIX, suffix: database }),
  keyOf: keyOf,
  build: build,
  max: MAX_WAYPOINTS,
  words: {
    invalid: "That can't be a waypoint.",
    full: function (most) { return "A course has " + most + " waypoints at most; remove one first."; },
    blocked: "This browser isn't letting the site save waypoints.",
  },
});

// --- The API ---------------------------------------------------------------------

export const ready = store.ready;

// Every waypoint, in the order they were picked.
export function list() {
  return store.list();
}

export function find(value) {
  return store.find(value);
}

// 1-based place of the waypoint with endpoint `value` in the list, or 0.
export function placeOf(value) {
  return store.list().findIndex(function (entry) { return entry.value === value; }) + 1;
}

// Adds `entry` ({name, value, kind?, url?}) after the last: a promise of {ok, entry} or {ok: false, reason}. Already a waypoint is ok.
export function add(entry) {
  return store.add(entry);
}

// A promise of whether the waypoint is gone.
export function remove(value) {
  return store.remove(value);
}

// A promise: every waypoint removed.
export async function clear() {
  for (const entry of store.list()) await store.remove(entry.value);
}

// Calls `fn` whenever the waypoints change, here or in another tab.
export function onChange(fn) {
  store.onChange(fn);
}

// The words on a waypoint button: add it, or remove it with its place in the course.
export function buttonWords(entry, place) {
  if (!entry) return "Waypoint";
  return place > 0 ? "Remove waypoint " + place : "Add waypoint";
}

// --- A Waypoint button -------------------------------------------------------------

// Makes `button` add, or remove, the waypoint `entryFn()` returns (null: nothing to add here). Returns refresh(), for when the entry
// changes. The button is hidden when the entry is not a place a course can go through. With `compact` it shows a flag (and the place,
// once a waypoint) with the words in its label and title, to sit by an object's name as the bookmark's star does.
export function toggleButton(button, entryFn, compact) {
  function refresh() {
    const entry = entryFn();
    const usable = isWaypoint(entry);
    button.hidden = !usable;
    if (!usable) return;
    const place = placeOf(entry.value);
    button.disabled = !store.available();
    button.setAttribute("aria-pressed", place > 0 ? "true" : "false");
    const words = buttonWords(entry, place);
    const title = place > 0 ? "Take " + entry.name + " out of the waypoints of the course" : "Make " + entry.name + " a stop of a course";
    button.textContent = compact ? FLAG + (place > 0 ? " " + place : "") : words;
    button.setAttribute("aria-label", compact ? words + ": " + entry.name : words);
    button.title = title;
  }
  button.addEventListener("click", async function () {
    const entry = entryFn();
    if (!isWaypoint(entry)) return;
    if (find(entry.value)) {
      await remove(entry.value);
    } else {
      const result = await add(entry);
      if (!result.ok) window.alert(result.reason);
    }
    refresh();
  });
  onChange(refresh);
  refresh();
  ready.then(refresh);
  return refresh;
}

// --- The Waypoints menu --------------------------------------------------------------

function el(tag, className, text) {
  const node = document.createElement(tag);
  if (className) node.className = className;
  if (text != null) node.textContent = text;
  return node;
}

function kindLabel(kind) {
  if (kind === "system") return "System";
  return kind.replace(/_/g, " ").replace(/^./, function (c) { return c.toUpperCase(); });
}

function renderMenu(menu) {
  const panel = menu.querySelector("[data-waypoints-panel]");
  const count = menu.querySelector("[data-waypoints-count]");
  const entries = list();
  if (count) count.textContent = entries.length ? " (" + entries.length + ")" : "";
  if (!panel) return;
  panel.textContent = "";
  if (!store.available()) {
    panel.appendChild(el("p", "hint", "Waypoints aren't available here."));
    return;
  }
  if (!entries.length) {
    panel.appendChild(el("p", "hint", "No waypoints yet. Select a star system or object and press Add waypoint."));
    return;
  }
  const ol = el("ol", "bookmarks-list waypoints-list");
  entries.forEach(function (entry, n) {
    const li = el("li", "bookmarks-item");
    const link = el("a", "bookmarks-open");
    if (entry.url) link.href = entry.url;
    link.appendChild(el("span", "bookmarks-name", (n + 1) + ". " + entry.name));
    link.appendChild(el("span", "bookmarks-kind", kindLabel(entry.kind)));
    li.appendChild(link);
    const actions = el("span", "bookmarks-actions");
    const removeButton = el("button", "bookmarks-action", "Remove");
    removeButton.type = "button";
    removeButton.setAttribute("aria-label", "Remove waypoint " + (n + 1) + ", " + entry.name);
    removeButton.addEventListener("click", async function () {
      await remove(entry.value);
      renderMenu(menu);
      const links = panel.querySelectorAll(".bookmarks-open");
      (links[Math.min(n, links.length - 1)] || menu.querySelector("summary")).focus();
    });
    actions.appendChild(removeButton);
    li.appendChild(actions);
    ol.appendChild(li);
  });
  panel.appendChild(ol);
  const clearButton = el("button", "bookmarks-action waypoints-clear", "Clear all");
  clearButton.type = "button";
  clearButton.addEventListener("click", async function () {
    await clear();
    renderMenu(menu);
    menu.querySelector("summary").focus();
  });
  panel.appendChild(clearButton);
}

function wireMenus() {
  document.querySelectorAll("[data-waypoints-menu]").forEach(function (menu) {
    renderMenu(menu);
    onChange(function () { renderMenu(menu); });
    ready.then(function () { renderMenu(menu); });
  });
}

wireMenus();
