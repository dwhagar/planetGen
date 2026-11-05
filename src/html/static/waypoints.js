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
// An entry is {name, kind, value, url?, position?, created}: `value` is the NAV endpoint ("system:12", "nebula:3"), the entry's key, and
// `position` (x, y, z in parsecs, galaxy frame) is where the Galaxy Map flags it when the card it was picked from knew it.
//
// Importing this module wires whatever the page has (no inline script, see the Content-Security-Policy in web/__init__.py):
// - [data-waypoint-toggle]: an Add waypoint button on a system, phenomenon or body page, its entry in data-waypoint-kind, -value, -name
//   and -url. Hidden until wired.
// - [data-waypoints-menu]: the Galaxy Map's Waypoints menu (a <details>), its [data-waypoints-panel] filled with the list and its
//   [data-waypoints-count] with the number of waypoints. Two or more waypoints plot a course (NAV.49): the Galaxy Map takes them in its
//   address, `?via=<a>,<b>,...`, and draws the course through them (web/nav_page.galaxy_course_through). The course stays: a Galaxy Map
//   opened with two or more waypoints and no course of its own is opened again with them plotted, until they are cleared (Clear all, or
//   Remove down to fewer than two, which opens the map again without the course).
// The database comes from the first element's data-bookmark-db, as the bookmarks' does.

const VERSION_QUERY = new URL(import.meta.url).search;
const { createStore, localBackend } = await import(`./savedlist.js${VERSION_QUERY}`);
const SC = await import(`./savedcourse.js${VERSION_QUERY}`);

export const MAX_WAYPOINTS = 20;
const FLAG = "\u2691";
const KEY_PREFIX = "planetgen.waypoints.";
// `<kind>:<printed object ID>` (objectref.js; API.23): hex parts, up to three, joined by dashes.
const ENDPOINT_RE = /^[a-z_]+:[0-9A-Fa-f]+(?:-[0-9A-Fa-f]+){0,2}$/;

// --- The store -------------------------------------------------------------------

function database() {
  const el = document.querySelector("[data-bookmark-db], [data-course-db]");
  return el ? el.getAttribute("data-bookmark-db") || el.getAttribute("data-course-db") || "" : "";
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
  if (Array.isArray(entry.position) && entry.position.length === 3 && entry.position.every(Number.isFinite)) saved.position = entry.position.slice();
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

// A promise: every waypoint removed, and the saved course they were opened from forgotten.
export async function clear() {
  setSource(null);
  await store.setAll([]);
}

// A promise of whether the waypoint with endpoint `value` now stands `by` places later (negative: earlier) in the course (NAV.58).
export function move(value, by) {
  return store.move(value, by);
}

// A promise of {ok, entry}: the waypoint with endpoint `value` is now `entry` ({name, value, kind?, url?, position?}), in the same place
// in the course (NAV.58: a waypoint dragged onto another object).
export function replace(value, entry) {
  return store.replace(value, entry);
}

// --- The saved course the waypoints came from (NAV.58) ------------------------------------

const SOURCE_PREFIX = "planetgen.waypoint-course.";

function sourceKeyName() {
  const name = database();
  return name ? SOURCE_PREFIX + name : null;
}

// The key of the saved course (savedcourse.js) the waypoints were opened from or last saved as, or null: Save over replaces that one.
export function source() {
  const key = sourceKeyName();
  if (!key) return null;
  try {
    return window.localStorage.getItem(key) || null;
  } catch (err) {
    return null;
  }
}

function setSource(courseKey) {
  const key = sourceKeyName();
  if (!key) return;
  try {
    if (courseKey) window.localStorage.setItem(key, courseKey);
    else window.localStorage.removeItem(key);
  } catch (err) {
    // Not kept: Save over is just not offered.
  }
}

export const courses = SC.courseStore(database);

// The saved courses that go through waypoints, in the order saved.
export function savedCourses() {
  return courses.list().filter(SC.isWaypointCourse);
}

// Saves the waypoints as a course named `name`, over the saved course `overKey` when given: a promise of {ok, record} or {ok: false,
// reason}. The saved course becomes the one the waypoints came from.
export async function saveCourse(name, overKey) {
  const entries = list();
  if (entries.length < 2) return { ok: false, reason: "A course needs two or more waypoints." };
  const record = SC.recordOfWaypoints(entries, name);
  const result = overKey ? await courses.replace(overKey, record) : await courses.add(record);
  if (!result.ok) return result;
  setSource(SC.courseKey(result.entry));
  return { ok: true, record: result.entry };
}

// Makes the waypoints those of the saved course `record` (one write), remembering it as their source: a promise of {ok}.
export async function openCourse(record) {
  const result = await store.setAll(record.stops);
  if (result.ok) setSource(SC.courseKey(record));
  return result;
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

function wireToggles() {
  document.querySelectorAll("[data-waypoint-toggle]").forEach(function (button) {
    const d = button.dataset;
    const entry = { name: d.waypointName || "", kind: d.waypointKind, value: d.waypointValue, url: d.waypointUrl || null };
    toggleButton(button, function () { return entry; });
  });
}

// --- The course through the waypoints (NAV.49) -------------------------------------------

// The endpoints a Galaxy Map address plots (`?via=`), in order.
export function viaOf(search) {
  const value = new URLSearchParams(search || "").get("via") || "";
  return value.split(",").map(function (part) { return part.trim(); }).filter(Boolean);
}

export function sameList(a, b) {
  return a.length === b.length && a.every(function (value, n) { return value === b[n]; });
}

// The address `href` (a URL) with the course through `values` plotted (`via`), or, for fewer than two, with no course through
// waypoints; the rest of the address is kept (the view, the layers). An explicit `course` goes: it is the other way to ask for one.
export function plotUrl(href, values) {
  const parsed = new URL(href, "http://localhost");
  parsed.searchParams.delete("course");
  if (values.length >= 2) parsed.searchParams.set("via", values.join(","));
  else parsed.searchParams.delete("via");
  return parsed.pathname + parsed.search.replace(/%3A/gi, ":").replace(/%2C/gi, ",") + parsed.hash;
}

// What the address says about the course at this visit: "none", "other" (a course of its own, or a pick or a generate run on the
// map: the waypoints leave the map alone), "plotted" (the course through exactly these waypoints), or "stale" (through others).
export function plotState(search, values) {
  const params = new URLSearchParams(search || "");
  if (params.get("course") || params.get("pick") || params.get("made")) return "other";
  const via = viaOf(search);
  if (!via.length) return "none";
  return sameList(via, values) ? "plotted" : "stale";
}

function here() {
  return window.location.pathname + window.location.search + window.location.hash;
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
    if (courses.available()) renderSaving(menu, panel, entries);
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
    [["Up", -1, "earlier"], ["Down", 1, "later"]].forEach(function (way) {
      const step = el("button", "bookmarks-action waypoints-move", way[0]);
      step.type = "button";
      step.setAttribute("aria-label", "Move waypoint " + (n + 1) + ", " + entry.name + ", " + way[2]);
      step.disabled = way[1] < 0 ? n === 0 : n === entries.length - 1;
      step.addEventListener("click", async function () {
        await move(entry.value, way[1]);
        renderMenu(menu);
        const items = panel.querySelectorAll("li.bookmarks-item");
        const again = items[n + way[1]];
        const same = again && again.querySelector(".waypoints-move:nth-of-type(" + (way[1] < 0 ? 1 : 2) + ")");
        const target = same && !same.disabled ? same : again && again.querySelector(".waypoints-move:not(:disabled)");
        (target || panel.querySelector(".bookmarks-open") || menu.querySelector("summary")).focus();
      });
      actions.appendChild(step);
    });
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
  const values = entries.map(function (entry) { return entry.value; });
  const state = plotState(window.location.search, values);
  if (values.length >= 2 && state !== "other") {
    if (state === "plotted") {
      panel.appendChild(el("p", "hint", "The course through these " + values.length + " waypoints is drawn on the map."));
    } else {
      const plot = el("a", "bookmarks-action waypoints-plot", state === "stale" ? "Plot the course again" : "Plot the course");
      plot.href = plotUrl(here(), values);
      panel.appendChild(plot);
    }
  } else if (values.length < 2 && state === "plotted") {
    panel.appendChild(el("p", "hint", "A course needs two or more waypoints."));
  }
  if (courses.available()) renderSaving(menu, panel, entries);
  const clearButton = el("button", "bookmarks-action waypoints-clear", "Clear all");
  clearButton.type = "button";
  clearButton.addEventListener("click", async function () {
    await clear();
    if (viaOf(window.location.search).length) {
      window.location.assign(plotUrl(here(), []));   // the course through them goes with them
      return;
    }
    renderMenu(menu);
    menu.querySelector("summary").focus();
  });
  panel.appendChild(clearButton);
}

// Saving the waypoints as a course, and the saved courses through waypoints (NAV.58): Save over the one these waypoints came from, Save as
// new, a line saying the course is saved, and Open for each saved course, which makes its stops the waypoints and plots them.
function renderSaving(menu, panel, entries) {
  const box = el("div", "waypoints-saving");
  const stopsKey = entries.length >= 2 ? SC.courseKey(SC.recordOfWaypoints(entries, "")) : null;
  const saved = stopsKey && courses.list().find(function (other) { return SC.courseKey(other) === stopsKey; });
  const from = courses.list().find(function (other) { return SC.courseKey(other) === source(); });
  if (saved) box.appendChild(el("p", "hint", "Saved as “" + saved.name + "”."));
  function done(result, button) {
    if (!result.ok) window.alert(result.reason);
    renderMenu(menu);
    (button.isConnected ? button : menu.querySelector("summary")).focus();
  }
  if (entries.length >= 2 && from && !saved) {
    const over = el("button", "bookmarks-action waypoints-save-over", "Save over “" + from.name + "”");
    over.type = "button";
    over.addEventListener("click", async function () { done(await saveCourse(from.name, SC.courseKey(from)), over); });
    box.appendChild(over);
  }
  if (entries.length >= 2 && !saved) {
    const fresh = el("button", "bookmarks-action waypoints-save-new", "Save as new course");
    fresh.type = "button";
    fresh.addEventListener("click", async function () {
      const name = window.prompt("Name this course:", entries[0].name + " to " + entries[entries.length - 1].name);
      if (name === null) return;
      done(await saveCourse(name, null), fresh);
    });
    box.appendChild(fresh);
  }
  const others = savedCourses().filter(function (record) { return SC.courseKey(record) !== stopsKey; });
  if (others.length) {
    box.appendChild(el("p", "hint", "Saved courses through waypoints:"));
    const ul = el("ul", "bookmarks-list");
    others.forEach(function (record) {
      const li = el("li", "bookmarks-item");
      li.appendChild(el("span", "bookmarks-name", record.name));
      const open = el("button", "bookmarks-action", "Open");
      open.type = "button";
      open.setAttribute("aria-label", "Open saved course " + record.name + " as the waypoints");
      open.addEventListener("click", async function () { await openAndPlot(record); });
      li.appendChild(open);
      ul.appendChild(li);
    });
    box.appendChild(ul);
  }
  panel.appendChild(box);
}

// Makes the stops of the saved course `record` the waypoints and shows the course through them on the Galaxy Map, opening that map when
// this is not it (`href` is the map's address; default: this page's own, a Galaxy Map).
export async function openAndPlot(record, href) {
  const result = await openCourse(record);
  if (!result.ok) {
    window.alert(result.reason);
    return;
  }
  rememberOpen();
  window.location.assign(plotUrl(href || here(), list().map(function (entry) { return entry.value; })));
}

// The course follows an edit (NAV.58): when the map shows the course through the waypoints and they change, the address is changed to plot
// the course through the new list, or to drop it when fewer than two are left; the menu comes back open as it was.
const OPEN_KEY = "planetgen.waypoints-open";

function rememberOpen() {
  try {
    window.sessionStorage.setItem(OPEN_KEY, "1");
  } catch (err) {
    // The menu just comes back closed.
  }
}

function wasOpen() {
  try {
    const open = window.sessionStorage.getItem(OPEN_KEY) === "1";
    window.sessionStorage.removeItem(OPEN_KEY);
    return open;
  } catch (err) {
    return false;
  }
}

function followEdit(menus) {
  const values = list().map(function (entry) { return entry.value; });
  const shown = viaOf(window.location.search).length > 0;
  const state = plotState(window.location.search, values);
  if (!shown || state === "other" || state === "plotted") return;
  if (menus.some(function (menu) { return menu.open; })) rememberOpen();
  window.location.assign(plotUrl(here(), values));
}

// The course stays until the waypoints are cleared: a Galaxy Map with two or more and no course through them is opened again with them
// plotted, and one showing a course through waypoints that are no longer two or more is opened again without it.
function keepCourse() {
  const values = list().map(function (entry) { return entry.value; });
  const state = plotState(window.location.search, values);
  if (state === "none" && values.length >= 2) window.location.replace(plotUrl(here(), values));
  else if (viaOf(window.location.search).length && values.length < 2 && state !== "other") window.location.replace(plotUrl(here(), []));
}

function wireMenus() {
  const menus = document.querySelectorAll("[data-waypoints-menu]");
  menus.forEach(function (menu) {
    renderMenu(menu);
    onChange(function () {
      renderMenu(menu);
      followEdit(Array.from(menus));
    });
    courses.onChange(function () { renderMenu(menu); });
    ready.then(function () { renderMenu(menu); });
    if (wasOpen()) menu.open = true;
  });
  if (menus.length) ready.then(keepCourse);
}

wireMenus();
wireToggles();
