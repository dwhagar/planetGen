// html/static/bookmarks.js
//
// Bookmarks (MAP.23, docs/design/galaxy-drilldown-navigation.md section
// 8.2): places saved per browser, in localStorage, one list per database
// ("planetgen.bookmarks.<db name>"), up to MAX_BOOKMARKS. Each entry is
// {name, kind, value, created}, plus the page it opens (url) and, for a
// generated sector, its id (sectorId, for the NAV page's system picker):
// - kind "stage": value is a Galaxy Map stage URL (path and query);
// - kind "sector": value is the sector's designation (or its page URL
//   when it has no galaxy address);
// - kind "system" or a phenomenon type ("nebula", ...): value is the NAV
//   endpoint ("system:12", "nebula:3", web/nav_page.endpoint).
// An entry is known by its kind and value. Every storage read and write
// is wrapped: with storage blocked (a private window, a strict setting)
// there are just no bookmarks, and saving says so.
//
// Importing this module wires whatever the page has (no inline script,
// see the Content-Security-Policy in web/__init__.py):
// - [data-bookmark-toggle]: a server-rendered ☆ Bookmark button on a
//   system, phenomenon or sector page, its entry in data-bookmark-kind,
//   -value, -name, -url and -sector-id. Hidden until wired.
// - [data-bookmarks-menu]: the Galaxy Map's Bookmarks menu (a <details>):
//   open, rename and delete; with data-bookmarks-keys="<map's id>", the
//   keys 1 to 9 open the first nine while that map (the element with that
//   id, the canvas or anything in it) has focus (MAP.81: Ctrl+1 to 9 are
//   the browser's own tab keys on Windows and Linux). The sector page
//   has one too while a NAV start or destination is being picked.
// - [data-bookmarks-nav]: the NAV page's Bookmarks select (a form), for
//   the endpoint in data-pick, keeping data-keep-name=data-keep-value;
//   [data-bookmarks-nav-group] around some is hidden while all are.
// While a NAV start or destination is being picked, a menu or select
// carries the pick (data-pick, data-keep-name, data-keep-value,
// data-nav-url, pickOf) and every bookmark keeps it (NAV.40, pickUrlOf).
// The database comes from the first element's data-bookmark-db. The
// Galaxy Map's breadcrumb ☆ (galaxystageview.js) uses toggleButton.

export const MAX_BOOKMARKS = 100;
const KEY_PREFIX = "planetgen.bookmarks.";
const CHANGE_EVENT = "planetgen-bookmarks-change";
const ENDPOINT_RE = /^[a-z_]+:\d+$/;

// --- Storage -------------------------------------------------------------------

function storageKey() {
  const el = document.querySelector("[data-bookmark-db]");
  const db = el ? el.getAttribute("data-bookmark-db") : "";
  return db ? KEY_PREFIX + db : null;
}

function isEntry(entry) {
  return !!entry && typeof entry === "object" && typeof entry.kind === "string" && typeof entry.value === "string"
    && typeof entry.name === "string" && entry.kind !== "" && entry.value !== "";
}

function read() {
  const key = storageKey();
  if (!key) return [];
  try {
    const parsed = JSON.parse(window.localStorage.getItem(key) || "[]");
    return Array.isArray(parsed) ? parsed.filter(isEntry).slice(0, MAX_BOOKMARKS) : [];
  } catch (e) {
    return [];
  }
}

// True when the list was saved.
function write(entries) {
  const key = storageKey();
  if (!key) return false;
  try {
    window.localStorage.setItem(key, JSON.stringify(entries));
  } catch (e) {
    return false;
  }
  window.dispatchEvent(new CustomEvent(CHANGE_EVENT));
  return true;
}

function indexOf(entries, kind, value) {
  return entries.findIndex(function (entry) { return entry.kind === kind && entry.value === value; });
}

// --- The API ---------------------------------------------------------------------

// Every bookmark, oldest first (the order the 1-9 keys and the menus use).
export function list() {
  return read();
}

export function find(kind, value) {
  const entries = read();
  const index = indexOf(entries, kind, value);
  return index < 0 ? null : entries[index];
}

// Saves `entry` ({name, kind, value, url?, sectorId?}): {ok, entry} or
// {ok: false, reason}. Already saved is ok, as it was.
export function add(entry) {
  if (!isEntry(entry)) return { ok: false, reason: "That can't be bookmarked." };
  const entries = read();
  const index = indexOf(entries, entry.kind, entry.value);
  if (index >= 0) return { ok: true, entry: entries[index] };
  if (entries.length >= MAX_BOOKMARKS) {
    return { ok: false, reason: "You have " + MAX_BOOKMARKS + " bookmarks, the most there can be; delete one first." };
  }
  const saved = { name: entry.name.trim() || entry.value, kind: entry.kind, value: entry.value, created: new Date().toISOString() };
  if (entry.url) saved.url = entry.url;
  if (entry.sectorId != null) saved.sectorId = entry.sectorId;
  entries.push(saved);
  if (!write(entries)) return { ok: false, reason: "This browser isn't letting the site save bookmarks." };
  return { ok: true, entry: saved };
}

export function rename(kind, value, name) {
  const entries = read();
  const index = indexOf(entries, kind, value);
  if (index < 0 || !name.trim()) return false;
  entries[index].name = name.trim();
  return write(entries);
}

export function remove(kind, value) {
  const entries = read();
  const index = indexOf(entries, kind, value);
  if (index < 0) return false;
  entries.splice(index, 1);
  return write(entries);
}

// Calls `fn` whenever the bookmarks change, here or in another tab.
export function onChange(fn) {
  window.addEventListener(CHANGE_EVENT, fn);
  window.addEventListener("storage", function (event) {
    if (event.key === null || event.key === storageKey()) fn();
  });
}

// The page a bookmark opens.
export function urlOf(entry) {
  if (entry.url) return entry.url;
  if (entry.kind === "stage") return entry.value;
  if (entry.kind === "sector") return entry.value.charAt(0) === "/" ? entry.value : "/galaxy?sector=" + encodeURIComponent(entry.value);
  return null;
}

// The NAV pick an element carries ({pick, keepName, keepValue, navUrl}:
// the endpoint being chosen, "from" or "to", and the other one already
// chosen), or null outside pick mode.
export function pickOf(element) {
  const pick = element ? element.getAttribute("data-pick") : null;
  if (pick !== "from" && pick !== "to") return null;
  return {
    pick: pick,
    keepName: element.getAttribute("data-keep-name") || (pick === "from" ? "to" : "from"),
    keepValue: element.getAttribute("data-keep-value") || "",
    navUrl: element.getAttribute("data-nav-url") || "/nav",
  };
}

// `url` (a path and query) with `params` ([[name, value]]) set in its
// query, replacing any of the same names.
function withParams(url, params) {
  const parsed = new URL(url, window.location.href);
  params.forEach(function (p) { parsed.searchParams.delete(p[0]); });
  params.forEach(function (p) { if (p[1]) parsed.searchParams.append(p[0], p[1]); });
  return parsed.pathname + parsed.search.replace(/%3A/gi, ":") + parsed.hash;
}

// The page a bookmark opens while `pick` (pickOf) is being chosen,
// keeping the other endpoint (NAV.40): a system or phenomenon is taken as
// that endpoint on the NAV page; a sector page, the Galaxy Map at a
// sector or a saved map view opens in pick mode. Without a pick, urlOf.
export function pickUrlOf(entry, pick) {
  if (!pick) return urlOf(entry);
  if (entry.kind !== "stage" && entry.kind !== "sector" && ENDPOINT_RE.test(entry.value)) {
    const ends = [[pick.pick, entry.value], [pick.keepName, pick.keepValue]];
    // The order nav_page.nav_url writes: from, then to.
    if (pick.pick === "to") ends.reverse();
    return withParams(pick.navUrl, ends);
  }
  const url = urlOf(entry);
  if (!url) return null;
  return withParams(url, [["pick", pick.pick], ["from", null], ["to", null], [pick.keepName, pick.keepValue]]);
}

export function kindLabel(kind) {
  if (kind === "stage") return "Map view";
  if (kind === "sector") return "Sector";
  if (kind === "system") return "System";
  return kind.replace(/_/g, " ").replace(/^./, function (c) { return c.toUpperCase(); });
}

// --- A ☆ button ----------------------------------------------------------------

// Makes `button` save, or offer to remove, the entry `entryFn()` returns
// (null: nothing to save here). `compact` shows just the star, with the
// words in its aria-label and title. Returns refresh(), for when the
// entry changes.
export function toggleButton(button, entryFn, compact) {
  function refresh() {
    const entry = entryFn();
    const saved = entry ? find(entry.kind, entry.value) : null;
    button.disabled = !entry || !storageKey();
    button.setAttribute("aria-pressed", saved ? "true" : "false");
    const words = saved ? "Bookmarked" : "Bookmark";
    const label = entry ? (saved ? "Bookmarked: " + saved.name + " (press to remove)" : "Bookmark " + entry.name) : "Bookmark";
    button.textContent = (saved ? "★" : "☆") + (compact ? "" : " " + words);
    button.setAttribute("aria-label", label);
    button.title = label;
  }
  button.addEventListener("click", function () {
    const entry = entryFn();
    if (!entry) return;
    const saved = find(entry.kind, entry.value);
    if (saved) {
      if (window.confirm("Remove the bookmark “" + saved.name + "”?")) remove(entry.kind, entry.value);
    } else {
      const result = add(entry);
      if (!result.ok) window.alert(result.reason);
    }
    refresh();
  });
  onChange(refresh);
  refresh();
  return refresh;
}

function wireToggles() {
  document.querySelectorAll("[data-bookmark-toggle]").forEach(function (button) {
    const d = button.dataset;
    const entry = { name: d.bookmarkName || "", kind: d.bookmarkKind, value: d.bookmarkValue, url: d.bookmarkUrl || null };
    if (d.bookmarkSectorId) entry.sectorId = Number(d.bookmarkSectorId);
    toggleButton(button, function () { return entry; }, false);
    button.hidden = false;
  });
}

// --- The Bookmarks menu ----------------------------------------------------------

function el(tag, className, text) {
  const node = document.createElement(tag);
  if (className) node.className = className;
  if (text != null) node.textContent = text;
  return node;
}

// Fills the menu's panel: each bookmark as a link (with its key, 1 to 9,
// for the first nine), Rename and Delete.
function renderMenu(menu) {
  const panel = menu.querySelector("[data-bookmarks-panel]");
  if (!panel) return;
  panel.textContent = "";
  const entries = read();
  if (!storageKey()) {
    panel.appendChild(el("p", "hint", "Bookmarks aren't available here."));
    return;
  }
  if (!entries.length) {
    panel.appendChild(el("p", "hint", "No bookmarks yet. Use ☆ on the map's breadcrumb, or on a system, phenomenon or sector page."));
    return;
  }
  const keys = menu.hasAttribute("data-bookmarks-keys");
  const pick = pickOf(menu);
  const ul = el("ul", "bookmarks-list");
  entries.forEach(function (entry, n) {
    const li = el("li", "bookmarks-item");
    const link = el("a", "bookmarks-open");
    const url = pickUrlOf(entry, pick);
    if (url) link.href = url;
    link.appendChild(el("span", "bookmarks-name", entry.name));
    const meta = el("span", "bookmarks-kind", kindLabel(entry.kind));
    if (keys && n < 9) meta.appendChild(el("kbd", null, String(n + 1)));
    link.appendChild(meta);
    li.appendChild(link);
    const actions = el("span", "bookmarks-actions");
    const renameButton = el("button", "bookmarks-action", "Rename");
    renameButton.type = "button";
    renameButton.setAttribute("aria-label", "Rename bookmark " + entry.name);
    renameButton.addEventListener("click", function () { startRename(menu, li, entry); });
    const deleteButton = el("button", "bookmarks-action", "Delete");
    deleteButton.type = "button";
    deleteButton.setAttribute("aria-label", "Delete bookmark " + entry.name);
    deleteButton.addEventListener("click", function () {
      remove(entry.kind, entry.value);
      // The change re-rendered the menu. Keep the keyboard in it: on the
      // entry that took this one's place, or the menu button when none
      // is left.
      renderMenu(menu);
      const links = panel.querySelectorAll(".bookmarks-open");
      (links[Math.min(n, links.length - 1)] || menu.querySelector("summary")).focus();
    });
    actions.appendChild(renameButton);
    actions.appendChild(deleteButton);
    li.appendChild(actions);
    ul.appendChild(li);
  });
  panel.appendChild(ul);
}

// Rename in place: a text box (Enter saves, Escape cancels).
function startRename(menu, li, entry) {
  li.textContent = "";
  const form = el("form", "bookmarks-rename");
  const label = el("label", "sr-only", "New name for " + entry.name);
  const input = el("input");
  input.type = "text";
  input.value = entry.name;
  input.maxLength = 120;
  input.id = "bookmarks-rename-" + Math.random().toString(36).slice(2);
  label.htmlFor = input.id;
  const save = el("button", "bookmarks-action", "Save");
  save.type = "submit";
  const cancel = el("button", "bookmarks-action", "Cancel");
  cancel.type = "button";
  function done() {
    renderMenu(menu);
    const links = menu.querySelectorAll(".bookmarks-open");
    const entries = read();
    const index = indexOf(entries, entry.kind, entry.value);
    if (links[index]) links[index].focus();
  }
  form.addEventListener("submit", function (event) {
    event.preventDefault();
    if (input.value.trim()) rename(entry.kind, entry.value, input.value);
    done();
  });
  cancel.addEventListener("click", done);
  input.addEventListener("keydown", function (event) {
    if (event.key === "Escape") {
      event.preventDefault();
      event.stopPropagation();
      done();
    }
  });
  form.appendChild(label);
  form.appendChild(input);
  form.appendChild(save);
  form.appendChild(cancel);
  li.appendChild(form);
  input.focus();
  input.select();
}

function isTextField(target) {
  if (!target || !target.tagName) return false;
  const tag = target.tagName.toLowerCase();
  return tag === "input" || tag === "textarea" || tag === "select" || target.isContentEditable;
}

function wireMenu(menu) {
  renderMenu(menu);
  onChange(function () {
    // Leave a rename in progress alone.
    if (!menu.querySelector(".bookmarks-rename")) renderMenu(menu);
  });
  menu.addEventListener("toggle", function () {
    if (menu.open) renderMenu(menu);
  });
  menu.addEventListener("keydown", function (event) {
    if (event.key === "Escape" && menu.open) {
      event.preventDefault();
      menu.open = false;
      menu.querySelector("summary").focus();
    }
  });
  // A click elsewhere closes it. The click's path, not its target: Rename
  // and Delete redraw the list, taking the clicked button out of it.
  document.addEventListener("click", function (event) {
    if (menu.open && event.composedPath().indexOf(menu) < 0) menu.open = false;
  });
  if (menu.hasAttribute("data-bookmarks-keys")) {
    const scope = document.getElementById(menu.getAttribute("data-bookmarks-keys")) || menu;
    document.addEventListener("keydown", function (event) {
      if (event.ctrlKey || event.altKey || event.metaKey || event.shiftKey || isTextField(event.target)) return;
      if (!event.target || !event.target.nodeType || !scope.contains(event.target)) return;
      const n = /^[1-9]$/.test(event.key) ? Number(event.key) : 0;
      const entry = n ? read()[n - 1] : null;
      const url = entry ? pickUrlOf(entry, pickOf(menu)) : null;
      if (!url) return;
      event.preventDefault();
      window.location.assign(url);
    });
  }
}

// --- The NAV page's Bookmarks select ----------------------------------------------

// The systems and phenomena (an endpoint each), the generated sectors
// (their system picker) and the saved map views (the Galaxy Map there,
// in pick mode), leaving out the endpoint already chosen.
function renderNavSelect(form) {
  const select = form.querySelector("select");
  if (!select) return;
  const keep = form.getAttribute("data-keep-value") || "";
  const entries = read();
  const places = entries.filter(function (e) { return e.kind !== "stage" && e.kind !== "sector" && ENDPOINT_RE.test(e.value) && e.value !== keep; });
  const sectors = entries.filter(function (e) { return e.kind === "sector" && Number.isInteger(e.sectorId); });
  const views = entries.filter(function (e) { return e.kind === "stage"; });
  select.textContent = "";
  [["Systems and phenomena", places, function (e) { return e.value; }],
    ["Sectors (then choose a system)", sectors, function (e) { return "sector:" + e.sectorId; }],
    ["Map views (then pick on the Galaxy Map)", views, function (e) { return "view:" + e.value; }]].forEach(function (group) {
    if (!group[1].length) return;
    const optgroup = el("optgroup");
    optgroup.label = group[0];
    group[1].forEach(function (entry) {
      const plain = entry.kind === "system" || entry.kind === "sector" || entry.kind === "stage";
      const option = el("option", null, entry.name + (plain ? "" : " (" + kindLabel(entry.kind) + ")"));
      option.value = group[2](entry);
      optgroup.appendChild(option);
    });
    select.appendChild(optgroup);
  });
  form.hidden = !(places.length || sectors.length || views.length);
  const group = form.closest("[data-bookmarks-nav-group]");
  if (group) group.hidden = !Array.prototype.some.call(group.querySelectorAll("[data-bookmarks-nav]"), function (f) { return !f.hidden; });
}

function wireNavSelect(form) {
  renderNavSelect(form);
  onChange(function () { renderNavSelect(form); });
  form.addEventListener("submit", function (event) {
    event.preventDefault();
    const chosen = form.querySelector("select").value;
    if (!chosen) return;
    if (chosen.indexOf("view:") === 0) {
      window.location.assign(pickUrlOf({ kind: "stage", value: chosen.slice("view:".length) }, pickOf(form)));
      return;
    }
    const pick = form.getAttribute("data-pick") === "to" ? "to" : "from";
    const params = [];
    const keepName = form.getAttribute("data-keep-name");
    const keepValue = form.getAttribute("data-keep-value");
    // The order nav_page.nav_url writes: from, then to.
    const pickParam = chosen.indexOf("sector:") === 0
      ? [pick + "_sector", chosen.slice("sector:".length)]
      : [pick, chosen];
    if (keepName && keepValue) params.push([keepName, keepValue]);
    params.push(pickParam);
    if (pick === "from") params.reverse();
    const query = params.map(function (p) { return p[0] + "=" + encodeURIComponent(p[1]).replace(/%3A/g, ":"); }).join("&");
    window.location.assign(form.getAttribute("data-nav-url") + "?" + query);
  });
}

wireToggles();
document.querySelectorAll("[data-bookmarks-menu]").forEach(wireMenu);
document.querySelectorAll("[data-bookmarks-nav]").forEach(wireNavSelect);
