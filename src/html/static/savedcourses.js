// html/static/savedcourses.js
//
// Saved courses on the NAV page (NAV.4, NAV.18, docs/TODO.md): the Save course button on a course, what has changed since it was
// saved, and the Courses list with open, rename, delete and the Direct / Route switch. The record and the store are savedcourse.js
// (NAV.17): one list per database in the browser ("planetgen.courses.<db name>"), through the saved-list store (savedlist.js), so a
// course is kept in this browser until accounts exist (NAV.19).
//
// Importing this module wires whatever the page has (no inline script, see the Content-Security-Policy in web/__init__.py):
// - [data-course-save]: on a course's page, a wrapper hidden until wired. Its data-course-db names the database, its
//   data-course-result holds what `/api/nav` answers for the course (the ends, `direct`, `route`, `adjusted`). It holds the
//   [data-course-save-button], the [data-course-status] line and the [data-course-changes] list.
// - [data-courses-list]: the Courses panel, hidden until wired; its [data-courses-panel] is filled with the saved courses.
//
// Every storage read and write is wrapped by the store: with storage blocked there are just no saved courses, and saving says so.

const VERSION_QUERY = new URL(import.meta.url).search;
const { courseStore, courseKey, recordOf, changesSince } = await import(`./savedcourse.js${VERSION_QUERY}`);
const { formatDistanceLy } = await import(`./unitladder.js${VERSION_QUERY}`);

// --- Words and links (pure) -------------------------------------------------------

export const FORMS = ["direct", "route"];

function stops(count) {
  return count + (count === 1 ? " stop" : " stops");
}

// The line a saved course shows for one of its forms, "direct" or "route" (a course with no route says so).
export function describeCourse(record, form) {
  if (form === "route") {
    if (!record.route) return "No route through known space.";
    const between = Math.max(record.route.path.length - 2, 0);
    return formatDistanceLy(record.route.distance_ly) + ", " + stops(between) + " between the ends";
  }
  const direct = record.direct;
  if (!direct) return "No direct course.";
  const pad = function (degrees) { return String(Math.round(((degrees % 360) + 360) % 360) % 360).padStart(3, "0"); };
  return formatDistanceLy(direct.distance_ly) + ", " + pad(direct.bearing_deg) + " mark " + pad(direct.mark_deg);
}

// The NAV page of a saved course: the course is planned again there and says what has changed.
export function courseUrl(record) {
  return "/nav?from=" + encodeURIComponent(record.from).replace(/%3A/g, ":") + "&to=" + encodeURIComponent(record.to).replace(/%3A/g, ":");
}

// The Galaxy Map showing the course.
export function mapUrl(record) {
  return "/galaxy?course=" + encodeURIComponent(record.from).replace(/%3A/g, ":") + "," + encodeURIComponent(record.to).replace(/%3A/g, ":");
}

// --- The store ----------------------------------------------------------------------

function database() {
  const el = document.querySelector("[data-course-db]");
  return el ? el.getAttribute("data-course-db") || "" : "";
}

export const store = courseStore(database);

function keyOfRecord(record) {
  return courseKey(record);
}

function el(tag, className, text) {
  const node = document.createElement(tag);
  if (className) node.className = className;
  if (text != null) node.textContent = text;
  return node;
}

// --- The Save course button ---------------------------------------------------------

function readResult(wrapper) {
  try {
    const result = JSON.parse(wrapper.getAttribute("data-course-result") || "null");
    return result && result.origin && result.destination ? result : null;
  } catch (err) {
    return null;
  }
}

function wireSave(wrapper) {
  const result = readResult(wrapper);
  const button = wrapper.querySelector("[data-course-save-button]");
  if (!result || !button) return;
  const status = document.querySelector("[data-course-status]");
  const changes = document.querySelector("[data-course-changes]");
  const key = result.origin.ref + ">" + result.destination.ref;

  function refresh() {
    const saved = store.find(key);
    button.disabled = !store.available();
    button.textContent = saved ? "Update saved course" : "Save course";
    if (status) {
      status.hidden = !saved;
      status.textContent = saved ? "Saved as “" + saved.name + "”." : "";
    }
    if (changes) {
      changes.textContent = "";
      const found = saved ? changesSince(saved, result) : [];
      found.forEach(function (sentence) { changes.appendChild(el("li", null, sentence)); });
      if (saved && !found.length) changes.appendChild(el("li", null, "Nothing has changed since you saved it."));
      changes.hidden = !saved;
    }
  }

  button.addEventListener("click", async function () {
    const saved = store.find(key);
    if (saved) {
      // Plan it again: the same name, today's course.
      await store.remove(key);
      const again = await store.add(recordOf(result, saved.name));
      if (!again.ok) window.alert(again.reason);
    } else {
      const fallback = (result.origin.name || result.origin.ref) + " to " + (result.destination.name || result.destination.ref);
      const name = window.prompt("Name this course:", fallback);
      if (name === null) return;
      const added = await store.add(recordOf(result, name));
      if (!added.ok) window.alert(added.reason);
    }
    refresh();
  });
  store.onChange(refresh);
  refresh();
  store.ready.then(refresh);
  wrapper.hidden = false;
}

// --- The Courses list ---------------------------------------------------------------

const shown = new Map();      // the form each course shows, by key, for this visit

function renderList(panel) {
  const target = panel.querySelector("[data-courses-panel]");
  if (!target) return;
  target.textContent = "";
  const entries = store.list();
  if (!store.available()) {
    target.appendChild(el("p", "hint", "Saved courses aren't available here."));
    return;
  }
  if (!entries.length) {
    target.appendChild(el("p", "hint", "No saved courses yet. Use Save course on a course."));
    return;
  }
  const ul = el("ul", "bookmarks-list");
  entries.forEach(function (entry, n) {
    const key = keyOfRecord(entry);
    const form = shown.get(key) || "direct";
    const li = el("li", "bookmarks-item");
    const link = el("a", "bookmarks-open");
    link.href = courseUrl(entry);
    link.appendChild(el("span", "bookmarks-name", entry.name));
    link.appendChild(el("span", "bookmarks-kind", entry.fromName + " → " + entry.toName));
    link.appendChild(el("span", "bookmarks-kind courses-form-text", describeCourse(entry, form)));
    li.appendChild(link);
    const actions = el("span", "bookmarks-actions");
    FORMS.forEach(function (name) {
      const choose = el("button", "bookmarks-action", name === "direct" ? "Direct" : "Route");
      choose.type = "button";
      choose.setAttribute("aria-pressed", name === form ? "true" : "false");
      choose.setAttribute("aria-label", "Show the " + name + " form of " + entry.name);
      choose.addEventListener("click", function () {
        shown.set(key, name);
        renderList(panel);
        const again = panel.querySelectorAll("li.bookmarks-item")[n];
        if (again) again.querySelectorAll("button[aria-pressed]")[FORMS.indexOf(name)].focus();
      });
      actions.appendChild(choose);
    });
    const map = el("a", "bookmarks-action", "On the map");
    map.href = mapUrl(entry);
    map.setAttribute("aria-label", "Show " + entry.name + " on the Galaxy Map");
    actions.appendChild(map);
    const renameButton = el("button", "bookmarks-action", "Rename");
    renameButton.type = "button";
    renameButton.setAttribute("aria-label", "Rename saved course " + entry.name);
    renameButton.addEventListener("click", function () { startRename(panel, li, entry); });
    const deleteButton = el("button", "bookmarks-action", "Delete");
    deleteButton.type = "button";
    deleteButton.setAttribute("aria-label", "Delete saved course " + entry.name);
    deleteButton.addEventListener("click", async function () {
      await store.remove(key);
      renderList(panel);
      const links = target.querySelectorAll(".bookmarks-open");
      (links[Math.min(n, links.length - 1)] || panel.querySelector("h2")).focus();
    });
    actions.appendChild(renameButton);
    actions.appendChild(deleteButton);
    li.appendChild(actions);
    ul.appendChild(li);
  });
  target.appendChild(ul);
}

// Rename in place: a text box (Enter saves, Escape cancels).
function startRename(panel, li, entry) {
  li.textContent = "";
  const form = el("form", "bookmarks-rename");
  const label = el("label", "sr-only", "New name for " + entry.name);
  const input = el("input");
  input.type = "text";
  input.value = entry.name;
  input.maxLength = 120;
  input.id = "courses-rename-" + Math.random().toString(36).slice(2);
  label.htmlFor = input.id;
  const save = el("button", "bookmarks-action", "Save");
  save.type = "submit";
  const cancel = el("button", "bookmarks-action", "Cancel");
  cancel.type = "button";
  function done() {
    renderList(panel);
    const links = panel.querySelectorAll(".bookmarks-open");
    const index = store.list().findIndex(function (other) { return keyOfRecord(other) === keyOfRecord(entry); });
    if (links[index]) links[index].focus();
  }
  form.addEventListener("submit", async function (event) {
    event.preventDefault();
    if (input.value.trim()) await store.rename(keyOfRecord(entry), input.value);
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

function wireList(panel) {
  renderList(panel);
  store.ready.then(function () { renderList(panel); });
  store.onChange(function () {
    // Leave a rename in progress alone.
    if (!panel.querySelector(".bookmarks-rename")) renderList(panel);
  });
  panel.hidden = false;
}

document.querySelectorAll("[data-course-save]").forEach(wireSave);
document.querySelectorAll("[data-courses-list]").forEach(wireList);
