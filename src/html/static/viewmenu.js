// html/static/viewmenu.js
//
// UX.278: the View menu's shared builder (design: design-drafts/view-menu-many-options-ux274-2026-11-05.md, sections 4 and 5). It makes
// the Menu of a map from the layers that map has (maplayers.js, MAP.217; each layer names its heading in UX.277): a `<details>` per
// heading, with a count in its summary ("Show: 4 of 12 on"), the layers under it as real checkboxes (grouped, with an all/some/none box
// for a group of two or more), a mark where a layer differs from its default, and the open or folded state of each heading kept in the
// visitor's browser. A map shows only the headings it has layers for, so a new layer appears in the right place by declaring a heading.
// Plain DOM calls only (no innerHTML). The System Map and a sector page use it; the Galaxy Map's own controls move onto it in UX.280.
// The search over the entries is viewmenusearch.js (UX.279). tests/js/viewmenu.test.mjs holds the plan and the DOM.

// A menu with more layers than this folds every heading but the first; a small one (the System Map's) stays open, since folding four
// switches away only costs a tap.
export const COLLAPSE_ABOVE = 12;

// The key under which the open or folded state of a heading is kept, per map mode: `planetgen.viewMenu.open.<mode>.<heading>`.
export const STORAGE_PREFIX = "planetgen.viewMenu.open.";

// How a layer is placed: the heading it names, or the first of `headings` ("Show") when it names none the menu knows (the guard test
// fails such a layer). `headings` is the registry's list (maplayers.js HEADINGS, in the order the menu shows them).
function headingOf(layer, headings) {
  return headings.indexOf(layer.heading) >= 0 ? layer.heading : headings[0];
}

// The layers laid out for the menu, as plain data: [{heading, groups: [{name, layers}]}] with the headings in the order of `headings` (only
// those with a layer), and inside a heading the groups in the order their first layer was registered. A layer with no group sits in a
// group with the name null, which has no label and no box.
export function menuPlan(layers, headings) {
  const plan = [];
  headings.forEach(function (heading) {
    const mine = layers.filter(function (layer) { return headingOf(layer, headings) === heading; });
    if (!mine.length) return;
    const groups = [];
    mine.forEach(function (layer) {
      const name = layer.group || null;
      let group = groups.find(function (g) { return g.name === name; });
      if (!group) {
        group = { name: name, layers: [] };
        groups.push(group);
      }
      group.layers.push(layer);
    });
    plan.push({ heading: heading, groups: groups });
  });
  return plan;
}

// "all", "some" or "none" for the layers of a group, given `isOn(id)`.
export function groupState(layers, isOn) {
  const on = layers.filter(function (layer) { return isOn(layer.id); }).length;
  return on === 0 ? "none" : on === layers.length ? "all" : "some";
}

// How many of `layers` are on, and how many differ from their default (a layer pulled forward differs too).
export function tally(layers, isOn, isPulled) {
  let on = 0;
  let changed = 0;
  layers.forEach(function (layer) {
    const shown = !!isOn(layer.id);
    if (shown) on += 1;
    if (differs(layer, shown, isPulled)) changed += 1;
  });
  return { on: on, total: layers.length, changed: changed };
}

// True when the layer is not as it starts: on when it starts off, off when it starts on, or pulled forward.
export function differs(layer, shown, isPulled) {
  return shown !== (layer.defaultOn !== false) || !!(isPulled && isPulled(layer.id));
}

// The words of a count: "4 of 12 on", with ", 2 changed" when something differs from its default.
export function countText(count) {
  return count.on + " of " + count.total + " on" + (count.changed ? ", " + count.changed + " changed" : "");
}

function storageOf(doc, options) {
  if (options && options.storage !== undefined) return options.storage;
  try { return doc.defaultView.localStorage; } catch (error) { return null; }
}

// Fills `container` with the View menu of `registry` and keeps it up to date as layers switch (the registry asks it to through
// `attachView`). Returns the refresh function. `options.storage` is where the open state is kept (default `localStorage`, which can be
// missing: then the state is not kept); `options.collapseAbove` replaces COLLAPSE_ABOVE.
export function renderViewMenu(registry, container, options) {
  const doc = container.ownerDocument || document;
  const storage = storageOf(doc, options);
  const collapseAbove = options && options.collapseAbove != null ? options.collapseAbove : COLLAPSE_ABOVE;
  const opened = new Map(); // heading -> open or not, as set by the visitor or read from storage
  const boxes = new Map(); // layer id -> its checkbox
  const forwards = new Map(); // layer id -> its Forward checkbox
  const rows = new Map(); // layer id -> its row
  const marks = new Map(); // layer id -> its "changed" mark
  const summaries = new Map(); // heading -> the count in its summary
  const groupViews = []; // [{layers, box, count}]
  let signature = null;
  let root = null;

  function keyOf(heading) {
    return STORAGE_PREFIX + registry.mode + "." + heading;
  }

  function remembered(heading) {
    try {
      const text = storage && storage.getItem(keyOf(heading));
      return text === "1" ? true : text === "0" ? false : null;
    } catch (error) {
      return null;
    }
  }

  function remember(heading, open) {
    opened.set(heading, open);
    try { if (storage) storage.setItem(keyOf(heading), open ? "1" : "0"); } catch (error) { /* the choice is not kept */ }
  }

  function startsOpen(heading, index, total) {
    if (opened.has(heading)) return opened.get(heading);
    const kept = remembered(heading);
    if (kept !== null) return kept;
    return total <= collapseAbove || index === 0;
  }

  function dataOf(element, layer) {
    if (layer.data) Object.keys(layer.data).forEach(function (name) { element.dataset[name] = layer.data[name]; });
  }

  function makeRow(layer) {
    const row = doc.createElement("div");
    row.className = "view-menu-row";
    row.dataset.row = layer.id;
    const label = doc.createElement("label");
    label.className = "view-menu-item";
    const box = doc.createElement("input");
    box.type = "checkbox";
    box.className = "view-menu-box";
    box.dataset.layer = layer.id;
    dataOf(box, layer);
    box.addEventListener("change", function () { registry.setOn(layer.id, box.checked); });
    const text = doc.createElement("span");
    text.className = "view-menu-label";
    text.textContent = layer.label;
    const mark = doc.createElement("span");
    mark.className = "view-menu-changed";
    mark.textContent = "changed";
    mark.title = "Not as it starts";
    mark.hidden = true;
    label.appendChild(box);
    label.appendChild(text);
    row.appendChild(label);
    row.appendChild(mark);
    boxes.set(layer.id, box);
    rows.set(layer.id, row);
    marks.set(layer.id, mark);
    if (layer.forward) {
      // MAP.306: the second switch of a row, pulled forward or not (pulling one forward shows it).
      const front = doc.createElement("label");
      front.className = "view-menu-forward";
      const frontBox = doc.createElement("input");
      frontBox.type = "checkbox";
      frontBox.dataset.forward = layer.id;
      frontBox.setAttribute("aria-label", "Pull " + layer.label.toLowerCase() + " forward");
      frontBox.addEventListener("change", function () { registry.setPulled(layer.id, frontBox.checked); });
      front.appendChild(frontBox);
      front.appendChild(doc.createTextNode(" Forward"));
      row.appendChild(front);
      forwards.set(layer.id, frontBox);
    }
    return row;
  }

  function makeGroup(group) {
    const wrap = doc.createElement("div");
    wrap.className = "view-menu-group";
    if (group.name) wrap.dataset.group = group.name;
    if (group.name) {
      wrap.setAttribute("role", "group");
      const head = doc.createElement("div");
      head.className = "view-menu-group-head";
      let box = null;
      if (group.layers.length > 1) {
        // The all/some/none box of a family: pressing it sets every layer of the group.
        box = doc.createElement("input");
        box.type = "checkbox";
        box.className = "view-menu-group-box";
        box.dataset.groupBox = group.name;
        box.setAttribute("aria-label", group.name + ": all or none");
        box.addEventListener("change", function () {
          const turn = groupState(group.layers, registry.isOn) !== "all";
          group.layers.forEach(function (layer) { registry.setOn(layer.id, turn); });
          box.checked = turn;
          box.indeterminate = false;
        });
        head.appendChild(box);
      }
      const name = doc.createElement("span");
      name.className = "view-menu-group-name";
      name.textContent = group.name;
      const count = doc.createElement("span");
      count.className = "view-menu-count";
      head.appendChild(name);
      if (group.layers.length > 1) head.appendChild(count);
      wrap.appendChild(head);
      groupViews.push({ layers: group.layers, box: box, count: group.layers.length > 1 ? count : null });
    }
    group.layers.forEach(function (layer) {
      wrap.appendChild(makeRow(layer));
      if (typeof layer.settings === "function") {
        // UX.277: a layer that needs a slider or a filter draws it under its row.
        const place = doc.createElement("div");
        place.className = "view-menu-settings";
        place.dataset.settings = layer.id;
        wrap.appendChild(place);
        layer.settings(place);
      }
    });
    return wrap;
  }

  function build(list) {
    container.textContent = "";
    boxes.clear();
    forwards.clear();
    rows.clear();
    marks.clear();
    summaries.clear();
    groupViews.length = 0;
    root = doc.createElement("div");
    root.className = "view-menu";
    const plan = menuPlan(list, registry.headings);
    plan.forEach(function (section, index) {
      const details = doc.createElement("details");
      details.className = "view-menu-heading";
      details.dataset.heading = section.heading;
      details.open = startsOpen(section.heading, index, list.length);
      const summary = doc.createElement("summary");
      const name = doc.createElement("span");
      name.className = "view-menu-heading-name";
      name.textContent = section.heading;
      const count = doc.createElement("span");
      count.className = "view-menu-count";
      summary.appendChild(name);
      summary.appendChild(count);
      details.appendChild(summary);
      section.groups.forEach(function (group) { details.appendChild(makeGroup(group)); });
      details.addEventListener("toggle", function () {
        if (opened.get(section.heading) !== details.open) remember(section.heading, details.open);
      });
      summaries.set(section.heading, { count: count, layers: [].concat.apply([], section.groups.map(function (g) { return g.layers; })) });
      root.appendChild(details);
    });
    container.appendChild(root);
  }

  function refresh() {
    const list = registry.layers().filter(function (layer) { return !layer.available || layer.available(); });
    const now = list.map(function (layer) { return layer.id; }).join("|");
    if (now !== signature) {
      signature = now;
      build(list);
    }
    list.forEach(function (layer) {
      const shown = registry.isOn(layer.id);
      const box = boxes.get(layer.id);
      if (!box) return;
      box.checked = shown;
      const title = typeof layer.title === "function" ? layer.title() : layer.title;
      if (title) rows.get(layer.id).title = title;
      const changed = differs(layer, shown, registry.isPulled);
      marks.get(layer.id).hidden = !changed;
      rows.get(layer.id).classList.toggle("view-menu-differs", changed);
      const front = forwards.get(layer.id);
      if (front) {
        front.checked = registry.isPulled(layer.id);
        front.disabled = !shown;
      }
    });
    summaries.forEach(function (view) {
      view.count.textContent = countText(tally(view.layers, registry.isOn, registry.isPulled));
    });
    groupViews.forEach(function (view) {
      const state = groupState(view.layers, registry.isOn);
      if (view.box) {
        view.box.checked = state === "all";
        view.box.indeterminate = state === "some";
      }
      if (view.count) view.count.textContent = countText(tally(view.layers, registry.isOn, registry.isPulled));
    });
  }

  refresh();
  registry.attachView(container, refresh);
  return refresh;
}
