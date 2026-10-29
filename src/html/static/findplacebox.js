// html/static/findplacebox.js
//
// MAP.336: the Galaxy Map's find-a-place box: one line, always on the map, `/` to reach it. As you type it lists, grouped, what the line
// can be (findplace.js chooses: coordinates, a sector designation, an object ID, charted regions, names the server finds, places
// visited, the way to the whole galaxy), and a result is gone to with Enter or a click: every result is a target for the map's one
// `goTo`. An empty line lists the places visited. The form and the list are in the page (`planetgen.web.maps.galaxymap3d`); this
// fills them.

const VERSION_QUERY = new URL(import.meta.url).search;
const FP = await import(`./findplace.js${VERSION_QUERY}`);

const TYPING_MS = 200;

// deps:
// - parse(text): galaxystages' parseAddress of the line ({sector, point}, {name} or {problem});
// - parseSector(designation): parseSectorDesignation;
// - sectorExists({shell, band, wedge}): whether a sector is a cell of the grid;
// - regions(): the charted regions' blobs;
// - recents: findplace.createRecents of the browser's storage;
// - locate(name): a Promise of the server's name matches;
// - home: {label, target} for the whole galaxy;
// - go(target): the map's goTo;
// - notice(text): the map's one-line notice.
// Returns {focus(), update(), submit(), query()}.
export function createFindBox(form, input, results, deps) {
  let groups = [];
  let names = null;
  let namesFor = "";
  let timer = null;
  let asking = 0;
  let active = -1;
  // The server's answer for the line last asked about, so a line that is typed and then submitted is asked once.
  let lookup = { text: "", promise: null };
  function ask(text) {
    if (lookup.text !== text || !lookup.promise) {
      const promise = deps.locate(text);
      lookup = { text: text, promise: promise };
      // A lookup that failed is asked again next time.
      promise.catch(function () { if (lookup.promise === promise) lookup = { text: "", promise: null }; });
    }
    return lookup.promise;
  }
  const flat = function () { return groups.reduce(function (all, g) { return all.concat(g.items); }, []); };

  function close() {
    results.hidden = true;
    results.textContent = "";
    active = -1;
    input.setAttribute("aria-expanded", "false");
    input.removeAttribute("aria-activedescendant");
  }

  function choose(item) {
    close();
    input.value = "";
    names = null;
    namesFor = "";
    deps.notice("");
    if (!deps.go(item.target) && item.target) deps.notice("Cannot go to " + item.label + " from here.");
    input.blur();
  }

  function render() {
    results.textContent = "";
    const items = flat();
    if (!items.length) {
      close();
      return;
    }
    let n = 0;
    const panel = document.createElement("div");
    panel.className = "galaxy-find-panel";
    results.appendChild(panel);
    groups.forEach(function (group) {
      const section = document.createElement("div");
      section.className = "galaxy-find-group";
      section.setAttribute("role", "group");
      const heading = document.createElement("h4");
      heading.textContent = group.title;
      heading.id = "galaxymap3d-find-" + group.key;
      section.setAttribute("aria-labelledby", heading.id);
      section.appendChild(heading);
      const list = document.createElement("ul");
      list.setAttribute("role", "presentation");
      group.items.forEach(function (item) {
        const row = document.createElement("li");
        row.setAttribute("role", "option");
        row.id = "galaxymap3d-find-option-" + n;
        row.dataset.index = String(n);
        row.dataset.kind = item.kind;
        const name = document.createElement("span");
        name.className = "galaxy-find-name";
        name.textContent = item.label;
        const detail = document.createElement("span");
        detail.className = "galaxy-find-detail";
        detail.textContent = item.detail || "";
        row.appendChild(name);
        row.appendChild(detail);
        row.addEventListener("mousedown", function (event) { event.preventDefault(); });
        row.addEventListener("click", function () { choose(item); });
        list.appendChild(row);
        n += 1;
      });
      section.appendChild(list);
      panel.appendChild(section);
    });
    results.hidden = false;
    input.setAttribute("aria-expanded", "true");
    mark(Math.min(active, items.length - 1));
  }

  function mark(index) {
    active = index;
    results.querySelectorAll('[role="option"]').forEach(function (row) {
      const on = Number(row.dataset.index) === index;
      row.classList.toggle("is-active", on);
      row.setAttribute("aria-selected", on ? "true" : "false");
      if (on) input.setAttribute("aria-activedescendant", row.id);
    });
    if (index < 0) input.removeAttribute("aria-activedescendant");
  }

  function compute() {
    const text = input.value.trim();
    let parsed = text ? deps.parse(text) : {};
    // Coordinates outside the grid are said to be, not flown to.
    if (parsed.point && parsed.sector && !deps.sectorExists(parsed.sector)) {
      parsed = { problem: "The point " + parsed.point.join(", ") + " pc is outside the galaxy." };
    }
    const id = text ? FP.objectIdOf(text, deps.parseSector) : null;
    groups = FP.buildGroups({
      text: text, parsed: parsed.problem ? {} : parsed, id: id, regions: deps.regions(), recents: deps.recents.list(),
      names: namesFor === text ? names : null, home: deps.home,
    });
    return { text: text, parsed: parsed, id: id };
  }

  // Refreshes the list for the line as it is now and, when the line may be a name the server knows, asks it.
  function update() {
    const now = compute();
    render();
    if (!now.id && FP.asksServer(now.parsed, now.text) && namesFor !== now.text) {
      const token = ++asking;
      ask(now.text).then(function (matches) {
        if (token !== asking || input.value.trim() !== now.text) return;
        names = matches;
        namesFor = now.text;
        compute();
        render();
        if (!groups.length) deps.notice("Nothing is named like " + now.text + ".");
        else deps.notice("");
      }, function () {
        if (token === asking) deps.notice("The lookup failed. Please try again shortly.");
      });
    }
    return now;
  }

  input.setAttribute("role", "combobox");
  input.setAttribute("aria-autocomplete", "list");
  input.setAttribute("aria-expanded", "false");
  input.setAttribute("aria-controls", results.id);
  results.setAttribute("role", "listbox");
  results.setAttribute("aria-label", "Places");

  input.addEventListener("input", function () {
    clearTimeout(timer);
    timer = setTimeout(function () { if (document.activeElement === input) update(); }, TYPING_MS);
  });
  input.addEventListener("focus", update);
  input.addEventListener("blur", function () { setTimeout(close, 150); });
  input.addEventListener("keydown", function (event) {
    const items = flat();
    if (event.key === "ArrowDown" || event.key === "ArrowUp") {
      if (!items.length) update();
      const count = flat().length;
      if (!count) return;
      event.preventDefault();
      mark(event.key === "ArrowDown" ? (active + 1) % count : (active - 1 + count) % count);
    } else if (event.key === "Escape" && !results.hidden) {
      event.preventDefault();
      event.stopPropagation();
      close();
    }
  });

  // Enter (the form's submit): the marked result, else the first; a line the server has not answered yet waits for its answer.
  function submit() {
    clearTimeout(timer);
    const now = update();
    if (!now.text && active < 0) {
      deps.notice("Type a designation, shell/band/wedge, x, y, z in pc, an object ID or a name.");
      return;
    }
    if (now.parsed.problem && !now.id && !flat().length) {
      deps.notice(now.parsed.problem);
      return;
    }
    const items = flat();
    if (items.length) {
      choose(items[active >= 0 ? active : 0]);
      return;
    }
    if (FP.asksServer(now.parsed, now.text)) {
      deps.notice("Looking up " + now.text + "…");
      const token = ++asking;
      ask(now.text).then(function (matches) {
        if (token !== asking) return;
        names = matches;
        namesFor = now.text;
        compute();
        const found = flat();
        if (!found.length) {
          deps.notice("Nothing is named like " + now.text + ".");
          return;
        }
        deps.notice("");
        // One match, or one exact name, is gone to; several are listed to pick from.
        const exact = found.filter(function (item) { return item.kind !== "home" && item.label.toLowerCase() === now.text.toLowerCase(); });
        if (found.length === 1 || exact.length === 1) choose(exact[0] || found[0]);
        else {
          deps.notice(found.length + " names match; pick one.");
          render();
        }
      }, function () {
        if (token === asking) deps.notice("The lookup failed. Please try again shortly.");
      });
    }
  }

  // `/` anywhere on the page, outside a field, reaches the box.
  document.addEventListener("keydown", function (event) {
    if (event.key !== "/" || event.ctrlKey || event.metaKey || event.altKey || form.hidden) return;
    const t = event.target;
    const typing = t && (t.isContentEditable || /^(INPUT|TEXTAREA|SELECT)$/.test(t.tagName));
    if (typing) return;
    event.preventDefault();
    input.focus();
    input.select();
  });

  return {
    focus: function () { input.focus(); },
    update: update,
    submit: submit,
    close: close,
    query: function () { return input.value.trim(); },
  };
}
