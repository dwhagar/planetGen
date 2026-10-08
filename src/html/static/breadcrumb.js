// html/static/breadcrumb.js
//
// The one breadcrumb every page and map uses (NAV.14, from MAP.93 and
// MAP.94): a single line at any width. When the steps don't all fit it
// shows the first, a "…" button (a menu of the steps it hides) and as many
// of the last ones as fit, then the current one, measured again whenever
// the line's width changes; the current step is cut short with an ellipsis
// only if even that is too long. A page can also give a Steps box (a
// <details> with a [data-steps-panel]) that lists every step, which the
// Galaxy Map shows at phone width in the line's place.
//
// The steps are `{label, last, href?}` (the picker's trail, picker.js):
// a step with an `href` is a link, any other a button that calls
// `onSelect(step)`; the last is the current one, as text.
//
// Plain DOM calls only: labels are database content.

// Draws the breadcrumb into `nav` (a <nav> holding or getting an <ol>).
// options: onSelect(step); steps (the Steps <details>, optional);
// trailing (an element, or a function giving one, kept after the list,
// e.g. the ☆ button).
export function createBreadcrumb(nav, options) {
  options = options || {};
  let items = [];
  let width = -1;
  // The line's "…" menu, while it has one.
  let moreMenu = null;

  // One step: the current one as text, a link or a button going there.
  function node(step, className) {
    if (step.last) {
      const here = document.createElement("span");
      here.textContent = step.label;
      here.setAttribute("aria-current", step.current || "location");
      return here;
    }
    if (step.href) {
      const link = document.createElement("a");
      link.href = step.href;
      link.className = className;
      link.textContent = step.label;
      return link;
    }
    const button = document.createElement("button");
    button.type = "button";
    button.className = className;
    button.textContent = step.label;
    button.addEventListener("click", function () {
      if (options.onSelect) options.onSelect(step);
    });
    return button;
  }

  // A <details> menu of `steps` (a summary showing `text`, labelled
  // `label`); a step taken closes it.
  function stepsMenu(steps, text, label, className) {
    const menu = document.createElement("details");
    menu.className = className;
    const summary = document.createElement("summary");
    summary.textContent = text;
    summary.setAttribute("aria-label", label);
    summary.title = label;
    menu.appendChild(summary);
    const list = document.createElement("ul");
    steps.forEach(function (step) {
      const item = document.createElement("li");
      const part = node(step, "");
      if (!step.last) part.addEventListener("click", function () { menu.open = false; });
      item.appendChild(part);
      list.appendChild(item);
    });
    menu.appendChild(list);
    return menu;
  }

  // The line's steps, those between the first and the last `keep` folded
  // into "…" (none folded while keep covers them all).
  function fill(list, keep) {
    list.textContent = "";
    const folded = keep >= items.length - 1 ? [] : items.slice(1, items.length - keep);
    const shown = folded.length ? [items[0], null].concat(items.slice(items.length - keep)) : items;
    moreMenu = null;
    shown.forEach(function (step) {
      const item = document.createElement("li");
      if (step) {
        item.appendChild(node(step, "crumb"));
      } else {
        item.className = "crumb-more";
        moreMenu = stepsMenu(folded, "…", folded.length === 1 ? "1 more step" : folded.length + " more steps", "crumb-menu");
        item.appendChild(moreMenu);
      }
      list.appendChild(item);
    });
  }

  function listOf() {
    let list = nav.querySelector("ol");
    if (!list) {
      list = document.createElement("ol");
      nav.insertBefore(list, nav.firstChild);
    }
    return list;
  }

  function layout() {
    const list = nav.querySelector("ol");
    if (!list) return;
    width = nav.clientWidth || 0;
    fill(list, items.length);
    // Hidden, or no layout (a phone, where the line gives way): nothing
    // to measure.
    if (!list.clientWidth) return;
    // Measured at full length: the current step is cut short only when
    // even the shortest line doesn't fit.
    list.classList.add("crumbs-measuring");
    for (let keep = items.length - 2; keep >= 1 && list.scrollWidth > list.clientWidth + 1; keep--) {
      fill(list, keep);
    }
    list.classList.remove("crumbs-measuring");
  }

  function renderSteps() {
    const box = options.steps;
    const panel = box && box.querySelector("[data-steps-panel]");
    if (!panel) return;
    panel.textContent = "";
    const list = document.createElement("ol");
    items.forEach(function (step) {
      const item = document.createElement("li");
      const part = node(step, "");
      if (!step.last) part.addEventListener("click", function () { box.open = false; });
      item.appendChild(part);
      list.appendChild(item);
    });
    panel.appendChild(list);
  }

  // Shows `steps`.
  function set(steps) {
    items = steps;
    listOf();
    const extra = typeof options.trailing === "function" ? options.trailing() : options.trailing;
    if (extra) nav.appendChild(extra);
    nav.classList.add("crumbs", "crumbs-ready");
    layout();
    renderSteps();
  }

  // The line is measured again when its width changes (a resize, the
  // panel shown).
  if (typeof ResizeObserver === "function") {
    new ResizeObserver(function () {
      if (nav.clientWidth !== width) layout();
    }).observe(nav);
  }

  // A "…" or Steps menu closes on a click elsewhere or Escape.
  function openMenus() {
    return [moreMenu, options.steps].filter(function (menu) { return menu && menu.open; });
  }
  document.addEventListener("click", function (event) {
    const path = event.composedPath ? event.composedPath() : [];
    openMenus().forEach(function (menu) {
      if (path.indexOf(menu) < 0) menu.open = false;
    });
  });
  document.addEventListener("keydown", function (event) {
    if (event.key !== "Escape") return;
    openMenus().forEach(function (menu) {
      menu.open = false;
      const summary = menu.querySelector("summary");
      if (summary && menu.contains(event.target)) summary.focus();
    });
  });

  return { set: set, layout: layout };
}

// The steps a server-rendered breadcrumb holds (its <li>s: a link, or the
// current page's text).
export function stepsOfList(nav) {
  const list = nav.querySelector("ol");
  const rows = list ? Array.from(list.children).filter(function (row) { return row.tagName === "LI"; }) : [];
  return rows.map(function (row, at) {
    const link = row.querySelector("a");
    return {
      label: (link || row).textContent.trim(),
      href: link ? link.getAttribute("href") : null,
      last: at === rows.length - 1,
      current: "page",
    };
  });
}
