// html/static/systemnav.js
//
// The NAV buttons of a body's info panel on the System Map, in both views
// (NAV.50). The server gives the System Map's root element URL templates
// (`planetgen/web/system_pages.py` `body_nav`), each with `{ref}` where the
// body's object reference goes (`star:3`, `planet:12`, `moon:40`,
// `belt:5`, `comet:8`):
// - data-nav-start, data-nav-end: "Start Here" and "End Here", which open
//   the NAV page with that end chosen;
// - data-nav-take and data-nav-label: while a start or destination is
//   being picked (`?pick=...`), the one button that returns to the NAV page
//   with the course.

const REF = /^(star|planet|moon|belt|comet):\d+$/;

function fill(template, ref) {
  return template.replace("{ref}", encodeURIComponent(ref).replace(/%3A/gi, ":"));
}

// The buttons for the body `ref` as [{label, href, primary}], read from
// `root`'s data attributes; none for anything that is not a body.
export function navActions(root, ref) {
  if (!root || !REF.test(String(ref || ""))) return [];
  const data = root.dataset;
  if (data.navTake) return [{ label: data.navLabel || "Use This", href: fill(data.navTake, ref), primary: true }];
  const actions = [];
  if (data.navStart) actions.push({ label: "Start Here", href: fill(data.navStart, ref) });
  if (data.navEnd) actions.push({ label: "End Here", href: fill(data.navEnd, ref) });
  return actions;
}

// Adds the buttons to `panel` as one row of links.
export function appendNavActions(panel, root, ref) {
  const actions = navActions(root, ref);
  if (!actions.length) return;
  const row = document.createElement("p");
  row.className = "page-actions map-info-actions";
  actions.forEach(function (action) {
    const link = document.createElement("a");
    link.className = action.primary ? "btn starmap-pick" : "btn btn-small btn-secondary";
    link.href = action.href;
    link.textContent = action.label;
    row.appendChild(link);
  });
  panel.appendChild(row);
}
