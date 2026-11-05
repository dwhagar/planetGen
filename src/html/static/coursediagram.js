// html/static/coursediagram.js
//
// NAV.69: links the course diagram beside the Galaxy Map (web/maps/navdiagram.py, in the course readout) to the map's course.
//
// - A plain click on a stop of the diagram selects that stop of the course and centres the map on it (`canvas.galaxyCourseSelect`,
//   in galaxymap3d.js) instead of leaving the page; a click with Ctrl, Shift, Alt or the Meta key, a middle click and "open in a new
//   tab" still follow the link to the stop's own page, and so does a keyboard Enter on a stop that is already selected.
// - A click on a stop's ring on the map selects it there (galaxymap3d.js), and the map says so with a `galaxycoursestop` event on its
//   canvas: the diagram marks the same stop (`is-selected`, `aria-current`), in each of its drawings, and so does the list of stops.
//
// The diagram is drawn on the server (web/maps/navdiagram.py), a drawing for each page size, each stop carrying its place in the course as
// `data-index`; nothing here draws anything.

const SELECTED_CLASS = "is-selected";

function stopsOf(doc) {
  return Array.from(doc.querySelectorAll(".nav-diagram-stop[data-index]"));
}

// Marks stop `index` of the course in every drawing of the diagram, and no other.
export function markStop(doc, index) {
  stopsOf(doc).forEach(function (el) {
    const chosen = el.dataset.index === String(index);
    el.classList.toggle(SELECTED_CLASS, chosen);
    if (chosen) el.setAttribute("aria-current", "true");
    else el.removeAttribute("aria-current");
  });
}

// Whether a click should select a stop of the course rather than follow its link: a plain left click.
export function isPlainClick(event) {
  return event.button === 0 && !event.ctrlKey && !event.metaKey && !event.shiftKey && !event.altKey;
}

export function connectCourseDiagram(canvasEl, doc) {
  doc.addEventListener("click", function (event) {
    const target = event.target && event.target.closest ? event.target.closest(".nav-diagram-stop[data-index]") : null;
    if (!target || !isPlainClick(event) || !canvasEl.galaxyCourseSelect) return;
    const index = Number(target.dataset.index);
    // A stop that is selected already opens its page, as a click on its ring on the map does.
    if (canvasEl.galaxyCourseSelected && canvasEl.galaxyCourseSelected() === index) return;
    event.preventDefault();
    canvasEl.galaxyCourseSelect(index);
  });
  canvasEl.addEventListener("galaxycoursestop", function (event) {
    markStop(doc, event.detail.index);
  });
}
