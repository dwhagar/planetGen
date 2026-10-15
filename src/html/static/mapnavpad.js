// static/mapnavpad.js
//
// The Galaxy Map's on-map controls (MAP.200): a slider for the zoom, a
// direction star of four arrows and a middle button that moves the center
// the view turns and zooms about, and a readout of where the center is while
// it moves. The map itself stays the stage view's (galaxystageview.js); this
// only drives its `zoomShare`, `moveCentre`, `resetCentre` and `centre`.

const STEP = 0.2;
const REPEAT_MS = 110;
const SHOW_MS = 1800;

const ARROWS = [
  { name: "up", label: "Move the center up", glyph: "▲", dx: 0, dy: 1, area: "up" },
  { name: "left", label: "Move the center left", glyph: "◀", dx: -1, dy: 0, area: "left" },
  { name: "reset", label: "Reset the center to the middle of the region being viewed", glyph: "●", area: "mid" },
  { name: "right", label: "Move the center right", glyph: "▶", dx: 1, dy: 0, area: "right" },
  { name: "down", label: "Move the center down", glyph: "▼", dx: 0, dy: -1, area: "down" },
];

function where(point) {
  return "Center " + point.map(function (v) { return Math.round(v * 10) / 10; }).join(", ") + " pc";
}

// Builds the controls inside `host` (the map's viewport) and returns
// {sync()}: call it whenever the camera moves.
export function createNavPad(host, stageView, cross) {
  const doc = host.ownerDocument;
  const pad = doc.createElement("div");
  pad.className = "starmap-navpad";
  pad.id = "galaxymap3d-navpad";

  const slider = doc.createElement("input");
  slider.type = "range";
  slider.className = "starmap-zoom-slider";
  slider.id = "galaxymap3d-zoom-slider";
  slider.min = "0";
  slider.max = "1000";
  slider.step = "1";
  slider.setAttribute("aria-label", "Zoom");
  slider.title = "Zoom";
  pad.appendChild(slider);

  const star = doc.createElement("div");
  star.className = "starmap-dirstar";
  star.setAttribute("role", "group");
  star.setAttribute("aria-label", "Move the center of the view");
  pad.appendChild(star);

  const readout = doc.createElement("div");
  readout.className = "starmap-centre-readout";
  readout.setAttribute("aria-live", "polite");
  readout.hidden = true;
  host.appendChild(readout);

  let hideTimer = null;
  function showCentre() {
    const point = stageView.centre();
    if (!point) return;
    readout.textContent = where(point);
    readout.hidden = false;
    if (cross) cross.classList.add("is-moving");
    clearTimeout(hideTimer);
    hideTimer = setTimeout(function () {
      readout.hidden = true;
      if (cross) cross.classList.remove("is-moving");
    }, SHOW_MS);
  }

  const buttons = {};
  let repeat = null;
  function stop() {
    clearInterval(repeat);
    repeat = null;
  }
  function act(arrow) {
    if (arrow.name === "reset") stageView.resetCentre();
    else stageView.moveCentre(arrow.dx * STEP, arrow.dy * STEP);
    showCentre();
  }
  ARROWS.forEach(function (arrow) {
    const button = doc.createElement("button");
    button.type = "button";
    button.className = "starmap-dir starmap-dir-" + arrow.name;
    button.dataset.dir = arrow.name;
    button.style.gridArea = arrow.area;
    button.setAttribute("aria-label", arrow.label);
    button.title = arrow.label;
    button.textContent = arrow.glyph;
    if (arrow.name === "reset") {
      button.addEventListener("click", function () { act(arrow); });
    } else {
      // A press steps and, held, keeps stepping; a key press (click with no
      // pointer) steps once.
      button.addEventListener("pointerdown", function (event) {
        if (event.button !== 0) return;
        act(arrow);
        stop();
        repeat = setInterval(function () { act(arrow); }, REPEAT_MS);
      });
      ["pointerup", "pointerleave", "pointercancel"].forEach(function (type) { button.addEventListener(type, stop); });
      button.addEventListener("click", function (event) { if (event.detail === 0) act(arrow); });
    }
    star.appendChild(button);
    buttons[arrow.name] = button;
  });

  // The slider follows the zoom (the wheel, the buttons) except while it is
  // being dragged.
  let held = false;
  slider.addEventListener("pointerdown", function () { held = true; });
  ["pointerup", "pointercancel"].forEach(function (type) { slider.addEventListener(type, function () { held = false; }); });
  slider.addEventListener("input", function () { stageView.setZoomShare(Number(slider.value) / 1000); });
  // The wheel over the controls zooms the map as it does anywhere else.
  pad.addEventListener("wheel", function (event) { if (stageView.onWheel(event)) event.preventDefault(); }, { passive: false });
  host.appendChild(pad);

  function sync() {
    const on = stageView.controlsEnabled();
    pad.classList.toggle("is-disabled", !on);
    slider.disabled = !on;
    Object.keys(buttons).forEach(function (name) { buttons[name].disabled = !on; });
    if (!held) slider.value = String(Math.round(stageView.zoomShare() * 1000));
    if (!readout.hidden) {
      const point = stageView.centre();
      if (point) readout.textContent = where(point);
    }
  }
  sync();
  return { sync: sync, element: pad };
}
