// static/mapnavpad.js
//
// The Galaxy Map's on-map controls (MAP.200): a direction star of four
// arrows and a middle button that moves the center the view turns and zooms
// about, and a readout of where the center is while it moves; and (MAP.336)
// the log ruler: one vertical scale for the camera's distance, the whole
// galaxy at the top and a moon at the bottom, with named stops that are each
// a Go to for that scale and a readout of the distance. The map itself stays
// the stage view's (galaxystageview.js); this only drives its `zoomShare`,
// `setZoomShare`, `scaleInfo`, `goTo`, `moveCentre`, `resetCentre` and `centre`.

const VERSION_QUERY = new URL(import.meta.url).search;
const R = await import(`./maplogruler.js${VERSION_QUERY}`);
const { formatDistancePc } = await import(`./unitladder.js${VERSION_QUERY}`);

const STEP = 0.2;
// Two stop names closer than this on the ruler (CSS pixels) would run into each other: the nearer one shows its tick only.
const LABEL_GAP_PX = 14;
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

  // The ruler (MAP.336): the zoom slider with the named stops beside it, and the readout of the distance.
  const ruler = doc.createElement("div");
  ruler.className = "starmap-ruler";
  ruler.id = "galaxymap3d-ruler";
  ruler.setAttribute("role", "group");
  ruler.setAttribute("aria-label", "Zoom scale: from the whole galaxy at the top to a moon at the bottom");
  const rulerBody = doc.createElement("div");
  rulerBody.className = "starmap-ruler-body";
  ruler.appendChild(rulerBody);
  const stopsEl = doc.createElement("div");
  stopsEl.className = "starmap-ruler-stops";
  rulerBody.appendChild(stopsEl);

  const slider = doc.createElement("input");
  slider.type = "range";
  slider.className = "starmap-zoom-slider";
  slider.id = "galaxymap3d-zoom-slider";
  slider.min = "0";
  slider.max = "1000";
  slider.step = "1";
  slider.setAttribute("aria-label", "Zoom");
  slider.title = "Zoom";
  rulerBody.appendChild(slider);

  // The distance readout goes in the map's lower-left stack, over the scale bar and under the legend, so it covers neither.
  const altitude = doc.createElement("div");
  altitude.className = "starmap-altitude";
  altitude.id = "galaxymap3d-altitude";
  const stack = host.querySelector(".map-corner-stack");
  const scale = stack ? stack.querySelector(".starmap-scale") : null;
  if (scale) stack.insertBefore(altitude, scale);
  else if (stack) stack.appendChild(altitude);
  else host.appendChild(altitude);

  let stops = [];
  let stopsFor = "";
  const stopButtons = {};
  function buildStops(info) {
    const key = info.farPc + "/" + info.floorPc;
    if (key === stopsFor) return;
    stopsFor = key;
    stops = R.stopsBetween(info.farPc, info.floorPc);
    stopsEl.textContent = "";
    stops.forEach(function (stop) {
      const button = doc.createElement("button");
      button.type = "button";
      button.className = "starmap-ruler-stop";
      button.dataset.stop = stop.key;
      const words = stop.label + ", about " + formatDistancePc(stop.pc);
      button.title = "Go to " + words + " from the middle of the map";
      button.setAttribute("aria-label", "Go to " + words);
      button.textContent = stop.label;
      button.addEventListener("click", function () {
        const centre = stageView.centre();
        if (centre) stageView.goTo({ center: centre, dist: stop.pc });
      });
      stopsEl.appendChild(button);
      stopButtons[stop.key] = button;
    });
    layoutRuler();
  }
  // Each stop at its share of the ruler's height, the thumb's own half-height in from each end; a name too close under the one above
  // keeps its tick and loses its words (the title and the label for a screen reader stay).
  const THUMB_PX = 14;
  function placeStops() {
    const height = rulerBody.clientHeight;
    if (!height) return;
    let lastShown = -Infinity;
    stops.forEach(function (stop) {
      const y = THUMB_PX / 2 + (height - THUMB_PX) * stop.share;
      const button = stopButtons[stop.key];
      button.style.top = y + "px";
      const bare = y - lastShown < LABEL_GAP_PX;
      button.classList.toggle("is-bare", bare);
      if (!bare) lastShown = y;
    });
  }
  // A map too narrow for the stops' names keeps their ticks (and so does one with a card standing beside the ruler, style.css).
  const NARROW_REM = 40;
  function layoutRuler() {
    const rem = parseFloat(getComputedStyle(doc.documentElement).fontSize) || 16;
    ruler.classList.toggle("is-slim", host.clientWidth < NARROW_REM * rem);
    placeStops();
  }
  if (typeof ResizeObserver === "function") {
    new ResizeObserver(placeStops).observe(rulerBody);
    new ResizeObserver(layoutRuler).observe(host);
  }

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
  [pad, ruler].forEach(function (control) {
    control.addEventListener("wheel", function (event) { if (stageView.onWheel(event)) event.preventDefault(); }, { passive: false });
  });
  host.appendChild(ruler);
  host.appendChild(pad);

  function syncReadout() {
    const info = stageView.scaleInfo();
    if (!info) {
      altitude.hidden = true;
      return;
    }
    buildStops(info);
    altitude.hidden = false;
    // Only a change is written: the camera's every frame would otherwise repaint the words when they have not moved.
    const text = R.readoutText(info.dist, info.nearest, formatDistancePc);
    if (altitude.textContent !== text) altitude.textContent = text;
    const scale = R.scaleNamed(info.dist, stops);
    const title = scale ? "About the scale of: " + scale.label : "";
    if (altitude.title !== title) altitude.title = title;
    Object.keys(stopButtons).forEach(function (key) {
      stopButtons[key].classList.toggle("is-near", !!scale && scale.key === key);
    });
  }

  function sync() {
    const on = stageView.controlsEnabled();
    pad.classList.toggle("is-disabled", !on);
    ruler.classList.toggle("is-disabled", !on);
    slider.disabled = !on;
    Object.keys(buttons).forEach(function (name) { buttons[name].disabled = !on; });
    if (!held) slider.value = String(Math.round(stageView.zoomShare() * 1000));
    syncReadout();
    if (!readout.hidden) {
      const point = stageView.centre();
      if (point) readout.textContent = where(point);
    }
  }
  sync();
  return { sync: sync, element: pad };
}
