// static/limitsliders.js
//
// The admin Generate page's limit sliders (web/templates/generate.html `limit_slider`, UX.116): each is a range input
// that moves a number box at its end, and the box shows the value and takes a typed one. The box is the form field, so
// without JavaScript it still works typed into. A slider's stops are in its data-stops (the range input holds an
// index); a typed value between two stops is kept as typed and the slider moves to the nearest, except for a limit
// that takes only presets (data-snap), where it snaps. A limit with data-follows tracks another box times
// data-follow-factor until the user sets it; "Follow the galaxy-wide limits" turns that following back on and a
// segment's Reset button returns its sliders to the defaults.

const VERSION_QUERY = new URL(import.meta.url).search;

const { formatNumber } = await import(`./numberformat.js${VERSION_QUERY}`);

// A value to six significant figures, so a stop's float noise (0.30000000000000004) does not reach a box.
const significant = (value) => {
  const number = Number(value);
  if (!number || !isFinite(number)) return number;
  const scale = 10 ** (5 - Math.floor(Math.log10(Math.abs(number))));
  return Math.round(number * scale) / scale;
};

function nearest(stops, value, scale) {
  let best = 0;
  let distance = Infinity;
  stops.forEach((stop, index) => {
    const gap = scale === "log" ? Math.abs(Math.log(value / stop)) : Math.abs(value - stop);
    if (gap < distance) { distance = gap; best = index; }
  });
  return best;
}

class LimitSlider {
  constructor(root) {
    this.root = root;
    this.range = root.querySelector("[data-limit-range]");
    this.box = root.querySelector("[data-limit-box]");
    this.stops = JSON.parse(root.dataset.stops);
    this.scale = root.dataset.scale;
    this.snap = "snap" in root.dataset;
    this.default = Number(root.dataset.default);
    this.touched = false;
    this.range.addEventListener("input", () => {
      this.put(this.stops[Number(this.range.value)], true);
    });
    this.box.addEventListener("input", () => {
      const value = Number(this.box.value);
      if (this.box.value !== "" && Number.isFinite(value) && value > 0) {
        this.range.value = nearest(this.stops, this.clamp(value), this.scale);
        this.touch();
        this.changed();
      }
    });
    this.box.addEventListener("change", () => {
      const value = Number(this.box.value);
      if (this.box.value === "" || !Number.isFinite(value)) { this.put(this.default, false); return; }
      this.put(this.snap ? this.stops[nearest(this.stops, this.clamp(value), this.scale)] : this.clamp(value), true);
    });
  }

  clamp(value) { return Math.min(Math.max(value, this.stops[0]), this.stops[this.stops.length - 1]); }

  get value() { return Number(this.box.value); }

  // Shows `value`; `byUser` marks a limit the user set, which stops it following another.
  put(value, byUser) {
    this.box.value = String(significant(value));
    this.range.value = nearest(this.stops, this.clamp(value), this.scale);
    if (byUser) this.touch();
    this.changed();
  }

  touch() {
    this.touched = true;
    this.root.dispatchEvent(new CustomEvent("limit-touched", { bubbles: true }));
  }

  changed() { this.root.dispatchEvent(new CustomEvent("limit-change", { bubbles: true })); }
}

const sliders = new Map();
for (const root of document.querySelectorAll("[data-limit-slider]")) {
  if (root.querySelector("[data-limit-range]") && root.querySelector("[data-limit-box]")) {
    sliders.set(root, new LimitSlider(root));
  }
}

// A limit that follows another (the core's, 110% of the galaxy-wide one) until it is set.
function follow(slider) {
  const source = document.getElementById(slider.root.dataset.follows);
  if (!source) return;
  const factor = Number(slider.root.dataset.followFactor);
  const update = () => {
    if (slider.touched) return;
    const value = Number(source.value);
    if (Number.isFinite(value) && value > 0) slider.put(value * factor, false);
  };
  source.closest("[data-limit-slider]").addEventListener("limit-change", update);
  slider.update = update;
}
for (const slider of sliders.values()) {
  if (slider.root.dataset.follows) follow(slider);
}

for (const segment of document.querySelectorAll("[data-limit-segment]")) {
  const members = [...segment.querySelectorAll("[data-limit-slider]")].map((root) => sliders.get(root)).filter(Boolean);
  const tick = segment.querySelector("[data-core-follow]");
  if (tick) {
    // Setting either core limit stops the following (and ticks the box off); ticking it follows again.
    segment.addEventListener("limit-touched", () => {
      for (const member of members) member.touched = true;
      tick.checked = false;
    });
    tick.addEventListener("change", () => {
      if (!tick.checked) return;
      for (const member of members) { member.touched = false; if (member.update) member.update(); }
    });
  }
  const reset = segment.querySelector("[data-limit-reset]");
  if (reset) {
    reset.addEventListener("click", () => {
      for (const member of members) {
        member.touched = false;
        member.put(member.default, false);
        if (member.update) member.update();
      }
      if (tick) tick.checked = true;
      const status = segment.querySelector("[data-limit-status]");
      if (status) status.textContent = "Reset to the defaults.";
    });
  }
}

// The core's floor under the galaxy-wide one gets a warning.
for (const warning of document.querySelectorAll("[data-core-warning]")) {
  const segment = warning.closest("[data-limit-segment]");
  const floor = segment.querySelector("[data-limit-slider][data-follows$='-floor-box'] [data-limit-box]");
  const source = floor && document.getElementById(floor.closest("[data-limit-slider]").dataset.follows);
  if (!floor || !source) continue;
  const check = () => { warning.hidden = !(Number(floor.value) < Number(source.value)); };
  segment.addEventListener("limit-change", check);
  source.closest("[data-limit-slider]").addEventListener("limit-change", check);
  check();
}

// What the cell star ceiling costs: a cell of that many systems takes this long and this much storage.
for (const note of document.querySelectorAll("[data-ceiling-cost]")) {
  const box = document.getElementById(note.dataset.ceilingBox);
  if (!box) continue;
  const perSystem = Number(note.dataset.secondsPerSystem);
  const bytes = Number(note.dataset.bytesPerSystem);
  const show = () => {
    const systems = Number(box.value);
    if (!(systems > 0)) { note.textContent = ""; return; }
    const seconds = systems * perSystem;
    const time = seconds < 90 ? `${Math.round(seconds)} seconds` : seconds < 5400 ? `${Math.round(seconds / 60)} minutes`
      : `${(seconds / 3600).toFixed(1)} hours`;
    const size = systems * bytes;
    const space = size < 1e9 ? `${Math.round(size / 1e6)} MB` : `${(size / 1e9).toFixed(1)} GB`;
    note.textContent = `A full cell makes about ${formatNumber(systems)} systems: roughly ${time} and ${space} each.`;
  };
  box.addEventListener("input", show);
  box.addEventListener("change", show);
  note.closest("[data-limit-segment]").addEventListener("limit-change", show);
  show();
}
