// static/leversliders.js
//
// The Generate page's galaxy measures as sliders (templates/partials/generate_fields.html `lever_input`; Boss, 2026-11-12).
// Each number or length row of the lever table (planetgen/galaxy/dials.py) is a tick box and a typed box side by side with
// a slider bar under them. Unticked is Auto: the box is blank and disabled (so the plan gets no flag) and the bar is
// greyed at the middle of what Auto draws (data-default). Ticked, the box takes the value and the bar moves it. The bar
// is sticky: dragged near that Auto mark it snaps onto it. The typed box is the form field (name = the lever), so
// the server still reads blank as Auto and the page works without this script, typed into.
//
// generatelevers.js turns a box off for a galaxy too light to have a halo (data-locked on the box and a "lever-lock"
// event); a locked measure unticks, greys and ignores the bar.

const VERSION_QUERY = new URL(import.meta.url).search;
const { fieldsReady } = await import(`./fieldsready.js${VERSION_QUERY}`);
await fieldsReady();

const POSITIONS = 1000;
const STICK = 14; // positions either side of the Auto mark the bar snaps from (1.4 percent of its length)

// Four significant figures for a bar value: the bar is coarse and float noise should never reach a box.
const clean = (value) => (value === 0 ? 0 : Number(Number(value).toPrecision(4)));

class LeverSlider {
  constructor(root) {
    this.root = root;
    this.tick = root.querySelector("[data-lever-set]");
    this.box = root.querySelector("[data-lever]");
    this.range = root.querySelector("[data-lever-range]");
    this.kind = root.dataset.kind;
    this.min = Number(root.dataset.min);
    this.max = Number(root.dataset.max);
    this.log = "log" in root.dataset;
    this.minAbs = "minAbs" in root.dataset ? Number(root.dataset.minAbs) : null;
    this.auto = Number(root.dataset.default);
    this.unit = root.dataset.unit || "";
    this.factors = root.dataset.factors ? JSON.parse(root.dataset.factors) : {};
    this.locked = !!this.box.dataset.locked;
    this.range.hidden = false;
    this.tick.addEventListener("sl-change", () => this.toggle());
    this.tick.addEventListener("change", () => this.toggle());
    this.range.addEventListener("input", () => this.drag());
    this.box.addEventListener("sl-input", () => this.typed());
    this.box.addEventListener("input", () => this.typed());
    root.addEventListener("lever-lock", (event) => this.lock(event.detail.off));
    this.sync(this.box.value.trim() !== "");
  }

  toPosition(value) {
    const clamped = Math.min(Math.max(value, this.min), this.max);
    const share = this.log
      ? (Math.log(clamped) - Math.log(this.min)) / (Math.log(this.max) - Math.log(this.min))
      : (clamped - this.min) / (this.max - this.min);
    return Math.round(share * POSITIONS);
  }

  fromPosition(position) {
    const share = position / POSITIONS;
    let value = this.log
      ? Math.exp(Math.log(this.min) + share * (Math.log(this.max) - Math.log(this.min)))
      : this.min + share * (this.max - this.min);
    if (this.kind === "int") value = Math.round(value);
    value = Math.min(Math.max(clean(value), this.min), this.max);
    if (this.minAbs !== null && Math.abs(value) < this.minAbs) value = value < 0 ? -this.minAbs : this.minAbs;
    return value;
  }

  // The box text for a value: a length with its unit (kiloparsecs from a thousand).
  text(value) {
    if (this.kind !== "length") return String(value);
    return value >= 1000 && this.factors.kpc ? `${clean(value / this.factors.kpc)}kpc` : `${value}${this.unit}`;
  }

  // The number a box holds in the row's unit, or null when it cannot be placed on the bar (empty, a scale length, junk).
  read() {
    const raw = this.box.value.trim();
    if (raw === "") return null;
    if (this.kind !== "length") {
      const value = Number(raw);
      return Number.isFinite(value) ? value : null;
    }
    const match = /^([+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?)\s*([A-Za-z]+)$/.exec(raw);
    if (!match || !(match[2] in this.factors)) return null;
    return Number(match[1]) * this.factors[match[2]];
  }

  // Shows the state the tick says: Auto (blank box, bar greyed at the Auto mark) or a value.
  sync(on) {
    this.tick.checked = on && !this.locked;
    this.paint();
    const value = this.read();
    this.range.value = this.toPosition(on && value !== null ? value : this.auto);
  }

  paint() {
    const on = this.tick.checked && !this.locked;
    this.box.disabled = !on;
    this.range.disabled = !on;
    this.tick.disabled = this.locked;
    this.root.classList.toggle("lever-slider-auto", !on);
  }

  announce() {
    this.box.dispatchEvent(new Event("input", { bubbles: true }));
  }

  toggle() {
    if (this.locked) return;
    if (this.tick.checked) {
      if (this.read() === null) this.box.value = this.text(this.auto);
      this.paint();
      this.range.value = this.toPosition(this.read() ?? this.auto);
    } else {
      this.box.value = "";
      this.paint();
      this.range.value = this.toPosition(this.auto);
    }
    this.announce();
  }

  drag() {
    let position = Number(this.range.value);
    const mark = this.toPosition(this.auto);
    const sticky = Math.abs(position - mark) <= STICK;
    if (sticky) position = mark;
    this.range.value = position;
    this.box.value = this.text(sticky ? this.auto : this.fromPosition(position));
    this.announce();
  }

  typed() {
    const value = this.read();
    if (value !== null && !this.range.matches(":active")) this.range.value = this.toPosition(value);
  }

  lock(off) {
    this.locked = !!off;
    if (off) {
      this.tick.checked = false;
      this.range.value = this.toPosition(this.auto);
    }
    this.paint();
  }
}

for (const root of document.querySelectorAll("[data-lever-slider]")) {
  if (root.querySelector("[data-lever-set]") && root.querySelector("[data-lever]") && root.querySelector("[data-lever-range]")) {
    root.leverSlider = new LeverSlider(root);
  }
}
