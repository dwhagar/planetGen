// static/fielddefaults.js
//
// The New galaxy page's other properties get the same Default toggle and ? help as the galaxy measures
// (leversliders.js; Boss, 2026-11-13): inside a [data-default-toggles] block each field becomes its title above, a
// ? button beside the title that opens the help text, and an unlabeled Default tick box right beside the box, ticked
// at the start. Ticked is Default: a field that starts blank (the dropdowns, the neighborhood boxes) is blank and
// disabled, so the server draws it; a field that starts with its usual value (prevalence, star mix, the limit
// sliders) shows that value, read-only, and its bar is disabled. Unticked frees it. The tick box has no name, so it
// adds nothing to the form, and without scripts the fields are typed into as before.

const VERSION_QUERY = new URL(import.meta.url).search;
const { fieldsReady } = await import(`./fieldsready.js${VERSION_QUERY}`);
await fieldsReady();

// The confirmation box is typed, and the compact limit has its own "use the star mass" tick box.
const SKIP_NAMES = new Set(["confirm", "compact_min_mass"]);
let serial = 0;

function titleRow(text, helpText) {
  const row = document.createElement("div");
  row.className = "default-title";
  const name = document.createElement("span");
  name.className = "default-title-name";
  name.textContent = text;
  row.append(name);
  let button = null;
  if (helpText) {
    button = document.createElement("button");
    button.type = "button";
    button.className = "lever-help-button";
    button.textContent = "?";
    button.setAttribute("aria-label", `About ${text}`);
    button.setAttribute("aria-expanded", "false");
    helpText.id = helpText.id || `default-help-${++serial}`;
    button.setAttribute("aria-controls", helpText.id);
    helpText.hidden = true;
    button.addEventListener("click", () => {
      helpText.hidden = !helpText.hidden;
      button.setAttribute("aria-expanded", String(!helpText.hidden));
    });
    row.append(button);
  }
  return row;
}

function tickBox(text, checked) {
  const box = document.createElement("sl-checkbox");
  box.checked = checked;
  const hidden = document.createElement("span");
  hidden.className = "sr-only";
  hidden.textContent = `Default for ${text}`;
  box.append(hidden);
  return box;
}

const announce = (element) => element.dispatchEvent(new Event("input", { bubbles: true }));

// A Shoelace input or dropdown with its own label and help text.
function enhanceField(field) {
  const control = field.querySelector(":scope > sl-input, :scope > sl-select");
  if (!control || SKIP_NAMES.has(control.getAttribute("name"))) return;
  const text = control.getAttribute("label") || control.getAttribute("name");
  const help = control.querySelector('[slot="help-text"]');
  if (help) { help.removeAttribute("slot"); }
  const auto = control.dataset.auto !== undefined ? control.dataset.auto : "";
  const prefilled = control.dataset.auto !== undefined;
  const isDefault = prefilled ? String(control.value) === String(auto) : String(control.value) === "";
  const tick = tickBox(text, isDefault);
  const row = document.createElement("div");
  row.className = "default-row";
  control.before(titleRow(text, help));
  control.before(row);
  row.append(tick, control);
  if (help) field.append(help);
  field.classList.add("default-field");
  const apply = (fromUser) => {
    if (tick.checked) {
      if (prefilled) { control.value = auto; control.readOnly = true; } else { control.value = ""; control.disabled = true; }
    } else {
      control.readOnly = false;
      control.disabled = false;
    }
    if (fromUser) announce(control);
  };
  tick.addEventListener("sl-change", () => apply(true));
  apply(false);
}

// A limit slider (limitsliders.js): a native bar and number box under a label.
function enhanceLimit(root) {
  const label = root.querySelector(":scope > label");
  const box = root.querySelector("[data-limit-box]");
  const range = root.querySelector("[data-limit-range]");
  if (!label || !box || !range) return;
  const text = label.textContent.trim();
  const help = root.querySelector(":scope > .field-hint");
  const auto = Number(root.dataset.auto);
  const tick = tickBox(text, Number(box.value) === auto);
  const row = root.querySelector(".limit-slider-row");
  label.replaceWith(titleRow(text, help));
  const head = document.createElement("div");
  head.className = "default-row";
  head.append(tick, box);
  row.before(head);
  root.classList.add("default-field");
  const apply = (fromUser) => {
    const on = tick.checked;
    box.readOnly = on;
    range.disabled = on;
    if (fromUser && on) {
      box.value = String(auto);
      box.dispatchEvent(new Event("change", { bubbles: true }));
    }
  };
  tick.addEventListener("sl-change", () => apply(true));
  apply(false);
}

for (const scope of document.querySelectorAll("[data-default-toggles]")) {
  for (const field of scope.querySelectorAll(".ui-field:not(.ui-field-check)")) {
    if (field.closest(".lever-slider, .default-field")) continue;
    enhanceField(field);
  }
  for (const root of scope.querySelectorAll("[data-limit-slider]:not([data-follows])")) {
    if (!root.closest(".default-field")) enhanceLimit(root);
  }
}

// The help paragraphs of a section (the text under a fold's title or a group's legend) go behind a ? button in the title,
// too, so no help text is left showing. A paragraph that reports something (a data-* attribute, a role, a warning)
// is not help and stays.
function groupHelp(scope) {
  for (const note of scope.querySelectorAll("fieldset > p.hint, fieldset > span.field-hint, details > p.hint")) {
    if (note.hidden || note.classList.contains("warning") || note.hasAttribute("role")
        || [...note.attributes].some((attribute) => attribute.name.startsWith("data-"))) continue;
    const parent = note.parentElement;
    const heading = parent.matches("fieldset") ? parent.querySelector(":scope > legend") : parent.querySelector(":scope > summary");
    if (!heading || heading.querySelector(".lever-help-button")) continue;
    const text = heading.textContent.trim().replace(/\s+/g, " ");
    const row = titleRow(text, note);
    const button = row.querySelector("button");
    // Inside a fold's title the button must not fold it.
    button.addEventListener("click", (event) => { event.preventDefault(); event.stopPropagation(); });
    heading.append(button);
  }
}

for (const scope of document.querySelectorAll("[data-default-toggles]")) groupHelp(scope);

// The sliders that follow another limit (the core's) keep their own "follow" tick box, so they get only the title and
// the ? help; a tick box's own help text moves behind a ? button after it.
for (const scope of document.querySelectorAll("[data-default-toggles]")) {
  for (const root of scope.querySelectorAll("[data-limit-slider][data-follows]")) {
    const label = root.querySelector(":scope > label");
    const help = root.querySelector(":scope > .field-hint");
    if (!label || !help) continue;
    label.replaceWith(titleRow(label.textContent.trim(), help));
  }
  for (const tick of scope.querySelectorAll("sl-checkbox")) {
    const help = tick.querySelector('[slot="help-text"]');
    if (!help || tick.closest(".default-row")) continue;
    help.removeAttribute("slot");
    const row = titleRow(tick.textContent.trim().slice(0, 60), help);
    const button = row.querySelector("button");
    button.setAttribute("aria-label", `About: ${tick.textContent.trim().slice(0, 60)}`);
    tick.after(button, help);
  }
}
