// static/settingspage.js
//
// The Admin settings page (web/templates/admin_settings.html, ADM.103): search that filters as you type, the
// "Show settings that cannot be changed" switch (remembered in this browser and carried in the address as
// ?show=fixed), a live count of the changes made, a duration box that says what it means in words, a colour picker
// that follows its box, and a rate's "No limit" box. Without script the same page works: the search and the switch
// are a plain GET form (its Apply button), and the server draws the page they ask for.

// Imported with this module's own `?v=` query, as systemmap.js explains.
const VERSION_QUERY = new URL(import.meta.url).search;
const { convert } = await import(`./unitregistry.js${VERSION_QUERY}`);

const STORE_KEY = "planetgen.settings.showFixed";

function stored() {
  try {
    return window.localStorage.getItem(STORE_KEY) === "1";
  } catch {
    return false;
  }
}

function remember(on) {
  try {
    window.localStorage.setItem(STORE_KEY, on ? "1" : "0");
  } catch {
    // private window or blocked storage: the page works without it
  }
}

// 21600 seconds -> "6 hours"; the same words the server writes (web/settings_page.py duration_text).
export function durationText(value, unit) {
  const seconds = convert(Number(value), unit === "hours" ? "h" : "s", "s");
  if (value === "" || !Number.isFinite(seconds) || seconds < 0) return "";
  if (seconds === 0) return "off";
  const parts = [];
  let rest = seconds;
  for (const [name, symbol] of [["day", "d"], ["hour", "h"], ["minute", "min"], ["second", "s"]]) {
    const count = Math.floor(convert(rest, "s", symbol));
    rest -= convert(count, symbol, "s");
    if (count >= 1) parts.push(`${count} ${name}${count === 1 ? "" : "s"}`);
    if (parts.length === 2) break;
  }
  return parts.length ? parts.join(" ") : `${seconds} seconds`;
}

function changed(control) {
  if (control.type === "checkbox" || control.type === "radio") return control.checked !== control.defaultChecked;
  if (control.tagName === "SELECT") return Array.from(control.options).some((o) => o.selected !== o.defaultSelected);
  return control.value !== control.defaultValue;
}

export function init(root) {
  const form = root.querySelector(".settings-form");
  const filter = root.querySelector(".settings-filter");
  const search = root.querySelector("[data-settings-search]");
  const toggle = root.querySelector("[data-show-fixed-toggle]");
  const apply = root.querySelector("[data-settings-apply]");
  const note = root.querySelector("[data-settings-note]");
  const rows = Array.from(root.querySelectorAll("[data-setting]"));
  const groups = Array.from(root.querySelectorAll("[data-group]"));
  const fixedCount = Number(root.dataset.fixedCount || 0);
  if (apply) apply.hidden = true;

  const setFixed = (on, save) => {
    root.dataset.showFixed = on ? "true" : "false";
    if (toggle) toggle.checked = on;
    if (save) {
      remember(on);
      const url = new URL(window.location.href);
      if (on) url.searchParams.set("show", "fixed");
      else url.searchParams.delete("show");
      window.history.replaceState(null, "", url);
    }
    refresh();
  };

  const visible = (row) => !row.classList.contains("setting-fixed") || root.dataset.showFixed === "true";

  function refresh() {
    const terms = (search ? search.value : "").toLowerCase().split(/\s+/).filter(Boolean);
    let shown = 0;
    let hiddenMatches = 0;
    for (const row of rows) {
      const match = terms.every((term) => row.dataset.search.includes(term));
      row.hidden = !match;
      if (match && visible(row)) shown += 1;
      else if (match) hiddenMatches += 1;
      if (match && terms.length) row.closest("details.settings-advanced")?.setAttribute("open", "");
    }
    for (const group of groups) {
      const any = Array.from(group.querySelectorAll("[data-setting]")).some((row) => !row.hidden && visible(row));
      group.hidden = !any;
      const link = root.querySelector(`[data-group-link="${group.dataset.group}"]`);
      if (link) link.hidden = !any;
    }
    if (!note) return;
    if (terms.length && !shown) {
      note.textContent = hiddenMatches
        ? `No setting you can change matches. ${hiddenMatches} that cannot be changed ${hiddenMatches === 1 ? "does" : "do"}.`
        : `No setting matches "${search.value}".`;
    } else if (terms.length && hiddenMatches) {
      note.textContent = `${hiddenMatches} more that cannot be changed match${hiddenMatches === 1 ? "es" : ""}.`;
    } else {
      note.textContent = "";
    }
  }

  if (search) {
    search.addEventListener("input", refresh);
    filter.addEventListener("submit", (event) => event.preventDefault());
  }
  if (toggle) toggle.addEventListener("change", () => setFixed(toggle.checked, true));
  const params = new URL(window.location.href).searchParams;
  if (fixedCount && !params.has("show") && stored()) setFixed(true, false);
  refresh();

  if (!form) return;
  const count = form.querySelector("[data-settings-changes]");
  const save = form.querySelector("[data-settings-save]");
  const idle = count ? count.textContent : "";

  function tally() {
    const touched = new Set();
    form.querySelectorAll("[data-setting]").forEach((row) => {
      const controls = Array.from(row.querySelectorAll("input:not([type=hidden]):not([data-color-picker]), select, textarea"));
      if (controls.some(changed)) touched.add(row);
    });
    if (count) count.textContent = touched.size ? `${touched.size} ${touched.size === 1 ? "change" : "changes"}. Nothing is saved until you confirm the next page.` : idle;
    if (save) save.disabled = touched.size === 0;
  }
  form.addEventListener("input", tally);
  form.addEventListener("change", tally);
  tally();

  form.querySelectorAll("[data-duration]").forEach((box) => {
    const out = box.parentElement.querySelector("[data-duration-text]");
    box.addEventListener("input", () => { if (out) out.textContent = durationText(box.value.trim(), box.dataset.duration); });
  });
  form.querySelectorAll("[data-color-picker]").forEach((picker) => {
    const text = picker.parentElement.querySelector("[data-color-text]");
    picker.addEventListener("input", () => { text.value = picker.value; tally(); });
    text.addEventListener("input", () => { if (/^#[0-9a-fA-F]{6}$/.test(text.value)) picker.value = text.value; });
  });
  form.querySelectorAll(".setting-rate").forEach((group) => {
    const off = group.querySelector("input[type=checkbox]");
    const countBox = group.querySelector("input[type=text]");
    off.addEventListener("change", () => { if (off.checked) countBox.value = ""; });
    countBox.addEventListener("input", () => { off.checked = countBox.value.trim() === ""; });
  });
}

document.querySelectorAll("[data-settings]").forEach(init);
