// static/adminpalette.js
//
// The admin palette and keyboard set (ADM.120, ADM.122), loaded on every page under the admin shell (base.html).
//
// Ctrl+K (or "/") opens a search over every admin page, task and setting (`/admin/palette.json`, the registry in
// web/adminregistry.py, loaded the first time it opens); the arrows and Enter pick one and open it. A task or setting
// only opens its page, card or row: nothing runs from the palette. Keys, all off while a field has focus:
//   Ctrl+K or /   open the palette          g then g/s/j/d   go to Generate, Settings, Jobs, Data
//   ?             list the keys             Ctrl+Enter       press the page's Review button (data-review)
//   Esc           close what is open
// A line that starts like a typed command ("add sectors radius 12", ADM.121) also lists what it means, from
// `/admin/generate/command`; it opens the task's Review step filled in, and nothing runs until Run is pressed there.
// Without script the shell's search box submits to /admin/search, which lists the same entries.

export const COMMAND_START = /^\s*(add|sectors|scatters?|new galaxy|wipe|check)\b/i;

export const GO_KEYS = { g: "generate", s: "settings", j: "jobs", d: "data" };

// The same ranking as web/adminregistry.py `search`: every word must appear; a word that starts the label ranks above
// one inside it, which ranks above one in the summary or keywords; then pages, tasks, settings; then list order.
export function rank(items, query, limit = 12) {
  const words = query.toLowerCase().split(/\s+/).filter(Boolean);
  const kinds = ["page", "task", "setting"];
  if (!words.length) return items.filter((item) => item.kind === "page").slice(0, limit);
  const scored = [];
  items.forEach((item, position) => {
    const label = item.label.toLowerCase();
    const text = [item.label, item.summary, (item.keywords || []).join(" "), item.kind].join(" ").toLowerCase();
    if (!words.every((word) => text.includes(word))) return;
    const score = words.reduce((sum, word) => sum + (label.startsWith(word) ? 0 : label.includes(word) ? 1 : 2), 0);
    scored.push([score, kinds.indexOf(item.kind), position, item]);
  });
  scored.sort((a, b) => a[0] - b[0] || a[1] - b[1] || a[2] - b[2]);
  return scored.slice(0, limit).map((entry) => entry[3]);
}

// True where a typed key would go into a field, so the single-key shortcuts leave it alone.
export function typing(target) {
  if (!target || !target.tagName) return false;
  const tag = target.tagName.toUpperCase();
  return tag === "INPUT" || tag === "TEXTAREA" || tag === "SELECT" || Boolean(target.isContentEditable);
}

// What a key press means: "open", "help", "close", "review", "go:<destination>" or "" (nothing). `armed` is true
// right after a "g".
export function meaning(event, armed) {
  const key = event.key;
  if ((event.ctrlKey || event.metaKey) && key.toLowerCase() === "k") return "open";
  if ((event.ctrlKey || event.metaKey) && key === "Enter") return "review";
  if (key === "Escape") return "close";
  if (typing(event.target) || event.ctrlKey || event.metaKey || event.altKey) return "";
  if (key === "/") return "open";
  if (key === "?") return "help";
  if (armed && GO_KEYS[key]) return `go:${GO_KEYS[key]}`;
  return "";
}

const KIND_NAMES = { page: "Page", task: "Task", setting: "Setting", command: "Command" };
const DANGER_NAMES = { rewrites: "Changes data", destroys: "Removes data" };

function el(tag, attrs = {}, text = "") {
  const node = document.createElement(tag);
  for (const [name, value] of Object.entries(attrs)) node.setAttribute(name, value);
  if (text) node.textContent = text;
  return node;
}

function init(shell) {
  let items = null;
  let loading = null;
  let armedUntil = 0;
  let active = 0;
  let shown = [];
  let command = null;
  let asked = 0;
  let sizeTimer = null;

  const dialog = el("dialog", { class: "admin-palette", "aria-label": "Admin palette" });
  const input = el("input", { type: "search", placeholder: "Find a page, task or setting", autocomplete: "off",
    "aria-label": "Find a page, task or setting", role: "combobox", "aria-expanded": "true", "aria-controls": "admin-palette-list" });
  const list = el("ul", { id: "admin-palette-list", role: "listbox" });
  const note = el("p", { class: "hint", role: "status" });
  dialog.append(input, list, note);

  const help = el("dialog", { class: "admin-palette admin-keys", "aria-label": "Keyboard shortcuts" });
  const keys = el("dl");
  [["Ctrl+K or /", "Open the palette"], ["g then g, s, j, d", "Go to Generate, Settings, Jobs, Data"],
    ["Ctrl+Enter", "Press Review on this page"], ["?", "This list"], ["Esc", "Close"]].forEach(([key, what]) => {
    keys.append(el("dt", {}, key), el("dd", {}, what));
  });
  help.append(el("h2", {}, "Keyboard"), keys);
  document.body.append(dialog, help);

  function load() {
    if (items) return Promise.resolve(items);
    loading = loading || fetch(shell.dataset.paletteUrl, { credentials: "same-origin", headers: { Accept: "application/json" } })
      .then((response) => (response.ok ? response.json() : Promise.reject(new Error(String(response.status)))))
      .then((body) => { items = body.items; return items; })
      .catch(() => { loading = null; note.textContent = "The list could not be loaded. Use the search box instead."; return []; });
    return loading;
  }

  function ask() {
    const text = input.value;
    command = null;
    clearTimeout(sizeTimer);
    if (!COMMAND_START.test(text) || !shell.dataset.commandUrl) { draw(); return; }
    const mine = ++asked;
    fetch(`${shell.dataset.commandUrl}?q=${encodeURIComponent(text)}`, { credentials: "same-origin", headers: { Accept: "application/json" } })
      .then((response) => (response.ok ? response.json() : null))
      .then((body) => {
        if (mine !== asked || !body) return;
        command = body.command ? { ...body.command, kind: "command", danger: "safe", keywords: [] }
          : body.error ? { error: body.error } : null;
        draw();
        if (command && command.size_url) askSize(command, mine);
      })
      .catch(() => {});
  }

  // ADM.126: a moment after the last key, the size the command would make is added to its row, from the same estimate
  // the Review step shows.
  function askSize(row, mine) {
    sizeTimer = setTimeout(() => {
      fetch(row.size_url, { credentials: "same-origin", headers: { Accept: "application/json" } })
        .then((response) => (response.ok ? response.json() : null))
        .then((body) => {
          if (mine !== asked || !body || body.state !== "ready") return;
          row.size = body.line;
          draw();
        })
        .catch(() => {});
    }, 500);
  }

  function draw() {
    shown = rank(items || [], input.value);
    if (command && command.url) shown = [command, ...shown];
    active = Math.min(active, Math.max(0, shown.length - 1));
    list.replaceChildren(...shown.map((item, index) => {
      const row = el("li", { role: "option", id: `admin-palette-${index}`, "aria-selected": String(index === active) });
      const link = el("a", { href: item.url });
      link.append(el("strong", {}, item.label), el("span", { class: "badge" }, KIND_NAMES[item.kind] || item.kind));
      if (DANGER_NAMES[item.danger]) link.append(el("span", { class: "badge" }, DANGER_NAMES[item.danger]));
      link.append(el("span", { class: "hint" }, item.size ? `${item.summary}. ${item.size}` : item.summary));
      row.append(link);
      return row;
    }));
    input.setAttribute("aria-activedescendant", shown.length ? `admin-palette-${active}` : "");
    if (command && command.error) note.textContent = command.error;
    else if (items) note.textContent = shown.length ? "" : "Nothing matches.";
  }

  function open() {
    if (help.open) help.close();
    if (!dialog.open) dialog.showModal();
    input.value = "";
    active = 0;
    command = null;
    load().then(draw);
    input.focus();
  }

  input.addEventListener("input", () => { active = 0; ask(); });
  input.addEventListener("keydown", (event) => {
    if (event.key === "ArrowDown" || event.key === "ArrowUp") {
      event.preventDefault();
      active = (active + (event.key === "ArrowDown" ? 1 : -1) + shown.length) % Math.max(1, shown.length);
      draw();
    } else if (event.key === "Enter" && shown[active]) {
      event.preventDefault();
      window.location.assign(shown[active].url);
    }
  });
  dialog.addEventListener("click", (event) => { if (event.target === dialog) dialog.close(); });
  help.addEventListener("click", (event) => { if (event.target === help) help.close(); });

  document.addEventListener("keydown", (event) => {
    const mean = meaning(event, Date.now() < armedUntil);
    if (!mean) {
      const plain = !typing(event.target) && !event.ctrlKey && !event.metaKey && !event.altKey;
      armedUntil = event.key === "g" && plain ? Date.now() + 1500 : 0;
      return;
    }
    armedUntil = 0;
    if (mean === "open") { event.preventDefault(); open(); }
    else if (mean === "help") { event.preventDefault(); if (!help.open) help.showModal(); }
    else if (mean === "close") { if (dialog.open) dialog.close(); if (help.open) help.close(); }
    else if (mean === "review") {
      const button = document.querySelector("[data-review]");
      if (button) { event.preventDefault(); button.click(); }
    } else if (mean.startsWith("go:")) {
      const link = shell.querySelector(`[data-destination="${mean.slice(3)}"]`);
      if (link) window.location.assign(link.href);
    }
  });
}

if (typeof document !== "undefined") {
  const shell = document.querySelector("[data-admin-shell]");
  if (shell && shell.dataset.paletteUrl) init(shell);
}
