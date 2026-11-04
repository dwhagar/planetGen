// static/preferencespage.js
//
// The Preferences page (web/templates/preferences.html, UX.236): each per-browser setting of mappreferences.js listed with its
// plain-words description and its default, a control that stores a change as it is made, and a Reset to defaults button. The
// settings live in this browser's `localStorage` (`planetgen.pref.<key>`) for now; the same keys move into the account when accounts
// land. Without script the page says so.

// Imported with this module's own `?v=` query, as systemmap.js explains.
const VERSION_QUERY = new URL(import.meta.url).search;
const { preferenceList, readSetting, writeSetting, clearSetting, resetSettings } = await import(`./mappreferences.js${VERSION_QUERY}`);

function storage() {
  try {
    const store = window.localStorage;
    store.getItem("planetgen.pref.probe");
    return store;
  } catch (error) {
    return null;
  }
}

function element(tag, attributes, text) {
  const node = document.createElement(tag);
  Object.keys(attributes || {}).forEach(function (name) { node.setAttribute(name, attributes[name]); });
  if (text !== undefined) node.textContent = text;
  return node;
}

function defaultText(entry) {
  if (entry.kind === "choice") return entry.choices[entry.fallback] ? entry.choices[entry.fallback].split(":")[0] : entry.fallback;
  return entry.fallback + " " + (entry.unit || "");
}

// The controls of one setting: radio buttons for words, a number box for a number. `changed` is told the key after a change is stored.
function controls(entry, store, say, changed) {
  const id = "pref-" + entry.key.replace(/\./g, "-");
  const box = element("div", { class: "preference-controls" });
  if (entry.kind === "choice") {
    const group = element("fieldset", { class: "preference-choice" });
    group.appendChild(element("legend", { class: "sr-only" }, entry.label));
    entry.choiceList = Object.keys(entry.choices);
    entry.choiceList.forEach(function (word) {
      const row = element("label", { class: "preference-option" });
      const radio = element("input", { type: "radio", name: id, value: word, "data-key": entry.key });
      radio.addEventListener("change", function () {
        if (!radio.checked) return;
        say(writeSetting(store, entry.key, word) ? "Saved." : "This browser would not store that.", entry.label);
        changed(entry.key);
      });
      row.appendChild(radio);
      row.appendChild(element("span", {}, " " + entry.choices[word]));
      group.appendChild(row);
    });
    box.appendChild(group);
    return box;
  }
  const label = element("label", { for: id }, entry.label + " (" + entry.min + " to " + entry.max + " " + entry.unit + ")");
  const input = element("input", {
    type: "number", id: id, min: String(entry.min), max: String(entry.max), step: "any", "data-key": entry.key, inputmode: "decimal",
  });
  input.addEventListener("change", function () {
    if (input.value.trim() === "") {
      clearSetting(store, entry.key);     // an emptied box means the default
      say("Back to the default.", entry.label);
    } else if (writeSetting(store, entry.key, input.value)) {
      say("Saved.", entry.label);
    } else {
      say("That is not between " + entry.min + " and " + entry.max + ", so it was not saved.", entry.label);
    }
    changed(entry.key);
  });
  box.appendChild(label);
  box.appendChild(input);
  return box;
}

function show(list, store, root) {
  list.forEach(function (entry) {
    const value = readSetting(store, entry.key);
    root.querySelectorAll('[data-key="' + entry.key + '"]').forEach(function (control) {
      if (control.type === "radio") control.checked = control.value === value;
      else control.value = value === entry.fallback && !hasStored(store, entry.key) ? "" : String(value);
    });
  });
}

function hasStored(store, key) {
  try { return store.getItem("planetgen.pref." + key) !== null; } catch (error) { return false; }
}

function build(root) {
  const store = storage();
  const list = preferenceList();
  const status = root.querySelector("[data-preferences-status]");
  const say = function (text, name) {
    if (status) status.textContent = name + ": " + text;
  };
  const note = root.querySelector("[data-preferences-note]");
  const holder = root.querySelector("[data-preferences-list]");
  if (!store) {
    if (note) {
      note.hidden = false;
      note.textContent = "This browser will not keep preferences (a private window, or storage blocked), so there is nothing to set here.";
    }
    return;
  }
  const resync = function () { show(list, store, root); };
  let currentGroup = "";
  let section = null;
  list.forEach(function (entry) {
    if (entry.group !== currentGroup) {
      currentGroup = entry.group;
      section = element("section", { class: "preference-group" });
      section.appendChild(element("h3", {}, currentGroup));
      holder.appendChild(section);
    }
    const item = element("div", { class: "preference", "data-preference": entry.key });
    item.appendChild(element("h4", {}, entry.label));
    item.appendChild(element("p", { class: "preference-about" }, entry.description));
    item.appendChild(controls(entry, store, say, resync));
    item.appendChild(element("p", { class: "preference-default" }, "Default: " + defaultText(entry)));
    section.appendChild(item);
  });
  const reset = root.querySelector("[data-preferences-reset]");
  if (reset) {
    reset.hidden = false;
    reset.addEventListener("click", function () {
      if (resetSettings(store)) {
        resync();
        if (status) status.textContent = "Every preference is back to its default.";
      } else if (status) {
        status.textContent = "This browser would not clear them.";
      }
    });
  }
  resync();
  if (note) note.hidden = true;
}

const root = document.querySelector("[data-preferences]");
if (root) build(root);
