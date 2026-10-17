// html/static/maplayers.js
//
// The layer interface the Galaxy Map, the Sector Map and the System Map share (MAP.217): one registry a
// toggled layer is declared in once, and the Menu entry, the URL key, the legend line, the picker layer and
// the blocker-fade rule are all made from it, so a new layer (a course, a gravity zone, labels, a region
// selection) is one declaration, not three sets of wiring.
//
// A layer is
//   id               unique name, also what the URL lists ("black-hole", "facility");
//   label            the Menu button's text;
//   modes            the maps it belongs to: any of "galaxy", "sector", "system";
//   group            the Menu heading it sits under ("Show on the map"); layers with no group sit apart;
//   title            optional hover text of the Menu button, or `function () -> string`;
//   available        optional `function () -> boolean`: false leaves the layer's button out of the Menu for now
//                    (a kind the open sector has none of), though the layer keeps its state and its URL entry;
//                    `registry.refresh()` asks again;
//   data             optional `{name: value}` put on the button's `dataset` (existing tests find `data-kind`);
//   defaultOn        shown until the visitor says otherwise (default true);
//   urlKey           the URL parameter that lists the layers of this key turned OFF ("hide"); layers that
//                    begin off need urlKeyOn, the parameter that lists the ones turned ON ("show");
//   legend           null, a string, `{text, swatch, swatchClass}`, or `function (on) -> the same | null`
//                    (`swatch` is a CSS color, `swatchClass` a CSS class that draws the mark);
//   pickerLayer      null, or the picker layer's declaration (mappick.js `addLayer`'s argument), added to a
//                    map's picker while the layer is on;
//   fadedWithSectors true when the blocker fade (MAP.149, nearfield.js) hides the layer with the sectors; a
//                    course layer says false (MAP.121's note);
//   draw             optional `function (scene, on)`, called by `draw(scene)` for every layer of the mode.
//
// The registry holds what is on and off by layer id, including ids not registered yet (a URL names a kind
// before any sector opens it), and writes the lists back to a URL's parameters. Plain DOM calls only, never
// markup from strings: labels can be database content (an object kind's name).
//
// No DOM at import, no three.js: tests/js/maplayers.test.mjs holds it to this description.

export const MODES = ["galaxy", "sector", "system"];

function entryOf(legend, on) {
  const value = typeof legend === "function" ? legend(on) : legend;
  if (!value) return null;
  if (typeof value === "string") return { text: value, swatch: null, swatchClass: null };
  return { text: value.text, swatch: value.swatch || null, swatchClass: value.swatchClass || null };
}

function listOf(value) {
  return String(value || "").split(",").map(function (name) { return name.trim(); }).filter(Boolean);
}

// A registry for `mode` ("galaxy", "sector" or "system"): the layers of the other modes are held but are
// not in its menu, legend, URL or draw.
//
// `options.forgetUnknown` makes a URL's ids that no layer is registered for count for nothing (a map that
// registers all its layers up front); by default they are held for a layer registered later.
export function createLayerRegistry(mode, options) {
  const forgetUnknown = !!(options && options.forgetUnknown);
  if (MODES.indexOf(mode) < 0) throw new Error("unknown map mode: " + mode);
  const layers = [];
  const state = new Map(); // id -> true (on) or false (off), only where it differs from the default
  const listeners = [];
  const origin = new Map(); // id -> the parameter that named it, for ids not registered (yet)
  const views = new Map(); // container -> its refresh, one menu or legend per container

  function find(id) {
    return layers.find(function (layer) { return layer.id === id; }) || null;
  }

  function defaultOf(layer) {
    return layer.defaultOn !== false;
  }

  function mine(layer) {
    return layer.modes.indexOf(mode) >= 0;
  }

  function changed(id) {
    views.forEach(function (refresh) { refresh(); });
    listeners.slice().forEach(function (listener) { listener(id, isOn(id)); });
  }

  // The layer's state: what the visitor or the URL set, else its default; false for an unknown id with
  // nothing set.
  function isOn(id) {
    if (state.has(id)) return state.get(id);
    const layer = find(id);
    return layer ? defaultOf(layer) : true;
  }

  function register(spec) {
    if (!spec || !spec.id || !spec.label) throw new Error("a layer needs an id and a label");
    if (find(spec.id)) throw new Error("layer " + spec.id + " is already registered");
    const modes = spec.modes || [mode];
    modes.forEach(function (name) {
      if (MODES.indexOf(name) < 0) throw new Error("layer " + spec.id + ": unknown map mode " + name);
    });
    if (spec.defaultOn === false && spec.urlKey && !spec.urlKeyOn) {
      throw new Error("layer " + spec.id + " begins off, so it needs urlKeyOn to say when it is on");
    }
    const layer = Object.assign({ group: null, title: null, defaultOn: true, urlKey: null, urlKeyOn: null,
      legend: null, pickerLayer: null, fadedWithSectors: true, draw: null, available: null, data: null }, spec, { modes: modes });
    layers.push(layer);
    if (state.get(layer.id) === defaultOf(layer)) state.delete(layer.id);
    return layer;
  }

  // Takes layers out, for ones made per sector (a kind that is not in the open sector has no button);
  // what was set for the id stays, so it applies again when the layer comes back.
  function unregister(id) {
    const index = layers.findIndex(function (layer) { return layer.id === id; });
    if (index >= 0) layers.splice(index, 1);
  }

  function setOn(id, on) {
    on = !!on;
    const layer = find(id);
    const before = isOn(id);
    if (layer && on === defaultOf(layer)) state.delete(id);
    else state.set(id, on);
    if (before !== isOn(id)) changed(id);
    return isOn(id);
  }

  function toggle(id) {
    return setOn(id, !isOn(id));
  }

  function forMode() {
    return layers.filter(mine);
  }

  // Reads the layers' parameters from a query string or URLSearchParams. What the URL does not name goes
  // back to its default. Ids not registered yet are kept, with the parameter that named them.
  function readUrl(source) {
    const params = typeof source === "string" ? new URLSearchParams(source) : source;
    const keysOff = new Set();
    const keysOn = new Set();
    layers.forEach(function (layer) {
      if (layer.urlKey) keysOff.add(layer.urlKey);
      if (layer.urlKeyOn) keysOn.add(layer.urlKeyOn);
    });
    const before = new Map(state);
    keysOff.forEach(function (key) { listOf(params.get(key)).forEach(function (id) { origin.set(id, { on: false, key: key }); }); });
    keysOn.forEach(function (key) { listOf(params.get(key)).forEach(function (id) { origin.set(id, { on: true, key: key }); }); });
    Array.from(state.keys()).forEach(function (id) {
      const layer = find(id);
      if (!layer || layer.urlKey || layer.urlKeyOn) state.delete(id); // a layer with no URL key keeps what was set
    });
    keysOff.forEach(function (key) { listOf(params.get(key)).forEach(function (id) { state.set(id, false); }); });
    keysOn.forEach(function (key) { listOf(params.get(key)).forEach(function (id) { state.set(id, true); }); });
    layers.forEach(function (layer) { if (state.get(layer.id) === defaultOf(layer)) state.delete(layer.id); });
    if (forgetUnknown) {
      Array.from(state.keys()).forEach(function (id) { if (!find(id)) state.delete(id); });
      origin.clear();
    }
    const ids = new Set(Array.from(before.keys()).concat(Array.from(state.keys())));
    ids.forEach(function (id) { if (before.get(id) !== state.get(id)) changed(id); });
  }

  // Sets or removes the layers' parameters on `params` (a URLSearchParams) and returns it.
  function writeUrl(params) {
    const lists = new Map();
    layers.forEach(function (layer) {
      if (layer.urlKey) lists.set(layer.urlKey, lists.get(layer.urlKey) || []);
      if (layer.urlKeyOn) lists.set(layer.urlKeyOn, lists.get(layer.urlKeyOn) || []);
    });
    state.forEach(function (value, id) {
      const layer = find(id);
      const key = layer ? (value ? layer.urlKeyOn : layer.urlKey) : (origin.has(id) && origin.get(id).on === value ? origin.get(id).key : null);
      if (key && lists.has(key)) lists.get(key).push(id);
    });
    lists.forEach(function (ids, key) {
      if (ids.length) params.set(key, ids.join(","));
      else params.delete(key);
    });
    return params;
  }

  // The legend lines of the layers that are on, in registration order: [{id, text, swatch}].
  function legend() {
    const lines = [];
    forMode().forEach(function (layer) {
      if (!isOn(layer.id)) return;
      const entry = entryOf(layer.legend, true);
      if (entry) lines.push({ id: layer.id, text: entry.text, swatch: entry.swatch, swatchClass: entry.swatchClass });
    });
    return lines;
  }

  // The picker layers to add to a map's picker: those of the layers that are on.
  function pickerLayers() {
    return forMode().filter(function (layer) { return layer.pickerLayer && isOn(layer.id); })
      .map(function (layer) { return layer.pickerLayer; });
  }

  // The share of a layer that shows with the blocker fade at `sectorShare` (0 to 1: how much of the sectors
  // shows). A layer the fade does not hide is whole while it is on.
  function visibleShare(id, sectorShare) {
    const layer = find(id);
    if (!layer || !isOn(id)) return 0;
    if (!layer.fadedWithSectors) return 1;
    return Math.min(1, Math.max(0, sectorShare));
  }

  function draw(scene) {
    forMode().forEach(function (layer) {
      if (layer.draw) layer.draw(scene, isOn(layer.id));
    });
  }

  // --- The Menu and the legend, as DOM ------------------------------------------------------------

  // Fills `container` with a button per layer of the mode that is available (grouped under headings), pressed
  // while the layer is on; a click toggles it. The buttons are made again when the set of available layers
  // changed and brought up to date otherwise. Returns the function that does so.
  function renderMenu(container) {
    const doc = container.ownerDocument || document;
    const buttons = new Map();
    let signature = null;
    function build(list) {
      container.textContent = "";
      buttons.clear();
      let group;
      list.forEach(function (layer) {
        if (layer.group !== group) {
          group = layer.group;
          if (group) {
            const heading = doc.createElement("span");
            heading.className = "galaxy-kinds-heading";
            heading.textContent = group;
            container.appendChild(heading);
          }
        }
        const button = doc.createElement("button");
        button.type = "button";
        button.className = "starmap-btn starmap-toggle";
        button.dataset.layer = layer.id;
        button.textContent = layer.label;
        if (layer.data) Object.keys(layer.data).forEach(function (name) { button.dataset[name] = layer.data[name]; });
        button.addEventListener("click", function () { toggle(layer.id); });
        buttons.set(layer.id, button);
        container.appendChild(button);
      });
    }
    function refresh() {
      const list = forMode().filter(function (layer) { return !layer.available || layer.available(); });
      const now = list.map(function (layer) { return layer.id; }).join("|");
      if (now !== signature) {
        signature = now;
        build(list);
      }
      buttons.forEach(function (button, id) {
        const layer = find(id);
        button.setAttribute("aria-pressed", String(isOn(id)));
        const title = layer && (typeof layer.title === "function" ? layer.title() : layer.title);
        if (title) button.title = title;
      });
    }
    refresh();
    views.set(container, refresh);
    return refresh;
  }

  // Fills `container` with a line per legend entry (a swatch, then the text) and hides it when there is none.
  function renderLegend(container) {
    const doc = container.ownerDocument || document;
    function refresh() {
      const lines = legend();
      container.textContent = "";
      lines.forEach(function (line) {
        const row = doc.createElement("span");
        row.className = "map-layer-legend-line";
        row.dataset.layer = line.id;
        if (line.swatch || line.swatchClass) {
          const swatch = doc.createElement("span");
          swatch.className = "map-layer-legend-swatch" + (line.swatchClass ? " " + line.swatchClass : "");
          swatch.setAttribute("aria-hidden", "true");
          if (line.swatch) swatch.style.background = line.swatch;
          row.appendChild(swatch);
        }
        row.appendChild(doc.createTextNode(line.text));
        container.appendChild(row);
      });
      container.hidden = lines.length === 0;
    }
    views.set(container, refresh);
    refresh();
    return refresh;
  }

  return {
    mode: mode,
    register: register,
    unregister: unregister,
    get: find,
    layers: forMode,
    isOn: isOn,
    setOn: setOn,
    toggle: toggle,
    onChange: function (listener) { listeners.push(listener); },
    // Asks every menu and legend drawn from this registry to bring itself up to date (an `available` or a
    // `title` changed without a layer being switched).
    refresh: function () { views.forEach(function (refresh) { refresh(); }); },
    readUrl: readUrl,
    writeUrl: writeUrl,
    legend: legend,
    pickerLayers: pickerLayers,
    visibleShare: visibleShare,
    draw: draw,
    renderMenu: renderMenu,
    renderLegend: renderLegend,
  };
}
