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
//   legend           null, a string, `{text, swatch, swatchClass}` or `{text, icon, token}`, or `function (on) -> the
//                    same | null` (`swatch` is a CSS color, `swatchClass` a CSS class that draws the mark; `icon` is a
//                    point-of-interest symbol of poiicons.js and `token` the CSS custom property that colours it, MAP.320);
//   pickerLayer      null, or the picker layer's declaration (mappick.js `addLayer`'s argument), added to a
//                    map's picker while the layer is on;
//   fadedWithSectors true when the blocker fade (MAP.149, nearfield.js) hides the layer with the sectors; a
//                    course layer says false (MAP.121's note);
//   draw             optional `function (scene, on)`, called by `draw(scene)` for every layer of the mode;
//   keyless          a sentence saying why the layer has no line in the key (a setting, names, nothing drawn): a layer
//                    with no `legend`, no key entry pointing at it (`registerKey`, `layer`) and no `keyless` shows up in
//                    `registry.unkeyed()`, which the browser test holds empty (MAP.234);
//   forward          true when the layer can also be pulled forward (MAP.306): a third state beside shown and hidden, which
//                    the layer's own code reads through `isPulled(id)` (it is drawn on top and named first). The Menu
//                    gives such a layer a second button; pulling a layer forward shows it, hiding it lets it go;
//   bulk             (a group) `registry.bulkGroup(group)` gives a Menu group one row of buttons that set every layer of it
//                    to shown, pulled forward or hidden at once.
//
// The registry holds what is on and off by layer id, including ids not registered yet (a URL names a kind
// before any sector opens it), and writes the lists back to a URL's parameters. Plain DOM calls only, never
// markup from strings: labels can be database content (an object kind's name).
//
// No DOM at import, no three.js: tests/js/maplayers.test.mjs holds it to this description.

export const MODES = ["galaxy", "sector", "system"];

const VERSION_QUERY = new URL(import.meta.url).search;
const SVG_NS = "http://www.w3.org/2000/svg";

function entryOf(legend, on) {
  const value = typeof legend === "function" ? legend(on) : legend;
  if (!value) return null;
  if (typeof value === "string") return { text: value, swatch: null, swatchClass: null, icon: null, token: null };
  return { text: value.text, swatch: value.swatch || null, swatchClass: value.swatchClass || null,
    icon: value.icon || null, token: value.token || null };
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
  const pulled = new Set(); // ids pulled forward (MAP.306); `front` in the URL lists them
  const bulkGroups = new Set(); // Menu groups with a row that sets all their layers at once
  const views = new Map(); // container -> its refresh, one menu or legend per container
  const keyEntries = []; // MAP.234: lines of the key for things drawn that are not a toggled layer, or several layers at once

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

  // True while the layer is on and has been pulled forward (a layer that is hidden is not pulled).
  function isPulled(id) {
    const layer = find(id);
    return pulled.has(id) && isOn(id) && !!(layer && layer.forward);
  }

  // Pulls a layer forward or lets it go; pulling one forward shows it.
  function setPulled(id, on) {
    const layer = find(id);
    if (on && !(layer && layer.forward)) return false;
    const before = isPulled(id);
    if (on) {
      pulled.add(id);
      setOn(id, true);
    } else {
      pulled.delete(id);
    }
    if (before !== isPulled(id)) changed(id);
    return isPulled(id);
  }

  // The state of a forward-capable layer as a word: "hidden", "shown" or "forward".
  function stateOf(id) {
    return !isOn(id) ? "hidden" : isPulled(id) ? "forward" : "shown";
  }

  // Sets every forward-capable layer of the mode in `group` to "hidden", "shown" or "forward" at once.
  function setGroup(group, word) {
    forMode().filter(function (layer) { return layer.group === group; }).forEach(function (layer) {
      if (word === "hidden") {
        pulled.delete(layer.id);
        setOn(layer.id, false);
      } else if (word === "forward" && layer.forward) {
        setPulled(layer.id, true);
      } else {
        pulled.delete(layer.id);
        setOn(layer.id, true);
        changed(layer.id);
      }
    });
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
      legend: null, keyless: null, pickerLayer: null, fadedWithSectors: true, draw: null, available: null, data: null }, spec, { modes: modes });
    layers.push(layer);
    if (state.get(layer.id) === defaultOf(layer)) state.delete(layer.id);
    return layer;
  }

  // MAP.234: a line of the key. `text` (a string or a function returning one), `swatch`, `swatchClass`, `icon` and `token` are as a
  // layer's legend entry; `drawn` is `function () -> boolean`, and the line shows only while it is true; `layer` names the layer the line
  // stands for, if any, so the layer needs no legend of its own.
  function registerKey(spec) {
    if (!spec || !spec.id || !spec.text || typeof spec.drawn !== "function") throw new Error("a key line needs an id, a text and a drawn()");
    if (keyEntries.some(function (entry) { return entry.id === spec.id; })) throw new Error("key line " + spec.id + " is already registered");
    const modes = spec.modes || [mode];
    modes.forEach(function (name) {
      if (MODES.indexOf(name) < 0) throw new Error("key line " + spec.id + ": unknown map mode " + name);
    });
    keyEntries.push(Object.assign({ layer: null, swatch: null, swatchClass: null, icon: null, token: null }, spec, { modes: modes }));
  }

  // The key lines of the mode: [{id, layer, drawn}].
  function keys() {
    return keyEntries.filter(function (entry) { return entry.modes.indexOf(mode) >= 0; }).map(function (entry) {
      return { id: entry.id, layer: entry.layer, drawn: !!entry.drawn() };
    });
  }

  // Brings every legend drawn from this registry up to date when what is drawn has changed since the last time (the key lines whose
  // `drawn()` changed with the view). True when it did.
  let keySignature = null;
  function refreshKey() {
    const now = keyEntries.filter(function (entry) { return entry.modes.indexOf(mode) >= 0 && entry.drawn(); })
      .map(function (entry) { return entry.id; }).join("|");
    if (now === keySignature) return false;
    keySignature = now;
    views.forEach(function (refresh) { refresh(); });
    return true;
  }

  // The ids of the layers of the mode that the key says nothing about: no legend, no key line for them and no reason given.
  function unkeyed() {
    return forMode().filter(function (layer) {
      return !layer.legend && !layer.keyless && !keyEntries.some(function (entry) { return entry.layer === layer.id && entry.modes.indexOf(mode) >= 0; });
    }).map(function (layer) { return layer.id; });
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
    if (!on) pulled.delete(id);    // hiding a layer lets it go from the front
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
    const pulledBefore = new Set(pulled);
    pulled.clear();
    listOf(params.get("front")).forEach(function (id) { pulled.add(id); });
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
    pulledBefore.forEach(function (id) { if (!pulled.has(id) && before.get(id) === state.get(id)) changed(id); });
    pulled.forEach(function (id) { if (!pulledBefore.has(id) && before.get(id) === state.get(id)) changed(id); });
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
    const front = Array.from(pulled).filter(function (id) { const layer = find(id); return !layer || (layer.forward && isOn(id)); });
    if (front.length) params.set("front", front.join(","));
    else params.delete("front");
    return params;
  }

  // The legend lines of the layers that are on, in registration order, then the key lines whose thing is drawn: [{id, text, swatch}].
  function legend() {
    const lines = [];
    forMode().forEach(function (layer) {
      if (!isOn(layer.id)) return;
      const entry = entryOf(layer.legend, true);
      if (entry) {
        lines.push({ id: layer.id, text: entry.text, swatch: entry.swatch, swatchClass: entry.swatchClass,
          icon: entry.icon, token: entry.token });
      }
    });
    keyEntries.forEach(function (entry) {
      if (entry.modes.indexOf(mode) < 0 || !entry.drawn()) return;
      lines.push({ id: entry.id, text: typeof entry.text === "function" ? entry.text() : entry.text, swatch: entry.swatch,
        swatchClass: entry.swatchClass, icon: entry.icon, token: entry.token });
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
    const forwards = new Map();
    let signature = null;
    function build(list) {
      container.textContent = "";
      buttons.clear();
      forwards.clear();
      let group;
      list.forEach(function (layer) {
        if (layer.group !== group) {
          group = layer.group;
          if (group) {
            const heading = doc.createElement("span");
            heading.className = "galaxy-kinds-heading";
            heading.textContent = group;
            container.appendChild(heading);
            if (bulkGroups.has(group)) {
              // MAP.306: one row of buttons that sets every layer of the group at once.
              const row = doc.createElement("span");
              row.className = "galaxy-kinds-bulk";
              [["shown", "All shown"], ["forward", "All forward"], ["hidden", "All hidden"]].forEach(function (pair) {
                const all = doc.createElement("button");
                all.type = "button";
                all.className = "starmap-btn";
                all.dataset.bulk = pair[0];
                all.dataset.group = group;
                all.textContent = pair[1];
                const name = group;
                all.addEventListener("click", function () { setGroup(name, pair[0]); });
                row.appendChild(all);
              });
              container.appendChild(row);
            }
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
        if (layer.forward) {
          // MAP.306: the second button of a row, pulled forward or not (it shows the layer too).
          const front = doc.createElement("button");
          front.type = "button";
          front.className = "starmap-btn starmap-toggle map-forward-toggle";
          front.dataset.forward = layer.id;
          front.textContent = "Forward";
          front.setAttribute("aria-label", "Pull " + layer.label.toLowerCase() + " forward");
          front.addEventListener("click", function () { setPulled(layer.id, !isPulled(layer.id)); });
          forwards.set(layer.id, front);
          container.appendChild(front);
        }
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
      forwards.forEach(function (button, id) {
        button.setAttribute("aria-pressed", String(isPulled(id)));
        button.disabled = !isOn(id);
        button.title = "Draw it on top of the other marks, name it first and keep its names when others are cut";
      });
    }
    refresh();
    views.set(container, refresh);
    return refresh;
  }

  // Fills `container` with a line per legend entry (a swatch, then the text) and hides it when there is none. With `options.collapsible`
  // (MAP.234) the lines sit in a "Key" the reader can fold away, which stays folded from one page to the next (`options.storageKey`).
  function renderLegend(container, options) {
    const doc = container.ownerDocument || document;
    const collapsible = !!(options && options.collapsible);
    const storageKey = (options && options.storageKey) || "planetgen.mapKey.folded";
    let fold = null;
    let body = container;
    if (collapsible) {
      fold = doc.createElement("details");
      fold.className = "map-key";
      const summary = doc.createElement("summary");
      summary.textContent = "Key";
      summary.title = "What the map draws; press to fold it away";
      body = doc.createElement("div");
      body.className = "map-key-lines";
      // TEST.519: the lines scroll inside their box, so the box takes the keyboard (arrow keys scroll it) and has a name.
      body.tabIndex = 0;
      body.setAttribute("role", "group");
      body.setAttribute("aria-label", "What the map draws");
      let folded = false;
      try { folded = doc.defaultView.localStorage.getItem(storageKey) === "1"; } catch (error) { /* no storage: it starts open */ }
      fold.open = !folded;
      fold.appendChild(summary);
      fold.appendChild(body);
      fold.addEventListener("toggle", function () {
        try { doc.defaultView.localStorage.setItem(storageKey, fold.open ? "0" : "1"); } catch (error) { /* the choice is not kept */ }
      });
    }
    function refresh() {
      const lines = legend();
      body.textContent = "";
      if (fold && fold.parentNode !== container) container.appendChild(fold);
      lines.forEach(function (line) {
        const row = doc.createElement("span");
        row.className = "map-layer-legend-line";
        row.dataset.layer = line.id;
        if (line.icon) {
          // MAP.320: the kind's own symbol (static/poi/<icon>.svg) in its token's colour; the words beside it name it.
          const symbol = doc.createElementNS(SVG_NS, "svg");
          symbol.setAttribute("class", "poi-icon map-layer-legend-icon poi-" + line.icon);
          symbol.setAttribute("aria-hidden", "true");
          symbol.setAttribute("focusable", "false");
          symbol.dataset.poi = line.icon;
          if (line.token) symbol.style.color = "var(" + line.token + ")";
          const use = doc.createElementNS(SVG_NS, "use");
          use.setAttribute("href", new URL("poi/" + line.icon + ".svg" + VERSION_QUERY + "#poi", import.meta.url).href);
          symbol.appendChild(use);
          row.appendChild(symbol);
        } else if (line.swatch || line.swatchClass) {
          const swatch = doc.createElement("span");
          swatch.className = "map-layer-legend-swatch" + (line.swatchClass ? " " + line.swatchClass : "");
          swatch.setAttribute("aria-hidden", "true");
          if (line.swatch) swatch.style.background = line.swatch;
          row.appendChild(swatch);
        }
        row.appendChild(doc.createTextNode(line.text));
        body.appendChild(row);
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
    registerKey: registerKey,
    keys: keys,
    refreshKey: refreshKey,
    unkeyed: unkeyed,
    unregister: unregister,
    get: find,
    layers: forMode,
    isOn: isOn,
    setOn: setOn,
    toggle: toggle,
    isPulled: isPulled,
    setPulled: setPulled,
    stateOf: stateOf,
    setGroup: setGroup,
    bulkGroup: function (group) { bulkGroups.add(group); },
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
