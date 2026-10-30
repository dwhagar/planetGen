// html/static/mappick.js
//
// The picking, hover and info-panel layer every 3D map shares (MAP.65):
// the Galaxy Map (galaxymap3d.js and its drill-down, galaxystageview.js)
// and the Sector Map's scene (sectorscene.js) hand it their objects and get back
// what is under the pointer, a tooltip that follows it, a ring around
// what is hovered or selected, and an info panel built the same way on
// every map (fields, the NAV links or the pick button, a page link,
// Generate buttons and the ☆ Bookmark button).
//
// Picking: a picker holds layers in priority order. Each layer offers its
// best candidate under a screen point, with its distance from the camera:
// - points: entries (anything with x, y, z) drawn a fixed number of
//   pixels across, picked on screen within `reach(entry)` pixels
//   (mapcore.nearestOnScreen);
// - meshes: three.js objects, raycast; `entryOf(hit)` says what a hit
//   means (null: nothing), and the optional `visible(hit, ray)` leaves out
//   a hit that is drawn too faint or is a back face (MAP.149);
// - custom: `pick(ctx)` returns {entry, distance} or null itself, given
//   {clientX, clientY, rect, ray, camera}.
// A layer marked `occludes` is solid: what sits behind its nearest hit
// can't be picked through it. Otherwise the first layer (by `priority`,
// then the order they were added) with a candidate wins, so a point of
// light beats the see-through block or cloud it sits in. `enabled()`
// leaves a layer out while it returns false. A map may give a layer more
// for itself: `tooltip(entry)` (its hover text), `hover(entry or null)`
// and `select(entry)`.
//
// Plain DOM calls only, never markup from strings: names are database
// content.

const VERSION_QUERY = new URL(import.meta.url).search;
const THREE = await import(`./vendor/three.module.min.js${VERSION_QUERY}`);
const { addField, makeRingTexture, nearestOnScreen } = await import(`./mapcore.js${VERSION_QUERY}`);
const B = await import(`./bookmarks.js${VERSION_QUERY}`);

// --- Picking -----------------------------------------------------------------------

// A picker for `camera` drawing on `canvasEl`.
export function createPicker(camera, canvasEl) {
  const layers = [];
  const raycaster = new THREE.Raycaster();
  const at = new THREE.Vector3();

  // Adds a layer below every one added before it with the same or a
  // lower `priority` (default 0; a higher number is picked later);
  // returns it.
  function addLayer(layer) {
    const rank = layer.priority || 0;
    let index = layers.length;
    while (index > 0 && (layers[index - 1].priority || 0) > rank) index -= 1;
    layers.splice(index, 0, layer);
    return layer;
  }

  function removeLayer(layer) {
    const index = layers.indexOf(layer);
    if (index >= 0) layers.splice(index, 1);
  }

  function candidate(layer, ctx) {
    if (layer.enabled && !layer.enabled()) return null;
    if (layer.points) {
      const found = nearestOnScreen(layer.points(), camera, ctx.rect, ctx.clientX, ctx.clientY, {
        reach: layer.reach, accept: layer.accept, lastWins: layer.lastWins,
      });
      if (!found) return null;
      at.set(found.entry.x, found.entry.y, found.entry.z);
      // An entry may stand in for another (a binary's companion for its primary).
      return { entry: found.entry.pickAs || found.entry, distance: camera.position.distanceTo(at), px: found.px };
    }
    if (layer.meshes) {
      const hits = raycaster.intersectObjects(layer.meshes(), false);
      for (let n = 0; n < hits.length; n++) {
        // `visible(hit, ray)` leaves out a hit that can't be picked (drawn
        // too faint, or a back face): what is behind it can.
        if (layer.visible && !layer.visible(hits[n], ctx.ray)) continue;
        const entry = layer.entryOf ? layer.entryOf(hits[n]) : hits[n].object;
        if (entry != null) return { entry: entry, distance: hits[n].distance, hit: hits[n] };
      }
      return null;
    }
    return layer.pick ? layer.pick(ctx) : null;
  }

  // What is under a screen point, as {layer, entry, distance}, or null.
  // `only(layer)` leaves out the layers it returns false for.
  function pick(clientX, clientY, only) {
    const rect = canvasEl.getBoundingClientRect();
    if (!rect.width || !rect.height) return null;
    const ndc = new THREE.Vector2(((clientX - rect.left) / rect.width) * 2 - 1, -((clientY - rect.top) / rect.height) * 2 + 1);
    raycaster.setFromCamera(ndc, camera);
    const ctx = { clientX: clientX, clientY: clientY, rect: rect, ray: raycaster.ray, camera: camera };
    const found = [];
    let solid = null;
    layers.forEach(function (layer) {
      if (only && !only(layer)) return;
      const c = candidate(layer, ctx);
      if (!c) return;
      c.layer = layer;
      found.push(c);
      if (layer.occludes && (!solid || c.distance < solid.distance)) solid = c;
    });
    for (let n = 0; n < found.length; n++) {
      const c = found[n];
      if (c.layer.occludes) return solid;
      if (!solid || c.distance <= solid.distance) return c;
    }
    return solid;
  }

  return { addLayer: addLayer, removeLayer: removeLayer, pick: pick };
}

// --- Tooltip ------------------------------------------------------------------------

// The hover tooltip `tipEl` (absolutely placed inside its parent, the
// map's viewport): show(text, clientX, clientY) puts it beside the
// pointer, kept inside the viewport; show("") or hide() hides it.
export function createTooltip(tipEl) {
  function hide() {
    if (tipEl) tipEl.hidden = true;
  }

  function show(text, clientX, clientY) {
    if (!tipEl) return;
    if (!text) {
      hide();
      return;
    }
    const rect = tipEl.parentElement.getBoundingClientRect();
    tipEl.textContent = text;
    tipEl.hidden = false;
    const x = Math.min(clientX - rect.left + 12, rect.width - tipEl.offsetWidth - 4);
    const y = Math.min(clientY - rect.top + 12, rect.height - tipEl.offsetHeight - 4);
    tipEl.style.left = Math.max(4, x) + "px";
    tipEl.style.top = Math.max(4, y) + "px";
  }

  // Beside a point in the scene (a choice picked by keyboard).
  function showAt(text, camera, canvasEl, x, y, z) {
    const point = new THREE.Vector3(x, y, z).project(camera);
    const rect = canvasEl.getBoundingClientRect();
    show(text, rect.left + ((point.x + 1) / 2) * rect.width, rect.top + ((1 - point.y) / 2) * rect.height);
  }

  return { show: show, showAt: showAt, hide: hide };
}

// --- Rings ----------------------------------------------------------------------------

// A ring in `color` around what is selected or hovered, added to `scene`:
// at(x, y, z, {px}) keeps it `px` pixels across on screen whatever its
// distance (a point of light); at(x, y, z, {radius}) makes it `radius`
// world units across its middle (a body or a cloud); hide() hides it.
// update() keeps a pixel-sized ring's size as the canvas changes: call it
// every frame. `options`: opacity, depthTest (default true), renderOrder.
export function createRing(scene, camera, canvasEl, color, options) {
  const o = options || {};
  const texture = makeRingTexture(color);
  function sprite(attenuate) {
    const s = new THREE.Sprite(new THREE.SpriteMaterial({
      map: texture, transparent: true, depthWrite: false, depthTest: o.depthTest !== false,
      opacity: o.opacity != null ? o.opacity : 1, sizeAttenuation: attenuate,
    }));
    s.visible = false;
    s.renderOrder = o.renderOrder || 6;
    scene.add(s);
    return s;
  }
  const screenRing = sprite(false);
  const worldRing = sprite(true);
  let px = 0;
  let target = null;

  function hide() {
    screenRing.visible = false;
    worldRing.visible = false;
    target = null;
  }

  // A sprite without size attenuation is scaled in units of the view's
  // height at distance 1, so a ring `px` pixels across needs
  // px / heightPx * 2 * tan(fov / 2).
  function update() {
    if (!screenRing.visible) return;
    const heightPx = canvasEl.clientHeight || 1;
    const scale = (px / heightPx) * 2 * Math.tan(THREE.MathUtils.degToRad(camera.fov) / 2);
    screenRing.scale.set(scale, scale, 1);
  }

  function at(x, y, z, size) {
    target = [x, y, z];
    if (size && size.radius != null) {
      screenRing.visible = false;
      worldRing.position.set(x, y, z);
      worldRing.scale.set(size.radius * 2, size.radius * 2, 1);
      worldRing.visible = true;
      return;
    }
    worldRing.visible = false;
    px = (size && size.px) || 20;
    screenRing.position.set(x, y, z);
    screenRing.visible = true;
    update();
  }

  // Changes a pixel-sized ring's size (a point that grows with the zoom).
  function setPx(next) {
    px = next;
    update();
  }

  return {
    at: at, hide: hide, update: update, setPx: setPx,
    visible: function () { return !!target; },
    position: function () { return target; },
    sprites: [screenRing, worldRing],
  };
}

// --- Info panel --------------------------------------------------------------------

function el(tag, className, text) {
  const node = document.createElement(tag);
  if (className) node.className = className;
  if (text != null) node.textContent = text;
  return node;
}

// A plain <a href> (GET, bookmarkable, opens in a new tab).
export function pageLink(href, label, className) {
  const link = el("a", className || "btn", label);
  link.href = href || "#";
  return link;
}

// The info panel `panelEl` (an <aside class="starmap-info">): show(spec)
// replaces what it shows with one object's, in this order:
// - title: the heading;
// - fields: [label, value] rows (an empty value leaves its row out);
// - more: [label, value] rows behind a collapsed "More" under the fields (MAP.325), for the facts a reader rarely needs;
// - after: more nodes under the rows (a <details> list, a note);
// - nav: [{label, href | onClick, primary}] -> the NAV buttons
//   (navpick.js): "Start Here" and "End Here", or, while a NAV start or
//   destination is picked, the one button for the end being chosen and no
//   link that would leave the course (NAV.30);
// - bookmark: the entry ☆ saves (bookmarks.js: {kind, value, name, url,
//   sectorId?}), beside the NAV links;
// - view: {href, label}, the "View" button to the item's own page
//   (UX.143), first among the actions, on every map and for every kind
//   (MAP.254);
// - links: [{href, label}], buttons: [{label, onClick}];
// - generate: Generate buttons (generatebuttons.js);
// - hint: a hint paragraph at the bottom.
// hint(text, keep) shows a hint alone, or (keep) under what is there.
//
// The panel is a popup inside the map's window, never a column beside it (MAP.256, MAP.288): the aside carries
// `data-overlay`, optionally naming (by id, space separated) the viewports it can sit in (the System Map has one
// for its diagram and one for its 3D view; the panel moves into the one that is showing), else the one it is in.
// style.css lays it over the map, in the corner farthest from where the reader last pressed (`data-corner`).
// `data-state` is "rest" (nothing picked, not drawn), "selected" or "closed" (the reader pressed Close).
// `reveal()` draws what a map wrote into the panel itself (a measurement) as the popup.
// MAP.347: the card's first view is its name, what kind of thing it is, this many of its values, and its buttons; the rest of what it says
// is behind the Expand button, in a region that scrolls on its own.
export const FIRST_FIELDS = 3;

// What a card's first view says of `spec`: its name, its kind (named by the spec or by its "Type" row, which then is not repeated as a row),
// and its rows with a value. The hover preview (createHoverPreview) takes the same.
export function firstViewOf(spec) {
  const fields = (spec.fields || []).filter(function (field) { return field[1] !== null && field[1] !== undefined && field[1] !== ""; });
  let kind = spec.kind || null;
  if (!kind) {
    const typed = fields.findIndex(function (field) { return field[0] === "Type"; });
    if (typed >= 0) kind = fields.splice(typed, 1)[0][1];
  }
  return { title: spec.title || "Unknown", kind: kind, fields: fields };
}

export function createInfoPanel(panelEl) {
  const hostIds = panelEl && panelEl.hasAttribute("data-overlay") ? (panelEl.dataset.overlay || "").split(/\s+/).filter(Boolean) : null;
  const hosts = !panelEl || !hostIds ? [] : (hostIds.length
    ? hostIds.map(function (id) { return document.getElementById(id); }).filter(Boolean) : [panelEl.parentElement]);
  const overlay = hosts.length > 0;

  // The corner of the map farthest from where the reader last pressed on it, so the popup does not cover what was just
  // picked (top right until the pointer has been used).
  let pressed = null;
  // The card's parts: `rest` holds what is behind Expand (null while the card has none), `expanded` is the reader's choice, kept from one
  // card to the next.
  let rest = null;
  let expandButton = null;
  let expanded = false;
  hosts.forEach(function (host) {
    host.addEventListener("pointerdown", function (event) {
      if (panelEl.contains(event.target)) return;
      const rect = host.getBoundingClientRect();
      pressed = { x: (event.clientX - rect.left) / (rect.width || 1), y: (event.clientY - rect.top) / (rect.height || 1) };
    }, true);
  });
  if (overlay) {
    panelEl.dataset.state = "rest";
    // What happens in the popup is not the map's to handle: a click on its buttons is not a click on the map.
    ["pointerdown", "mousedown", "click", "dblclick", "touchstart"].forEach(function (type) {
      panelEl.addEventListener(type, function (event) { event.stopPropagation(); }, { passive: true });
    });
    // The wheel scrolls the expanded region while it has more to show that way (MAP.295); at the end of its scroll, or when the card
    // has none, the turn goes on to the map under it, so a popup never takes the zoom away.
    panelEl.addEventListener("wheel", function (event) {
      event.stopPropagation();
      const room = rest && !rest.hidden && (event.deltaY < 0 ? rest.scrollTop > 0 : rest.scrollTop + rest.clientHeight < rest.scrollHeight - 1);
      if (room) return;
      const host = panelEl.parentElement;
      const canvas = host && host.querySelector("canvas");
      if (!canvas) return;
      event.preventDefault();
      canvas.dispatchEvent(new WheelEvent("wheel", {
        bubbles: true, cancelable: true, clientX: event.clientX, clientY: event.clientY,
        deltaX: event.deltaX, deltaY: event.deltaY, deltaMode: event.deltaMode,
        ctrlKey: event.ctrlKey, shiftKey: event.shiftKey, altKey: event.altKey, metaKey: event.metaKey,
      }));
    }, { passive: false });
  }

  // The popup sits in the viewport that is showing, in the corner away from the last press.
  function place() {
    const host = (hosts.length > 1 && hosts.find(function (h) { return h.getClientRects().length > 0; })) || hosts[0];
    if (panelEl.parentElement !== host) host.appendChild(panelEl);
    const left = pressed && pressed.x > 0.5 ? "left" : "right";
    const top = pressed && pressed.y > 0.5 ? "top" : (pressed ? "bottom" : "top");
    panelEl.dataset.corner = top + "-" + left;
    panelEl.dataset.state = "selected";
  }

  function closeButton() {
    const close = el("button", "btn btn-small btn-secondary map-info-close", "Close");
    close.type = "button";
    close.addEventListener("click", function () { panelEl.dataset.state = "closed"; });
    return close;
  }

  // Draws what a map wrote into the panel itself (the System Map's measurement) as the popup, with its Close.
  function reveal() {
    if (!overlay) return;
    rest = null;
    expandButton = null;
    place();
    if (!panelEl.querySelector(".map-info-close")) {
      const row = el("p", "page-actions map-info-actions");
      row.appendChild(closeButton());
      panelEl.appendChild(row);
    }
    fit();
  }

  // MAP.347: the card says its first view and nothing more until the reader presses Expand; what is behind it scrolls inside the card
  // (`.map-info-rest`), never the card. An expanded card stands beside the zoom controls at the map's whole height (`data-beside`,
  // style.css) when the window is wide enough for it, else it takes the room above them. MAP.332: a card keeps clear of the zoom
  // controls and the legend in the map's lower corners.
  const BESIDE_REM = 24;
  function fit() {
    if (!overlay || panelEl.dataset.state !== "selected") return;
    delete panelEl.dataset.beside;
    delete panelEl.dataset.expanded;
    panelEl.style.removeProperty("--info-pad-clear");
    const has = !!rest && Array.prototype.some.call(rest.children, function (node) { return !node.hidden; });
    if (expandButton) {
      expandButton.hidden = !has;
      expandButton.textContent = expanded ? "Collapse" : "Expand";
      expandButton.setAttribute("aria-expanded", expanded && has ? "true" : "false");
    }
    if (rest) rest.hidden = !(expanded && has);
    if (!rest || rest.hidden) {
      clearOfPress();
      return;
    }
    panelEl.dataset.expanded = "";
    const view = panelEl.parentElement;
    const rem = parseFloat(getComputedStyle(document.documentElement).fontSize) || 16;
    if (view && view.clientWidth >= BESIDE_REM * rem) {
      // The zoom controls: the direction star's pad and (MAP.336) the log ruler above it, whichever reaches farther in.
      const reach = [".starmap-navpad", ".starmap-ruler"].map(function (selector) {
        const control = view.querySelector(selector);
        return control ? view.getBoundingClientRect().right - control.getBoundingClientRect().left + 8 : 0;
      });
      if (Math.max.apply(null, reach) > 0) panelEl.style.setProperty("--info-pad-clear", Math.ceil(Math.max.apply(null, reach)) + "px");
      panelEl.dataset.beside = "";
    }
    clearOfPress();
  }
  // The card sits on the side of the map away from the press, but a card wide enough, or standing beside the zoom controls, can still lie
  // over the place pressed (a press near the middle of a short map): then it takes the other side, if that leaves the place clear, so a
  // second click of a double click, or the star just picked, is not under the card's own buttons.
  function clearOfPress() {
    const host = panelEl.parentElement;
    const was = panelEl.dataset.corner;
    if (!pressed || !host || !was) return;
    const box = host.getBoundingClientRect();
    const x = box.left + pressed.x * box.width;
    const y = box.top + pressed.y * box.height;
    const covers = function () {
      const r = panelEl.getBoundingClientRect();
      return x >= r.left && x <= r.right && y >= r.top && y <= r.bottom;
    };
    if (!covers()) return;
    panelEl.dataset.corner = was.replace(/(left|right)$/, function (side) { return side === "left" ? "right" : "left"; });
    if (covers()) panelEl.dataset.corner = was;
  }
  // Opening or closing a More changes what the card has to hold: it is fitted again (a `toggle` does not bubble).
  if (overlay) panelEl.addEventListener("toggle", fit, true);
  if (overlay && typeof ResizeObserver !== "undefined") {
    hosts.forEach(function (host) { new ResizeObserver(fit).observe(host); });
  }

  let bookmarkButton = null;
  let bookmarkEntry = null;
  let refreshBookmark = null;

  function bookmarkToggle(entry) {
    bookmarkEntry = entry;
    if (!bookmarkButton) {
      bookmarkButton = el("button", "btn btn-small btn-secondary btn-bookmark map-info-bookmark");
      bookmarkButton.type = "button";
      refreshBookmark = B.toggleButton(bookmarkButton, function () { return bookmarkEntry; });
    } else {
      refreshBookmark();
    }
    return bookmarkButton;
  }

  // Every action in one wrapping row, so the space between two buttons
  // is the same across and down (--btn-gap): the View button, the pick button, the Go to buttons and Expand (MAP.347; the ☆ sits by the name). The
  // other NAV buttons and the page links are a second row at the top of what is behind Expand (`extra`), with the same spacing.
  function actions(spec, secondary) {
    const row = el("p", "page-actions map-info-actions");
    const extra = secondary || row;
    if (spec.view) row.appendChild(pageLink(spec.view.href, spec.view.label, "btn btn-small map-info-view"));
    (spec.nav || []).forEach(function (action) {
      const into = action.primary ? row : extra;
      if (action.href) {
        into.appendChild(pageLink(action.href, action.label, action.primary ? "btn starmap-pick" : "btn btn-small btn-secondary"));
        return;
      }
      const node = el("button", action.primary ? "btn starmap-pick" : "btn btn-small btn-secondary", action.label);
      node.type = "button";
      node.addEventListener("click", action.onClick);
      into.appendChild(node);
    });
    if (spec.bookmark && !overlay) row.appendChild(bookmarkToggle(spec.bookmark));
    (spec.links || []).forEach(function (link) { extra.appendChild(pageLink(link.href, link.label, "btn btn-small btn-secondary")); });
    (spec.buttons || []).forEach(function (button) {
      const node = el("button", "btn btn-small btn-secondary", button.label);
      node.type = "button";
      node.addEventListener("click", button.onClick);
      row.appendChild(node);
    });
    if (overlay) {
      expandButton = el("button", "btn btn-small btn-secondary map-info-expand", "Expand");
      expandButton.type = "button";
      expandButton.setAttribute("aria-expanded", "false");
      expandButton.hidden = true;
      expandButton.addEventListener("click", function () {
        expanded = !expanded;
        fit();
      });
      row.appendChild(expandButton);
    }
    if (row.childNodes.length) panelEl.appendChild(row);
  }

  function show(spec) {
    if (!panelEl) return;
    panelEl.textContent = "";
    bookmarkEntry = null;
    rest = null;
    expandButton = null;
    if (overlay) place();
    // The name, with the ☆ and the Close beside it so the buttons below keep to one row.
    const view = firstViewOf(spec);
    const head = el("div", "map-info-head");
    head.appendChild(el("h3", null, view.title));
    if (overlay && spec.bookmark) head.appendChild(bookmarkToggle(spec.bookmark));
    if (overlay) head.appendChild(closeButton());
    panelEl.appendChild(head);
    // The kind: named by the card or by its "Type" row, which then is not repeated as a row.
    const fields = view.fields;
    if (view.kind) panelEl.appendChild(el("p", "map-info-kind", String(view.kind)));
    const first = overlay ? fields.slice(0, FIRST_FIELDS) : fields;
    if (first.length) {
      const dl = document.createElement("dl");
      first.forEach(function (field) { addField(dl, field[0], field[1]); });
      panelEl.appendChild(dl);
    }
    // What the card says beyond its first view sits in one region, hidden until Expand (MAP.347); a card of a map with no popup shows it all.
    const behind = overlay ? el("div", "map-info-rest") : panelEl;
    const extra = overlay ? el("p", "page-actions map-info-actions map-info-secondary") : null;
    if (overlay) {
      behind.hidden = true;
      rest = behind;
    }
    actions(spec, extra);
    if (overlay && extra.childNodes.length) behind.appendChild(extra);
    if (overlay && fields.length > FIRST_FIELDS) {
      const dl = document.createElement("dl");
      fields.slice(FIRST_FIELDS).forEach(function (field) { addField(dl, field[0], field[1]); });
      behind.appendChild(dl);
    }
    // MAP.325: the facts a reader rarely needs sit behind a More that starts collapsed.
    const more = (spec.more || []).filter(function (field) { return field[1] !== null && field[1] !== undefined && field[1] !== ""; });
    if (more.length) {
      const fold = el("details", "map-info-more");
      fold.appendChild(el("summary", null, "More"));
      const list = document.createElement("dl");
      more.forEach(function (field) { addField(list, field[0], field[1]); });
      fold.appendChild(list);
      behind.appendChild(fold);
    }
    (spec.after || []).forEach(function (node) { behind.appendChild(node); });
    if (spec.generate) behind.appendChild(spec.generate);
    if (spec.hint) behind.appendChild(el("p", "hint", spec.hint));
    if (overlay) panelEl.appendChild(behind);
    fit();
  }

  function hint(text, keep) {
    if (!panelEl) return;
    if (!keep) {
      panelEl.textContent = "";
      bookmarkEntry = null;
      rest = null;
      expandButton = null;
      if (overlay) panelEl.dataset.state = "rest";
    }
    (keep && rest ? rest : panelEl).appendChild(el("p", "hint", text));
    fit();
  }

  return { show: show, hint: hint, reveal: reveal, refit: fit, element: panelEl };
}

// --- Hover preview -----------------------------------------------------------------

// MAP.348: the small card that a hovered star, system, region or cloud shows: the card's first view (its name, kind and up to three
// values, no buttons), in a corner of the map `host`. It is never under the pointer or the hovered thing, never over the open card or the
// map's controls, never takes a pointer event, and shows after the pointer has rested on the same thing for `delay` ms.
// show(key, spec, avoid): `key` names the thing (the same key keeps the preview as it is), `spec` is the card's spec,
// `avoid()` lists the client points {x, y} it must stay clear of (the pointer, the thing's spot). hide() puts it away.
const PREVIEW_CLEAR_SELECTORS = [
  ".starmap-info[data-state=\"selected\"]", ".map-corner-stack", ".map-layer-legend", ".starmap-navpad", ".starmap-ruler", ".starmap-scale",
  ".galaxy-menu-panel", ".map-tooltip:not([hidden])",
];
const PREVIEW_EDGE = 8;
const PREVIEW_CLEARANCE = 28;

export function createHoverPreview(host, options) {
  if (!host) return { show: function () {}, hide: function () {}, element: null };
  const delay = options && options.delay != null ? options.delay : 180;
  const node = el("aside", "map-hover-preview");
  node.hidden = true;
  node.setAttribute("aria-hidden", "true");    // the tooltip carries the words for assistive technology
  host.appendChild(node);
  let shownKey = null;
  let timer = null;

  function hide() {
    if (timer) window.clearTimeout(timer);
    timer = null;
    shownKey = null;
    node.hidden = true;
  }

  function overlaps(a, b) {
    return a.left < b.right && a.right > b.left && a.top < b.bottom && a.bottom > b.top;
  }

  // The first corner, nearest the top and farthest from the points to avoid, that nothing else of the map is under; false when none is.
  function place(points) {
    const box = host.getBoundingClientRect();
    const own = node.getBoundingClientRect();
    const width = own.width;
    const height = own.height;
    const blocked = [];
    PREVIEW_CLEAR_SELECTORS.forEach(function (selector) {
      host.querySelectorAll(selector).forEach(function (other) {
        if (other === node || other.hidden || !other.getClientRects().length) return;
        blocked.push(other.getBoundingClientRect());
      });
    });
    points.forEach(function (point) {
      blocked.push({ left: point.x - PREVIEW_CLEARANCE, right: point.x + PREVIEW_CLEARANCE, top: point.y - PREVIEW_CLEARANCE, bottom: point.y + PREVIEW_CLEARANCE });
    });
    const middle = points.length ? points.reduce(function (sum, point) { return sum + point.x; }, 0) / points.length : box.left + box.width / 2;
    const farLeft = Math.abs(box.left + PREVIEW_EDGE + width / 2 - middle) >= Math.abs(box.right - PREVIEW_EDGE - width / 2 - middle);
    const sides = farLeft ? ["left", "right"] : ["right", "left"];
    const corners = [["top", sides[0]], ["top", sides[1]], ["bottom", sides[0]], ["bottom", sides[1]]];
    for (let i = 0; i < corners.length; i++) {
      const left = corners[i][1] === "left" ? box.left + PREVIEW_EDGE : box.right - PREVIEW_EDGE - width;
      const top = corners[i][0] === "top" ? box.top + PREVIEW_EDGE : box.bottom - PREVIEW_EDGE - height;
      const rect = { left: left, right: left + width, top: top, bottom: top + height };
      if (rect.left < box.left || rect.right > box.right || rect.top < box.top || rect.bottom > box.bottom) continue;
      if (blocked.some(function (other) { return overlaps(rect, other); })) continue;
      node.style.left = (left - box.left) + "px";
      node.style.top = (top - box.top) + "px";
      node.dataset.corner = corners[i].join("-");
      return true;
    }
    return false;
  }

  function draw(spec, avoid) {
    const view = firstViewOf(spec);
    node.textContent = "";
    node.appendChild(el("h3", null, view.title));
    if (view.kind) node.appendChild(el("p", "map-info-kind", String(view.kind)));
    const rows = view.fields.slice(0, FIRST_FIELDS);
    if (rows.length) {
      const dl = document.createElement("dl");
      rows.forEach(function (field) { addField(dl, field[0], field[1]); });
      node.appendChild(dl);
    }
    node.hidden = false;
    node.style.visibility = "hidden";    // measured in place before it is seen
    const fits = place(avoid());
    node.style.visibility = "";
    if (!fits) node.hidden = true;
  }

  function show(key, spec, avoid) {
    if (key === shownKey) return;
    hide();
    shownKey = key;
    timer = window.setTimeout(function () {
      timer = null;
      draw(spec(), avoid);
    }, delay);
  }

  return { show: show, hide: hide, element: node };
}

// The one info panel for `panelEl` (so its ☆ button is made once), or
// null when there is no such element.
const panels = new WeakMap();
export function infoPanelOf(panelEl) {
  if (!panelEl) return null;
  if (!panels.has(panelEl)) panels.set(panelEl, createInfoPanel(panelEl));
  return panels.get(panelEl);
}

// The bookmark entry for a system or phenomenon: `endpoint` is its NAV
// endpoint ("system:12", "nebula:3"), whose kind is the part before the
// colon.
export function endpointBookmark(endpoint, name, url) {
  if (!endpoint || endpoint.indexOf(":") < 1) return null;
  return { kind: endpoint.split(":")[0], value: endpoint, name: name || endpoint, url: url || null };
}
