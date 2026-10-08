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
//   means (null: nothing);
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
      return { entry: found.entry, distance: camera.position.distanceTo(at), px: found.px };
    }
    if (layer.meshes) {
      const hits = raycaster.intersectObjects(layer.meshes(), false);
      for (let n = 0; n < hits.length; n++) {
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
// - after: more nodes under the rows (a <details> list, a note);
// - nav: {from, to} -> "Nav from here" and "Nav to here" links, or, while
//   a NAV start or destination is picked, {pick, pickLabel} -> one pick
//   button and no link that would leave the course (NAV.30);
// - bookmark: the entry ☆ saves (bookmarks.js: {kind, value, name, url,
//   sectorId?}), beside the NAV links;
// - links: [{href, label}], buttons: [{label, onClick}];
// - generate: Generate buttons (generatebuttons.js);
// - hint: a hint paragraph at the bottom.
// hint(text, keep) shows a hint alone, or (keep) under what is there.
export function createInfoPanel(panelEl) {
  let bookmarkButton = null;
  let bookmarkEntry = null;
  let refreshBookmark = null;

  function bookmarkToggle(entry) {
    bookmarkEntry = entry;
    if (!bookmarkButton) {
      bookmarkButton = el("button", "btn btn-small btn-bookmark map-info-bookmark");
      bookmarkButton.type = "button";
      refreshBookmark = B.toggleButton(bookmarkButton, function () { return bookmarkEntry; }, false);
    } else {
      refreshBookmark();
    }
    return bookmarkButton;
  }

  // Every action in one wrapping row, so the space between two buttons
  // is the same across and down (--btn-gap): the pick button or the NAV
  // links, the ☆, page links and buttons.
  function actions(spec) {
    const row = el("p", "page-actions map-info-actions");
    const nav = spec.nav;
    if (nav && nav.pick) {
      row.appendChild(pageLink(nav.pick, nav.pickLabel, "btn starmap-pick"));
    } else if (nav) {
      [["from", "Nav from here"], ["to", "Nav to here"]].forEach(function (pair) {
        if (nav[pair[0]]) row.appendChild(pageLink(nav[pair[0]], pair[1], "btn btn-small"));
      });
    }
    if (spec.bookmark) row.appendChild(bookmarkToggle(spec.bookmark));
    (spec.links || []).forEach(function (link) { row.appendChild(pageLink(link.href, link.label, "btn btn-small")); });
    (spec.buttons || []).forEach(function (button) {
      const node = el("button", "btn btn-small", button.label);
      node.type = "button";
      node.addEventListener("click", button.onClick);
      row.appendChild(node);
    });
    if (row.childNodes.length) panelEl.appendChild(row);
  }

  function show(spec) {
    if (!panelEl) return;
    panelEl.textContent = "";
    bookmarkEntry = null;
    panelEl.appendChild(el("h3", null, spec.title || "Unknown"));
    if (spec.fields && spec.fields.length) {
      const dl = document.createElement("dl");
      spec.fields.forEach(function (field) { addField(dl, field[0], field[1]); });
      panelEl.appendChild(dl);
    }
    (spec.after || []).forEach(function (node) { panelEl.appendChild(node); });
    actions(spec);
    if (spec.generate) panelEl.appendChild(spec.generate);
    if (spec.hint) hint(spec.hint, true);
  }

  function hint(text, keep) {
    if (!panelEl) return;
    if (!keep) {
      panelEl.textContent = "";
      bookmarkEntry = null;
    }
    panelEl.appendChild(el("p", "hint", text));
  }

  return { show: show, hint: hint, element: panelEl };
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
