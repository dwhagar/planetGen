// html/static/galaxysector.js
//
// The Galaxy Map's last drill-down stage (MAP.66): one sector opened in
// its place in the galaxy. The sector's own scene (sectorscene.js, the same
// stars, phenomena and bodies the sector page draws) is added to the
// Galaxy Map's scene at the sector's position, its picker layers join the
// map's picker (mappick.js), and a click, hover or tooltip on one of its
// entries goes through the same info panel and rings as anywhere else.
// galaxystageview.js decides when a sector is open and flies the camera to
// it; this file owns what is loaded.
//
// host:
// - THREE, scene, camera, canvasEl, picker (the map's);
// - accentColor, lightBackground, pixelRatio();
// - fetchScene(id): a Promise of GET /sector/<id>/scene's JSON;
// - showInfo(spec): the info panel (mappick.js's spec);
// - closeness(): how near the view is to the sector's fit (1 at the fit,
//   2.5 at the closest), which grows the points a little (default 1);
// - hideCell(bounds or null): the galaxy's own tile stars and clouds are
//   left out of this cell (the sector's scene draws them), null to show
//   them again;
// - viewport (optional): the map's viewport element, where the
//   screen-reader list of the open sector's entries goes;
// - opened() / closed() (optional): called when a sector's scene is added
//   to or removed from the map.

const VERSION_QUERY = new URL(import.meta.url).search;
const { createRing } = await import(`./mappick.js${VERSION_QUERY}`);
const {
  buildSectorScene, entryLabel, infoSpec, POINT_CLOSE_GROWTH, tooltipText,
} = await import(`./sectorscene.js${VERSION_QUERY}`);

// The corners of the sector's cell as a flat [x, y, z, ...] list, for
// fitting the camera round it: a cube `halfEdge` either side of `center`.
export function cubeCorners(center, halfEdge) {
  const points = [];
  [-1, 1].forEach(function (sx) {
    [-1, 1].forEach(function (sy) {
      [-1, 1].forEach(function (sz) {
        points.push(center[0] + sx * halfEdge, center[1] + sy * halfEdge, center[2] + sz * halfEdge);
      });
    });
  });
  return points;
}

export function createSectorStage(host) {
  const selectionRing = createRing(host.scene, host.camera, host.canvasEl, host.accentColor, { depthTest: false, renderOrder: 5 });
  const hoverRing = createRing(host.scene, host.camera, host.canvasEl, host.accentColor, { depthTest: false, renderOrder: 5, opacity: 0.45 });
  // {id, data, sector (the scene), layers, selected, hovered, center, halfEdge}
  let open = null;
  let token = 0;
  // "Mark rogue planets" (MAP.46): kept across sectors opened.
  let roguesMarked = false;

  // How much the points have grown with the view closing in (as the
  // Sector Map's own zoom does): their own size at the fit, POINT_CLOSE_GROWTH
  // times it at the closest zoom.
  function pointSizeScale() {
    const near = host.closeness ? host.closeness() : 1;
    const share = Math.max(0, Math.min(1, (near - 1) / 1.5));
    return 1 + (POINT_CLOSE_GROWTH - 1) * share;
  }

  function ringAround(ring, entry) {
    ring.at(entry.x, entry.y, entry.z, open.sector.ringSize(entry));
    return entry.light ? entry : null;
  }

  function layerFor(layer, state) {
    return Object.assign({}, layer, {
      tooltip: tooltipText,
      hover: function (entry) {
        if (!entry) {
          hoverRing.hide();
          state.hovered = null;
          return;
        }
        state.hovered = ringAround(hoverRing, entry);
      },
      select: function (entry) {
        host.showInfo(infoSpec(entry, state.data));
        state.selected = ringAround(selectionRing, entry);
      },
    });
  }

  function close() {
    token += 1;
    if (!open) return;
    open.layers.forEach(host.picker.removeLayer);
    if (open.list && open.list.parentNode) open.list.parentNode.removeChild(open.list);
    open.sector.dispose();
    selectionRing.hide();
    hoverRing.hide();
    host.hideCell(null);
    open = null;
    if (host.closed) host.closed();
  }

  function selectEntry(entry) {
    if (open && entry) open.layers[0].select(entry);
  }

  // A visually hidden button per entry, for the keyboard and screen
  // readers (a canvas has no focusable parts), in the scene's order.
  function entryList(entries) {
    if (!host.viewport) return null;
    const list = document.createElement("ul");
    list.className = "starmap-sr-list sr-only";
    entries.forEach(function (entry) {
      const item = document.createElement("li");
      const button = document.createElement("button");
      button.type = "button";
      button.textContent = entryLabel(entry);
      button.addEventListener("click", function () { selectEntry(entry); });
      item.appendChild(button);
      list.appendChild(item);
    });
    host.viewport.appendChild(list);
    return list;
  }

  // Loads sector `id` (its cell in the galaxy's polar frame, `bounds`, as
  // galaxystages.js's blocks have it) and puts it in the scene: resolves with
  // {center, halfEdge, fitPoints} (the galaxy's frame, parsecs), or
  // null when it is replaced meanwhile, can't be placed, or fails.
  function openSector(id, bounds) {
    const mine = ++token;
    return host.fetchScene(id).then(function (data) {
      if (mine !== token || !data || !data.centerPc || !(data.halfEdgePc > 0)) return null;
      if (open) {
        const keep = token;
        close();
        token = keep;
      }
      // The galaxy shows the neighbors and the cell's frame itself.
      data.neighbors = [];
      data.outline = null;
      data.compass = null;
      const half = data.sceneHalfPx || 160;
      const state = { id: id, data: data, selected: null, hovered: null };
      state.center = data.centerPc;
      state.halfEdge = data.halfEdgePc;
      state.sector = buildSectorScene(data, {
        origin: data.centerPc, unit: data.halfEdgePc / half, flipY: true,
        accentColor: host.accentColor, lightBackground: host.lightBackground, pixelRatio: host.pixelRatio(),
        sizeScale: pointSizeScale,
      });
      host.scene.add(state.sector.group);
      state.sector.setRoguesMarked(roguesMarked);
      state.layers = state.sector.layers.map(function (layer) { return host.picker.addLayer(layerFor(layer, state)); });
      state.list = entryList(state.sector.entries);
      open = state;
      host.hideCell(bounds || null);
      if (host.opened) host.opened();
      return { center: state.center, halfEdge: state.halfEdge, fitPoints: cubeCorners(state.center, state.halfEdge) };
    });
  }

  // Every frame: the points' growth and the rings' size.
  function update() {
    if (!open) return;
    open.sector.update();
    [[selectionRing, open.selected], [hoverRing, open.hovered]].forEach(function (pair) {
      const size = pair[1] && open.sector.ringSize(pair[1]);
      if (size && size.px != null) pair[0].setPx(size.px);
    });
    selectionRing.update();
    hoverRing.update();
  }

  return {
    open: openSector,
    close: close,
    update: update,
    isOpen: function () { return !!open; },
    sectorId: function () { return open ? open.id : null; },
    select: selectEntry,
    // The open sector's entry with this key ("rogue_planet:12"), or null.
    entryByKey: function (key) { return open ? open.sector.entryByKey.get(key) || null : null; },
    setRoguesMarked: function (on) {
      roguesMarked = on;
      if (open) open.sector.setRoguesMarked(on);
    },
    roguesMarked: function () { return roguesMarked; },
    entries: function () { return open ? open.sector.entries : []; },
  };
}
