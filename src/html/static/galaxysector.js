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
// - fetchUnchartedScene(at) (optional): a Promise of the scene JSON of the
//   cell `at` ({ring, layer, slot}) nothing was generated in, which holds
//   what the scatters left there (MAP.162); without it such a cell can't
//   be opened;
// - showInfo(spec): the info panel (mappick.js's spec);
// - closeness(): how near the view is to the sector's fit (1 at the fit,
//   2.5 at the closest), which grows the points a little (default 1);
// - hideCell(bounds or null): the galaxy's own tile stars and clouds are
//   left out of this cell (the sector's scene draws them), null to show
//   them again;
// - viewport (optional): the map's viewport element, where the
//   screen-reader list of the open sector's entries goes;
// - opened() / closed() (optional): called when a sector's scene is added
//   to or removed from the map;
// - openSystem(entry) (optional): opens a star's system in place (MAP.125);
//   the star's panel offers it when given;
// - kindsChanged() (optional): called when a kind of object is shown or
//   hidden other than by setKindHidden (a selection showing its kind).

const VERSION_QUERY = new URL(import.meta.url).search;
const { createRing } = await import(`./mappick.js${VERSION_QUERY}`);
const {
  buildSectorScene, entryLabel, infoSpec, kindOf, KINDS, STAR_CLASSES, starClassOf, POINT_CLOSE_GROWTH, tooltipText,
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

export const DEFAULT_HIDDEN_KINDS = ["roguePlanet"];

export function createSectorStage(host) {
  const selectionRing = createRing(host.scene, host.camera, host.canvasEl, host.accentColor, { depthTest: false, renderOrder: 5 });
  const hoverRing = createRing(host.scene, host.camera, host.canvasEl, host.accentColor, { depthTest: false, renderOrder: 5, opacity: 0.45 });
  // {id, data, sector (the scene), layers, selected, hovered, center, halfEdge}
  let open = null;
  let token = 0;
  // "Mark rogue planets" (MAP.46): kept across sectors opened. On to begin
  // with (MAP.137), so a rogue planet is easy to find once it is shown.
  let roguesMarked = true;
  // The kinds of object left off the map (MAP.79), kept across sectors:
  // rogue planets to begin with (MAP.137).
  const hiddenKinds = new Set(DEFAULT_HIDDEN_KINDS);
  const markedKinds = new Set();
  // MAP.123: star classes left off and the dimmest star shown (L☉).
  const hiddenClasses = new Set();
  let minLuminosity = 0;

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
          state.hoveredEntry = null;
          return;
        }
        state.hoveredEntry = entry;
        state.hovered = ringAround(hoverRing, entry);
      },
      select: function (entry) {
        // The selected one again clears it: a nebula can cover the whole
        // sector, so there is nowhere else to click (MAP.113).
        if (state.selectedEntry === entry) {
          deselect();
          return;
        }
        const spec = infoSpec(entry, state.data, host.navPick);
        // A system opens in the map itself (MAP.125), unless a NAV end is being picked.
        if (host.openSystem && entry.endpoint && entry.href && !(host.navPick && host.navPick.active())) {
          spec.buttons = (spec.buttons || []).concat([{ label: "Open system here", onClick: function () { host.openSystem(entry); } }]);
        }
        host.showInfo(spec);
        state.selectedEntry = entry;
        state.selected = ringAround(selectionRing, entry);
      },
    });
  }

  // Clears the selection; true when there was one.
  function deselect() {
    if (!open || !open.selectedEntry) return false;
    selectionRing.hide();
    open.selected = open.selectedEntry = null;
    if (host.deselected) host.deselected();
    return true;
  }

  function close() {
    token += 1;
    if (!open) return;
    open.layers.forEach(host.picker.removeLayer);
    if (open.list && open.list.parentNode) open.list.parentNode.removeChild(open.list);
    if (open.badge && open.badge.parentNode) open.badge.parentNode.removeChild(open.badge);
    open.sector.dispose();
    selectionRing.hide();
    hoverRing.hide();
    host.hideCell(null);
    open = null;
    if (host.closed) host.closed();
  }

  // A selection of something hidden shows its kind again first.
  function selectEntry(entry) {
    if (!open || !entry) return;
    // Already selected: it stays so (a click on it clears it).
    if (open.selectedEntry === entry) return;
    const kind = kindOf(entry);
    if (hiddenKinds.has(kind)) {
      setKindHidden(kind, false);
      if (host.kindsChanged) host.kindsChanged();
    }
    if (kind === "star" && open.sector.hides(entry)) {
      // Picked by name (a Contents row): its class and the floor give way.
      setStarClassHidden(starClassOf(entry), false);
      if ((entry.luminositySol || 0) < minLuminosity) setMinLuminosity(0);
      if (host.kindsChanged) host.kindsChanged();
    }
    open.layers[0].select(entry);
  }

  function setKindMarked(kind, marked) {
    if (marked) markedKinds.add(kind);
    else markedKinds.delete(kind);
    if (open) open.sector.setKindMarked(kind, marked);
  }

  function setKindHidden(kind, hidden) {
    if (hidden) hiddenKinds.add(kind);
    else hiddenKinds.delete(kind);
    if (!open) return;
    open.sector.setKindHidden(kind, hidden);
    if (!hidden) {
      showListed();
      return;
    }
    // What was selected or hovered goes with its kind.
    if (open.selectedEntry && kindOf(open.selectedEntry) === kind) {
      selectionRing.hide();
      open.selected = open.selectedEntry = null;
    }
    if (open.hoveredEntry && kindOf(open.hoveredEntry) === kind) {
      hoverRing.hide();
      open.hovered = open.hoveredEntry = null;
    }
    showListed();
  }

  // The screen-reader list offers only what is on the map.
  function showListed() {
    if (!open || !open.listed) return;
    open.listed.forEach(function (pair) { pair[1].hidden = open.sector.hides(pair[0]); });
  }

  // What was selected or hovered goes when it is hidden.
  function dropHidden() {
    if (open.selectedEntry && open.sector.hides(open.selectedEntry)) {
      selectionRing.hide();
      open.selected = open.selectedEntry = null;
    }
    if (open.hoveredEntry && open.sector.hides(open.hoveredEntry)) {
      hoverRing.hide();
      open.hovered = open.hoveredEntry = null;
    }
    showListed();
  }

  // MAP.123: hides or shows the stars of one class.
  function setStarClassHidden(starClass, hidden) {
    if (hidden) hiddenClasses.add(starClass);
    else hiddenClasses.delete(starClass);
    if (!open) return;
    open.sector.setStarClassHidden(starClass, hidden);
    dropHidden();
  }

  // MAP.123: shows only stars at least this luminous (L☉; 0 shows all).
  function setMinLuminosity(value) {
    minLuminosity = value > 0 ? value : 0;
    if (!open) return;
    open.sector.setMinLuminosity(minLuminosity);
    dropHidden();
  }

  // A visually hidden button per entry, for the keyboard and screen
  // readers (a canvas has no focusable parts), in the scene's order.
  function entryList(entries) {
    if (!host.viewport) return null;
    const list = document.createElement("ul");
    list.className = "starmap-sr-list sr-only";
    list.listed = [];
    entries.forEach(function (entry) {
      const item = document.createElement("li");
      list.listed.push([entry, item]);
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
  //
  // `at` ({ring, layer, slot}) opens the cell nothing was generated in
  // instead (MAP.162, `id` is then null): it opens only when the scatters
  // left something in it, and carries an "Uncharted" badge in its frame.
  function openSector(id, bounds, at) {
    const mine = ++token;
    const fetched = at ? host.fetchUnchartedScene(at) : host.fetchScene(id);
    return fetched.then(function (data) {
      if (mine !== token || !data || !data.centerPc || !(data.halfEdgePc > 0)) return null;
      if (data.uncharted && !(data.stars.length || data.clouds.length)) return null;
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
      const state = { id: id, data: data, selected: null, hovered: null, uncharted: !!data.uncharted };
      state.center = data.centerPc;
      state.halfEdge = data.halfEdgePc;
      state.sector = buildSectorScene(data, {
        origin: data.centerPc, unit: data.halfEdgePc / half, flipY: true,
        accentColor: host.accentColor, lightBackground: host.lightBackground, pixelRatio: host.pixelRatio(),
        sizeScale: pointSizeScale,
        // How many nebulae are drawn from their shape now (read by the browser tests).
        onShape: function (count) { host.canvasEl.dataset.sectorNebulaMeshes = String(count); },
      });
      host.scene.add(state.sector.group);
      state.sector.setRoguesMarked(roguesMarked);
      hiddenKinds.forEach(function (kind) { state.sector.setKindHidden(kind, true); });
      markedKinds.forEach(function (kind) { state.sector.setKindMarked(kind, true); });
      hiddenClasses.forEach(function (c) { state.sector.setStarClassHidden(c, true); });
      if (minLuminosity > 0) state.sector.setMinLuminosity(minLuminosity);
      state.layers = state.sector.layers.map(function (layer) { return host.picker.addLayer(layerFor(layer, state)); });
      state.list = entryList(state.sector.entries);
      state.listed = state.list ? state.list.listed : null;
      open = state;
      if (state.uncharted) state.badge = unchartedBadge(data);
      showListed();
      host.hideCell(bounds || null);
      if (host.opened) host.opened();
      return { center: state.center, halfEdge: state.halfEdge, fitPoints: cubeCorners(state.center, state.halfEdge) };
    });
  }

  // The "Uncharted" mark on the open sector's frame (the map's viewport),
  // naming the sector and what was left in it.
  function unchartedBadge(data) {
    if (!host.viewport) return null;
    const badge = document.createElement("div");
    badge.className = "sector-uncharted-badge";
    badge.setAttribute("role", "note");
    const count = data.stars.length + data.clouds.length;
    badge.textContent = "Uncharted sector " + (data.designation || "") + ": not generated yet, "
      + count + (count === 1 ? " object" : " objects") + " placed";
    host.viewport.appendChild(badge);
    return badge;
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
    // Whether the open sector is one nothing was generated in (MAP.162), and
    // what the scatters left in it: {stars, objects, designation}.
    isUncharted: function () { return !!(open && open.uncharted); },
    unchartedSummary: function () {
      return open && open.uncharted
        ? { stars: open.data.stars.length, objects: open.data.clouds.length, designation: open.data.designation }
        : null;
    },
    select: selectEntry,
    deselect: deselect,
    // The open sector's entry with this key ("rogue_planet:12"), or null.
    entryByKey: function (key) { return open ? open.sector.entryByKey.get(key) || null : null; },
    setRoguesMarked: function (on) {
      roguesMarked = on;
      if (open) open.sector.setRoguesMarked(on);
    },
    roguesMarked: function () { return roguesMarked; },
    // The kinds of object in the open sector, [{kind, label, count}] (MAP.79).
    kinds: function () { return open ? open.sector.kinds() : []; },
    setKindHidden: setKindHidden,
    kindHidden: function (kind) { return hiddenKinds.has(kind); },
    hiddenKinds: function () { return Array.from(hiddenKinds); },
    setKindMarked: setKindMarked,
    kindMarked: function (kind) { return markedKinds.has(kind); },
    markedKinds: function () { return Array.from(markedKinds); },
    // The kinds named by a URL's `mark` (unknown names ignored), before any sector opens.
    setMarkedKinds: function (kinds) {
      let changed = false;
      KINDS.forEach(function (pair) {
        const want = kinds.indexOf(pair[0]) >= 0;
        if (want === markedKinds.has(pair[0])) return;
        setKindMarked(pair[0], want);
        changed = true;
      });
      if (changed && host.kindsChanged) host.kindsChanged();
    },
    setStarClassHidden: setStarClassHidden,
    starClassHidden: function (c) { return hiddenClasses.has(c); },
    hiddenClasses: function () { return Array.from(hiddenClasses); },
    starClasses: function () { return open ? open.sector.starClasses() : []; },
    luminosityRange: function () { return open ? open.sector.luminosityRange() : null; },
    // The classes named by a URL's `stars` (unknown names ignored), before any sector opens.
    setHiddenClasses: function (classes) {
      let changed = false;
      STAR_CLASSES.forEach(function (c) {
        const want = classes.indexOf(c) >= 0;
        if (want === hiddenClasses.has(c)) return;
        setStarClassHidden(c, want);
        changed = true;
      });
      if (changed && host.kindsChanged) host.kindsChanged();
    },
    setMinLuminosity: setMinLuminosity,
    minLuminosity: function () { return minLuminosity; },
    // The kinds named by a URL's `hide` (unknown names ignored), before any sector opens.
    setHiddenKinds: function (kinds) {
      let changed = false;
      KINDS.forEach(function (pair) {
        const want = kinds.indexOf(pair[0]) >= 0;
        if (want === hiddenKinds.has(pair[0])) return;
        setKindHidden(pair[0], want);
        changed = true;
      });
      if (changed && host.kindsChanged) host.kindsChanged();
    },
    entries: function () { return open ? open.sector.entries : []; },
    // The entry picked in the open sector, if any.
    selectedEntry: function () { return open ? open.selectedEntry : null; },
  };
}
