// html/static/galaxysystem.js
//
// A star system opened in place on the Galaxy Map (MAP.125): below the
// sector stage in the zoom, the system's own scene (systemview3d.js, the
// same bodies and orbits the system page draws) is added to the map's scene
// at the system's position, its bodies join the map's picker (mappick.js),
// and a click, hover or tooltip on one goes through the map's info panel and
// rings. galaxystageview.js decides when a system is open and flies the
// camera there; this file owns what is loaded.
//
// The scene is built in its own units (systemscale.js) inside a group
// positioned at the system and scaled so the system's real extent is what the
// scene's units add up to in parsecs. Every body hangs from its parent with
// small local coordinates, so the camera-relative maths three.js does keeps a
// moon steady next to a planet at galactic distances (floating origin).
//
// host:
// - THREE, scene, camera, canvasEl, picker (the map's); accentColor;
// - fetchScene(href): a Promise of GET /system/<id>/scene's JSON;
// - showInfo(spec): the info panel (mappick.js's spec);
// - flyTo(entry): the map's camera flies to a body of the open system;
// - viewport (optional): the map's viewport element, where the
//   screen-reader list of the open system's bodies goes;
// - picked(ref) (optional): a body was picked, or (null) the pick cleared;
// - opened() / closed() (optional).

const VERSION_QUERY = new URL(import.meta.url).search;
const { createRing } = await import(`./mappick.js${VERSION_QUERY}`);
const { buildSystemScene } = await import(`./systemview3d.js${VERSION_QUERY}`);
const { createLayout, layoutPositions, MODE_COMPRESSED, MODE_TRUE, UNIT_KM } = await import(`./systemscale.js${VERSION_QUERY}`);
const { createOrbitClock } = await import(`./orbitclock.js${VERSION_QUERY}`);
const { relativeAt } = await import(`./orbitpositions.js${VERSION_QUERY}`);
const { formatDistanceM } = await import(`./distance.js${VERSION_QUERY}`);

const KM_PER_PC = 3.0856775814913673e13;
const AU_KM = 149597870.7;
const PICK_REACH_PX = 16;
// How far the compressed layout reaches (its outer orbit), at least, in units.
const MIN_EXTENT_KM = 5 * AU_KM;

// The corners of a cube `half` either side of `center`, for fitting the camera.
export function cubeAround(center, half) {
  const points = [];
  [-1, 1].forEach(function (sx) {
    [-1, 1].forEach(function (sy) {
      [-1, 1].forEach(function (sz) {
        points.push(center[0] + sx * half, center[1] + sy * half, center[2] + sz * half);
      });
    });
  });
  return points;
}

// The real size of a system: its outermost planet or belt, or its
// heliopause, in km (at least MIN_EXTENT_KM).
export function realExtentKm(scene) {
  let far = MIN_EXTENT_KM;
  scene.planets.forEach(function (p) { far = Math.max(far, p.orbit.distance_km * 1.3); });
  scene.belts.forEach(function (b) { far = Math.max(far, b.outer_km * 1.1); });
  scene.stars.forEach(function (s) { if (s.orbit) far = Math.max(far, s.orbit.distance_km * 1.5); });
  return far;
}

export function createSystemStage(host) {
  const THREE = host.THREE;
  const selectionRing = createRing(host.scene, host.camera, host.canvasEl, host.accentColor, { depthTest: false, renderOrder: 5 });
  const hoverRing = createRing(host.scene, host.camera, host.canvasEl, host.accentColor, { depthTest: false, renderOrder: 5, opacity: 0.45 });
  // {href, id, name, data, mode, center, built, wrapper, world (entries in
  // the map's world), unitPc, radiusPc, clock, layer, selected, hovered}
  let open = null;
  let token = 0;

  // The scale of the open system, set by its mode: true scale is one unit a
  // million km; compressed fits the layout's reach to the system's real extent.
  function measure(data, mode) {
    const layout = createLayout(data, mode);
    if (mode === MODE_TRUE) {
      const unitPc = UNIT_KM / KM_PER_PC;
      return { layout: layout, unitPc: unitPc, radiusPc: (realExtentKm(data) / UNIT_KM) * unitPc };
    }
    // The layout's reach in units: the furthest body drawn at the epoch.
    const placed = relativeAt(data, 0);
    let reach = 1;
    const probe = layoutPositions(layout, placed);
    Object.keys(probe).forEach(function (ref) { reach = Math.max(reach, Math.hypot(...probe[ref])); });
    const radiusPc = realExtentKm(data) / KM_PER_PC;
    return { layout: layout, unitPc: radiusPc / (reach * 1.15), radiusPc: radiusPc };
  }

  function rebuild(state) {
    if (state.built) {
      state.wrapper.remove(state.built.group);
      state.built.dispose();
    }
    const m = measure(state.data, state.mode);
    state.unitPc = m.unitPc;
    state.radiusPc = m.radiusPc;
    state.built = buildSystemScene(state.data, m.layout, { years: state.clock.years(), glow: false });
    state.wrapper.scale.setScalar(m.unitPc);
    state.wrapper.add(state.built.group);
    state.world = state.built.entries.map(function (e) {
      return { ref: e.ref, key: e.ref, name: e.name, kind: e.kind, body: e.body, local: e, x: 0, y: 0, z: 0, r: 0 };
    });
    syncWorld(state);
  }

  function syncWorld(state) {
    state.world.forEach(function (w) {
      w.x = state.center[0] + state.unitPc * w.local.x;
      w.y = state.center[1] + state.unitPc * w.local.y;
      w.z = state.center[2] + state.unitPc * w.local.z;
      w.r = state.unitPc * w.local.radius;
    });
  }

  function infoSpecFor(entry, state) {
    const body = entry.body;
    const fields = [["Kind", entry.kind]];
    if (body.planet_class) fields.push(["Class", body.planet_class]);
    if (body.star_type) fields.push(["Star type", body.star_type]);
    if (body.radius_km) fields.push(["Radius", formatDistanceM(body.radius_km * 1000)]);
    if (body.orbit && body.orbit.distance_km) {
      fields.push(["Orbit", formatDistanceM(body.orbit.distance_km * 1000)]);
      if (body.orbit.period_years) fields.push(["Period", body.orbit.period_years.toPrecision(3) + " years"]);
    }
    return {
      title: entry.name,
      fields: fields,
      buttons: [
        { label: "Fly to", onClick: function () { host.flyTo(entry); } },
        {
          label: state.mode === MODE_TRUE ? "Compressed scale" : "True scale",
          onClick: function () { host.setMode(state.mode === MODE_TRUE ? MODE_COMPRESSED : MODE_TRUE); },
        },
      ],
    };
  }

  function ringAround(ring, entry) {
    ring.at(entry.x, entry.y, entry.z, { px: 22 });
    return entry;
  }

  function layerFor(state) {
    return host.picker.addLayer({
      points: function () { return state.world; },
      reach: function () { return PICK_REACH_PX; },
      tooltip: function (entry) { return entry.name + ", " + entry.kind; },
      hover: function (entry) {
        if (!entry) {
          hoverRing.hide();
          state.hovered = null;
          return;
        }
        state.hovered = ringAround(hoverRing, entry);
      },
      select: function (entry) {
        if (state.selected === entry) {
          deselect();
          return;
        }
        host.showInfo(infoSpecFor(entry, state));
        state.selected = ringAround(selectionRing, entry);
        state.built.highlight(entry.ref);
        if (host.picked) host.picked(entry.ref);
      },
    });
  }

  // A visually hidden button per body, for the keyboard and screen readers
  // (a canvas has no focusable parts).
  function bodyList(state) {
    if (!host.viewport) return null;
    const list = document.createElement("ul");
    list.className = "starmap-sr-list starmap-sr-system sr-only";
    state.world.forEach(function (entry) {
      const item = document.createElement("li");
      const button = document.createElement("button");
      button.type = "button";
      button.textContent = entry.name + " (" + entry.kind + ")";
      button.addEventListener("click", function () { state.layer.select(entry); });
      item.appendChild(button);
      list.appendChild(item);
    });
    host.viewport.appendChild(list);
    return list;
  }

  function deselect() {
    if (!open || !open.selected) return false;
    selectionRing.hide();
    open.selected = null;
    open.built.highlight(null);
    if (host.deselected) host.deselected();
    if (host.picked) host.picked(null);
    return true;
  }

  function close() {
    token += 1;
    if (!open) return;
    host.picker.removeLayer(open.layer);
    if (open.list && open.list.parentNode) open.list.parentNode.removeChild(open.list);
    host.scene.remove(open.wrapper);
    open.built.dispose();
    selectionRing.hide();
    hoverRing.hide();
    open = null;
    if (host.closed) host.closed();
  }

  // Loads the system at `href` (its page, whose /scene this fetches) at
  // `center` in the map's world (parsecs) and puts it in the scene:
  // resolves with {center, radiusPc, fitPoints}, or null when it is
  // replaced meanwhile or fails to load.
  function openSystem(href, name, center, mode) {
    const mine = ++token;
    return host.fetchScene(href).then(function (data) {
      if (mine !== token || !data) return null;
      if (open) {
        const keep = token;
        close();
        token = keep;
      }
      const state = {
        href: href, id: data.system.id, name: name || data.system.name, data: data, mode: mode || MODE_COMPRESSED,
        center: center.slice(), built: null, wrapper: new THREE.Group(), world: [], selected: null, hovered: null,
        clock: createOrbitClock({ epochUnix: data.epoch_unix, nowMs: Date.now() }),
      };
      state.wrapper.position.set(center[0], center[1], center[2]);
      host.scene.add(state.wrapper);
      open = state;
      rebuild(state);
      state.layer = layerFor(state);
      state.list = bodyList(state);
      if (host.opened) host.opened();
      return { center: state.center, radiusPc: state.radiusPc, fitPoints: cubeAround(state.center, state.radiusPc) };
    });
  }

  // Changes the scale mode; resolves with the new fit like openSystem.
  function setMode(mode) {
    if (!open || open.mode === mode) return null;
    open.mode = mode;
    const keep = open.selected ? open.selected.ref : null;
    selectionRing.hide();
    open.selected = null;
    hoverRing.hide();
    open.hovered = null;
    if (open.list && open.list.parentNode) open.list.parentNode.removeChild(open.list);
    rebuild(open);
    open.list = bodyList(open);
    if (keep) {
      const again = open.world.find(function (w) { return w.ref === keep; });
      if (again) {
        open.selected = ringAround(selectionRing, again);
        open.built.highlight(again.ref);
      }
    }
    return { center: open.center, radiusPc: open.radiusPc, fitPoints: cubeAround(open.center, open.radiusPc) };
  }

  // Every frame: the bodies move on with the clock, the rings follow.
  function update(nowMs) {
    if (!open) return;
    open.clock.tick(nowMs);
    open.built.update(open.clock.years());
    syncWorld(open);
    [[selectionRing, open.selected], [hoverRing, open.hovered]].forEach(function (pair) {
      if (pair[1]) pair[0].at(pair[1].x, pair[1].y, pair[1].z, { px: 22 });
    });
    selectionRing.update();
    hoverRing.update();
  }

  return {
    open: openSystem,
    close: close,
    update: update,
    setMode: setMode,
    deselect: deselect,
    isOpen: function () { return !!open; },
    systemId: function () { return open ? open.id : null; },
    name: function () { return open ? open.name : null; },
    mode: function () { return open ? open.mode : MODE_COMPRESSED; },
    // The open system's body with this ref ("planet:12"), or null.
    entry: function (ref) { return open ? open.world.find(function (w) { return w.ref === ref; }) || null : null; },
    entries: function () { return open ? open.world : []; },
    select: function (ref) {
      const entry = open && open.world.find(function (w) { return w.ref === ref; });
      if (entry) open.layer.select(entry);
    },
    selected: function () { return open && open.selected ? open.selected.ref : null; },
    // The opacity of a body's orbit line, for tests (MAP.126).
    trailOpacity: function (ref) { return open ? open.built.trailOpacity(ref) : null; },
    // How far a body's orbit line reaches from what it goes round, for tests (MAP.136).
    trailRadius: function (ref) { return open ? open.built.trailRadius(ref) : null; },
    // Where a body is in the map's world now: {center, radius}, for flying to it.
    where: function (ref) {
      const entry = open && open.world.find(function (w) { return w.ref === ref; });
      return entry ? { center: [entry.x, entry.y, entry.z], radius: entry.r } : null;
    },
  };
}
