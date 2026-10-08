// html/static/systemview3d.js
//
// The 3D system view (MAP.72, MAP.73): one star system from
// GET /api/systems/<id>/scene drawn with three.js on the shared map
// controller (mapcontrol.js), picker and info panel (mappick.js).
//
// buildSystemScene() makes the THREE.Group: stars as spheres with the glow
// shell of bodyRendering.js and a light each, planets, moons and comets as
// lit spheres, every orbit as a line (a 3D ellipse for a comet), belts as
// particle rings, and the heliopause as a faint sphere. Each body is placed
// at its own position in doubles and every orbit line and belt hangs from
// its parent's position with small local coordinates, so the camera-relative
// maths three.js does keeps a moon steady next to a planet 40 AU out
// (floating origin). Where bodies are, and how large, is systemscale.js's
// layout (true or compressed); when they are is orbitpositions.js and the
// clock.
//
// createSystemView() adds the camera: drag to turn, shift-drag or right-drag
// to pan, wheel and pinch to zoom, W/A/S/D and the arrow keys to fly, a click
// to select, a double click to fly to a body, "follow" to ride along with
// one, and reset. Labels are plain elements laid over the canvas, the
// nearest bodies first and none over another.

const VERSION_QUERY = new URL(import.meta.url).search;
const THREE = await import(`./vendor/three.module.min.js${VERSION_QUERY}`);
const { makeGlowMaterial } = await import(`./bodyRendering.js${VERSION_QUERY}`);
const { createPointerControl, orbitByDrag, orbitByKey, panInScreenPlane, wheelPixels } = await import(`./mapcontrol.js${VERSION_QUERY}`);
const { fitRendererToCanvas, watchResize } = await import(`./mapcore.js${VERSION_QUERY}`);
const { createPicker, createTooltip, infoPanelOf } = await import(`./mappick.js${VERSION_QUERY}`);
const { orbitPath, relativeAt } = await import(`./orbitpositions.js${VERSION_QUERY}`);
const { createLayout, layoutPositions, MODE_COMPRESSED } = await import(`./systemscale.js${VERSION_QUERY}`);

const AU_KM = 149597870.7;
const SPHERE_SEGMENTS = [28, 20];
const BELT_PARTICLES = 700;
const PICK_REACH_PX = 16;
const MIN_PHI = 0.05;
const MAX_PHI = Math.PI - 0.05;

// A small repeatable random stream (mulberry32), so a belt looks the same
// every time it is drawn.
function randomStream(seed) {
  let a = (seed >>> 0) || 1;
  return function () {
    a = (a + 0x6d2b79f5) >>> 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

function colorOf(body) {
  return new THREE.Color(body.color || "#ffffff");
}

function lineFrom(layout, around, rel, color, opacity) {
  const points = rel.map((p) => new THREE.Vector3(...layout.place(around, p)));
  const geometry = new THREE.BufferGeometry().setFromPoints(points);
  const material = new THREE.LineBasicMaterial({ color: color, transparent: true, opacity: opacity, depthWrite: false });
  return new THREE.Line(geometry, material);
}

// Every body of the scene as {ref, kind, name, body, parent (ref or
// "barycenter")}, parents before children.
function listBodies(scene) {
  const out = [];
  scene.stars.forEach((s) => out.push({ ref: s.ref, kind: "star", name: s.name, body: s, parent: s.orbit ? s.orbit.around : "barycenter" }));
  scene.planets.forEach((p) => {
    out.push({ ref: p.ref, kind: "planet", name: p.name, body: p, parent: p.orbit.around });
    p.moons.forEach((m) => out.push({ ref: m.ref, kind: "moon", name: m.name, body: m, parent: p.ref }));
  });
  scene.comets.forEach((c) => out.push({ ref: c.ref, kind: "comet", name: c.name, body: c, parent: c.orbit.around }));
  return out;
}

// The scene's THREE.Group for `layout`. Returns {group, entries, update(years
// relative), dispose()}; `entries` are the pickable bodies with world x, y, z
// and radius kept current by update().
export function buildSystemScene(scene, layout, options) {
  const opts = options || {};
  const group = new THREE.Group();
  const bodies = listBodies(scene);
  const entries = [];
  const byRef = {};
  const objects = {};
  const hangers = {}; // parent ref -> orbit lines and belts that follow it
  const lights = [];
  const trails = {}; // ref -> its orbit line

  function hang(parent, object) {
    group.add(object);
    (hangers[parent] = hangers[parent] || []).push(object);
  }

  group.add(new THREE.AmbientLight(0xffffff, 0.35));

  bodies.forEach((item) => {
    const body = item.body;
    const radius = Math.max(layout.radiusOf(item.ref), 1e-9);
    const color = colorOf(body);
    const mesh = new THREE.Mesh(
      new THREE.SphereGeometry(1, SPHERE_SEGMENTS[0], SPHERE_SEGMENTS[1]),
      item.kind === "star"
        ? new THREE.MeshBasicMaterial({ color: color })
        : new THREE.MeshLambertMaterial({ color: color, emissive: color, emissiveIntensity: 0.12 }));
    mesh.scale.setScalar(radius);
    mesh.userData.ref = item.ref;
    group.add(mesh);
    objects[item.ref] = mesh;
    if (item.kind === "star") {
      // The Galaxy Map draws on a logarithmic depth buffer, which the glow
      // shader doesn't write, so a system drawn there has none (opts.glow false).
      if (opts.glow !== false) {
        const glow = new THREE.Mesh(new THREE.SphereGeometry(1, 20, 14), makeGlowMaterial(THREE, color, 2.2, 0.9, 2.2));
        mesh.add(glow);
        glow.scale.setScalar(2.2);
      }
      const light = new THREE.PointLight(color, 2.2, 0, 0);
      group.add(light);
      lights.push({ ref: item.ref, light: light });
    }
    const entry = {
      ref: item.ref, key: item.ref, name: item.name, kind: item.kind, body: body, x: 0, y: 0, z: 0, radius: radius,
    };
    entries.push(entry);
    byRef[item.ref] = entry;

    if (item.kind === "planet" || item.kind === "moon" || item.kind === "comet"
        || (item.kind === "star" && body.orbit && body.orbit.around !== "barycenter")) {
      const path = orbitPath(body.orbit);
      const around = body.orbit.around;
      const line = lineFrom(layout, around, path, color, item.kind === "moon" ? 0.45 : 0.6);
      line.userData.ref = item.ref;
      line.userData.base = { color: line.material.color.clone(), opacity: line.material.opacity };
      trails[item.ref] = line;
      hang(around, line);
    }
  });

  // A dot for every body, a few pixels across whatever its size, so a body
  // too small to see at true scale is still found; the sphere covers it once
  // it is bigger than that.
  const dotPositions = new Float32Array(entries.length * 3);
  const dotColors = new Float32Array(entries.length * 3);
  entries.forEach((entry, n) => {
    const c = colorOf(entry.body);
    dotColors.set([c.r, c.g, c.b], n * 3);
  });
  const dotGeometry = new THREE.BufferGeometry();
  dotGeometry.setAttribute("position", new THREE.BufferAttribute(dotPositions, 3));
  dotGeometry.setAttribute("color", new THREE.BufferAttribute(dotColors, 3));
  const dots = new THREE.Points(dotGeometry, new THREE.PointsMaterial({
    size: 5, sizeAttenuation: false, vertexColors: true, transparent: true, opacity: 0.95, depthWrite: false,
  }));
  dots.frustumCulled = false;
  group.add(dots);

  scene.belts.forEach((belt) => {
    const random = randomStream(belt.id * 2654435761);
    const positions = new Float32Array(BELT_PARTICLES * 3);
    for (let n = 0; n < BELT_PARTICLES; n += 1) {
      const km = belt.inner_km + (belt.outer_km - belt.inner_km) * random();
      const angle = random() * Math.PI * 2;
      const lift = (random() - 0.5) * 0.08 * km;
      const placed = layout.place(belt.around, [km * Math.cos(angle), km * Math.sin(angle), lift]);
      positions.set(placed, n * 3);
    }
    const geometry = new THREE.BufferGeometry();
    geometry.setAttribute("position", new THREE.BufferAttribute(positions, 3));
    const points = new THREE.Points(geometry, new THREE.PointsMaterial({
      color: 0xb8a58a, size: 1.6, sizeAttenuation: false, transparent: true, opacity: 0.8, depthWrite: false,
    }));
    points.userData.ref = belt.ref;
    hang(belt.around, points);
  });

  if (scene.system.heliopause_au) {
    const radius = layout.distance("barycenter", scene.system.heliopause_au * AU_KM);
    const sphere = new THREE.Mesh(
      new THREE.SphereGeometry(radius, 36, 24),
      new THREE.MeshBasicMaterial({ color: 0x6f8fff, transparent: true, opacity: 0.05, side: THREE.BackSide, depthWrite: false }));
    sphere.userData.ref = scene.system.ref;
    group.add(sphere);
  }

  // Moves everything to `years` after the epoch (or to `relative`, the
  // answer of relativeAt, when given).
  function update(years, relative) {
    const at = relative || relativeAt(scene, years);
    const placed = layoutPositions(layout, at);
    for (const ref of Object.keys(placed)) {
      const p = placed[ref];
      const entry = byRef[ref];
      objects[ref].position.set(p[0], p[1], p[2]);
      entry.x = p[0];
      entry.y = p[1];
      entry.z = p[2];
    }
    entries.forEach((entry, n) => {
      dotPositions[n * 3] = entry.x;
      dotPositions[n * 3 + 1] = entry.y;
      dotPositions[n * 3 + 2] = entry.z;
    });
    dotGeometry.attributes.position.needsUpdate = true;
    lights.forEach((l) => l.light.position.copy(objects[l.ref].position));
    for (const parent of Object.keys(hangers)) {
      const p = parent === "barycenter" ? [0, 0, 0] : placed[parent];
      hangers[parent].forEach((object) => object.position.set(p[0], p[1], p[2]));
    }
    return placed;
  }

  // MAP.126: the selected body's path is drawn bright and full, in the frame
  // of what it goes round (the lines hang from their parent), and the others
  // fade back; null puts them all as they were.
  function highlight(ref) {
    Object.keys(trails).forEach((key) => {
      const line = trails[key];
      const base = line.userData.base;
      if (!ref) {
        line.material.color.copy(base.color);
        line.material.opacity = base.opacity;
      } else if (key === ref) {
        line.material.color.copy(base.color).lerp(new THREE.Color("#ffffff"), 0.55);
        line.material.opacity = 1;
      } else {
        line.material.color.copy(base.color);
        line.material.opacity = base.opacity * 0.3;
      }
    });
  }

  function dispose() {
    group.traverse((object) => {
      if (object.geometry) object.geometry.dispose();
      if (object.material) object.material.dispose();
    });
  }

  update(opts.years || 0);
  return { group: group, entries: entries, byRef: byRef, objects: objects, update: update, highlight: highlight,
    // A body's orbit line opacity (null: it has none), for tests.
    trailOpacity: (ref) => (trails[ref] ? trails[ref].material.opacity : null), dispose: dispose };
}

// --- The view -------------------------------------------------------------------------

function easeInOut(t) {
  return t * t * (3 - 2 * t);
}

// options: canvasEl, labelsEl (an element over the canvas), tooltipEl,
// infoEl (the info panel), scene (the /scene JSON), clock (orbitclock.js),
// mode (systemscale.js), onSelect(entry or null), onMode(mode).
export function createSystemView(options) {
  const canvasEl = options.canvasEl;
  const sceneData = options.scene;
  const renderer = new THREE.WebGLRenderer({ canvas: canvasEl, antialias: true, alpha: true });
  renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
  const stage = new THREE.Scene();
  const camera = new THREE.PerspectiveCamera(50, 1, 1e-6, 1e7);
  const target = new THREE.Vector3();
  const state = { theta: 0.6, phi: 1.15, distance: 1, mode: options.mode || MODE_COMPRESSED };
  const follow = { ref: null };
  let flight = null;
  let built = null;
  let layout = null;
  let selected = null;
  let hovered = null;
  let fitDistance = 1;
  let years = 0;
  let disposed = false;
  const picker = createPicker(camera, canvasEl);
  const tooltip = createTooltip(options.tooltipEl);
  const info = infoPanelOf(options.infoEl);
  const labels = new Map();

  function sceneRadius() {
    let far = 1;
    built.entries.forEach((e) => { far = Math.max(far, Math.hypot(e.x, e.y, e.z)); });
    return far;
  }

  function placeCamera() {
    const sinPhi = Math.sin(state.phi);
    camera.position.set(
      target.x + state.distance * sinPhi * Math.cos(state.theta),
      target.y + state.distance * Math.cos(state.phi),
      target.z + state.distance * sinPhi * Math.sin(state.theta));
    camera.near = Math.max(state.distance * 1e-4, 1e-9);
    camera.far = Math.max(state.distance * 1e5, 10);
    camera.updateProjectionMatrix();
    camera.lookAt(target);
  }

  function rebuild() {
    if (built) {
      stage.remove(built.group);
      built.dispose();
    }
    layout = createLayout(sceneData, state.mode);
    built = buildSystemScene(sceneData, layout, { years: years });
    stage.add(built.group);
    labels.forEach((el) => el.remove());
    labels.clear();
    built.entries.forEach((entry) => {
      const el = document.createElement("span");
      el.className = "sysview-label sysview-label-" + entry.kind;
      el.textContent = entry.name;
      options.labelsEl.appendChild(el);
      labels.set(entry.ref, el);
    });
    if (selected) {
      selected = built.byRef[selected.ref] || null;
      built.highlight(selected ? selected.ref : null);
    }
  }

  function fit() {
    fitDistance = Math.max(sceneRadius() * 2.4, 1);
    return fitDistance;
  }

  function reset() {
    follow.ref = null;
    target.set(0, 0, 0);
    state.distance = fit();
    state.theta = 0.6;
    state.phi = 1.15;
    placeCamera();
  }

  function flyTo(ref) {
    const entry = built.byRef[ref];
    if (!entry) return;
    flight = {
      start: performance.now(), duration: 700, ref: ref,
      from: { target: target.clone(), distance: state.distance },
      distance: Math.max(entry.radius * 6, 1e-6),
    };
  }

  function stepFlight(now) {
    if (!flight) return;
    const entry = built.byRef[flight.ref];
    const t = Math.min(1, (now - flight.start) / flight.duration);
    const k = easeInOut(t);
    const goal = new THREE.Vector3(entry.x, entry.y, entry.z);
    target.lerpVectors(flight.from.target, goal, k);
    state.distance = Math.exp(Math.log(flight.from.distance) * (1 - k) + Math.log(flight.distance) * k);
    if (t >= 1) {
      follow.ref = flight.ref;
      flight = null;
    }
  }

  function infoSpec(entry) {
    const body = entry.body;
    const fields = [["Kind", entry.kind]];
    if (body.planet_class) fields.push(["Class", body.planet_class]);
    if (body.star_type) fields.push(["Star type", body.star_type]);
    if (body.radius_km) fields.push(["Radius", Math.round(body.radius_km).toLocaleString("en-US") + " km"]);
    if (body.orbit && body.orbit.distance_km) {
      fields.push(["Orbit", (body.orbit.distance_km / AU_KM).toPrecision(3) + " AU"]);
      if (body.orbit.period_years) fields.push(["Period", body.orbit.period_years.toPrecision(3) + " years"]);
    }
    return {
      title: entry.name,
      fields: fields,
      buttons: [
        { label: "Fly to", onClick: function () { flyTo(entry.ref); } },
        { label: "Follow", onClick: function () { follow.ref = entry.ref; flyTo(entry.ref); } },
      ],
    };
  }

  function select(entry) {
    selected = entry;
    if (built) built.highlight(entry ? entry.ref : null);
    if (entry) info.show(infoSpec(entry));
    if (options.onSelect) options.onSelect(entry);
  }

  picker.addLayer({
    points: () => (built ? built.entries : []),
    reach: () => PICK_REACH_PX,
  });

  const gesture = { distance0: 1 };
  const clampPhi = (phi) => Math.max(MIN_PHI, Math.min(MAX_PHI, phi));
  const control = createPointerControl(canvasEl, {
    attach: true,
    dragClickPx: 5,
    isPan: (event) => event.button === 2 || event.shiftKey,
    onDrag: function (dx, dy, pan) {
      follow.ref = pan ? null : follow.ref;
      if (pan) panInScreenPlane(THREE, camera, target, dx, dy, state.distance / Math.max(canvasEl.clientHeight, 1));
      else orbitByDrag(state, dx, dy, 0.006, clampPhi);
      flight = null;
    },
    pinch: {
      canStart: () => true,
      start: function () { gesture.distance0 = state.distance; },
      move: function (ratio) { state.distance = gesture.distance0 * ratio; },
    },
    onClick: function (event) {
      const hit = picker.pick(event.clientX, event.clientY);
      select(hit ? hit.entry : null);
    },
    onHover: function (event) {
      const hit = picker.pick(event.clientX, event.clientY);
      hovered = hit ? hit.entry : null;
      tooltip.show(hovered ? hovered.name : "", event.clientX, event.clientY);
    },
    onLeave: function () { hovered = null; tooltip.hide(); },
  });
  void control;
  canvasEl.addEventListener("contextmenu", (event) => event.preventDefault());
  canvasEl.addEventListener("dblclick", function (event) {
    const hit = picker.pick(event.clientX, event.clientY);
    if (hit) {
      select(hit.entry);
      flyTo(hit.entry.ref);
    }
  });
  canvasEl.addEventListener("wheel", function (event) {
    event.preventDefault();
    const pixels = wheelPixels(event, canvasEl.clientHeight, 200);
    state.distance *= Math.exp(pixels * 0.0015);
  }, { passive: false });
  canvasEl.tabIndex = 0;
  canvasEl.addEventListener("keydown", function (event) {
    const stepMove = state.distance * 0.04;
    if (orbitByKey(state, event.key, 0.08, clampPhi)) {
      event.preventDefault();
    } else if (event.key === "+" || event.key === "=") {
      state.distance *= 0.85;
      event.preventDefault();
    } else if (event.key === "-" || event.key === "_") {
      state.distance /= 0.85;
      event.preventDefault();
    } else if ("wasdWASD".indexOf(event.key) >= 0 && event.key.length === 1) {
      const forward = new THREE.Vector3().subVectors(target, camera.position).normalize();
      const right = new THREE.Vector3().crossVectors(forward, camera.up).normalize();
      const k = event.key.toLowerCase();
      if (k === "w") target.addScaledVector(forward, stepMove);
      if (k === "s") target.addScaledVector(forward, -stepMove);
      if (k === "d") target.addScaledVector(right, stepMove);
      if (k === "a") target.addScaledVector(right, -stepMove);
      follow.ref = null;
      event.preventDefault();
    } else if (event.key === "Home" || event.key === "0") {
      reset();
      event.preventDefault();
    }
  });

  function layoutLabels() {
    const rect = canvasEl.getBoundingClientRect();
    const taken = [];
    const order = built.entries.slice().sort((a, b) => {
      const da = camera.position.distanceToSquared(new THREE.Vector3(a.x, a.y, a.z));
      const db = camera.position.distanceToSquared(new THREE.Vector3(b.x, b.y, b.z));
      return da - db;
    });
    order.forEach((entry) => {
      const el = labels.get(entry.ref);
      const point = new THREE.Vector3(entry.x, entry.y, entry.z).project(camera);
      const x = ((point.x + 1) / 2) * rect.width;
      const y = ((1 - point.y) / 2) * rect.height;
      const behind = point.z > 1 || point.z < -1;
      const box = { x: x + 8, y: y - 8, w: entry.name.length * 7 + 6, h: 16 };
      const clash = taken.some((t) => box.x < t.x + t.w && t.x < box.x + box.w && box.y < t.y + t.h && t.y < box.y + box.h);
      const show = !behind && !clash && x > -50 && x < rect.width + 50 && y > -20 && y < rect.height + 20;
      el.hidden = !show;
      if (show) {
        el.style.transform = "translate(" + Math.round(box.x) + "px," + Math.round(box.y) + "px)";
        taken.push(box);
        el.classList.toggle("sysview-label-selected", entry === selected);
      }
    });
  }

  let frame = 0;
  function render(now) {
    if (disposed) return;
    frame = requestAnimationFrame(render);
    if (options.clock) years = options.clock.tick(Date.now());
    built.update(years);
    stepFlight(now);
    if (follow.ref && built.byRef[follow.ref] && !flight) {
      const e = built.byRef[follow.ref];
      target.set(e.x, e.y, e.z);
    }
    placeCamera();
    fitRendererToCanvas(renderer, camera, canvasEl);
    renderer.render(stage, camera);
    layoutLabels();
  }

  const stopWatching = watchResize(canvasEl, function () { fitRendererToCanvas(renderer, camera, canvasEl); });

  rebuild();
  reset();
  frame = requestAnimationFrame(render);

  return {
    setMode: function (mode) {
      state.mode = mode;
      rebuild();
      reset();
      if (options.onMode) options.onMode(layout);
    },
    mode: () => state.mode,
    note: () => layout.note,
    reset: reset,
    flyTo: flyTo,
    follow: function (ref) { follow.ref = ref; if (ref) flyTo(ref); },
    select: function (ref) { select(built.byRef[ref] || null); },
    selected: () => selected,
    entries: () => built.entries,
    // The bodies by name, for the screen-reader list.
    list: () => built.entries.map((e) => ({ ref: e.ref, name: e.name, kind: e.kind })),
    camera: camera,
    dispose: function () {
      disposed = true;
      cancelAnimationFrame(frame);
      if (stopWatching) stopWatching();
      built.dispose();
      renderer.dispose();
    },
  };
}
