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
// one, and reset. Names are the label engine's (maplabels.js, MAP.231): laid over the canvas, the nearest
// bodies first and none over another.

const VERSION_QUERY = new URL(import.meta.url).search;
const THREE = await import(`./vendor/three.module.min.js${VERSION_QUERY}`);
const { makeGlowMaterial, paintedTexture } = await import(`./bodyRendering.js${VERSION_QUERY}`);
const Spin = await import(`./spin.js${VERSION_QUERY}`);
const { makeCache, makeQueue, tierForDiameter, TIER_SIZES } = await import(`./bodytexture.js${VERSION_QUERY}`);
const { paintSolid } = await import(`./solidpaint.js${VERSION_QUERY}`);
const { paintGiant } = await import(`./giantpaint.js${VERSION_QUERY}`);
const { isLumpy, lumpyGeometry } = await import(`./lumpymesh.js${VERSION_QUERY}`);
const { paintStar } = await import(`./starpaint.js${VERSION_QUERY}`);
const { createPointerControl, orbitByDrag, orbitByKey, panInScreenPlane, wheelPixels } = await import(`./mapcontrol.js${VERSION_QUERY}`);
const { fitRendererToCanvas, watchResize } = await import(`./mapcore.js${VERSION_QUERY}`);
const { watchWindowFit } = await import(`./mapfit.js${VERSION_QUERY}`);
const { createPicker, createTooltip, infoPanelOf } = await import(`./mappick.js${VERSION_QUERY}`);
const { navActions, nearbyLinks, viewLink } = await import(`./systemnav.js${VERSION_QUERY}`);
const { fitDistance: framingDistance } = await import(`./flytofit.js${VERSION_QUERY}`);
const { orbitPath, nodePaths, relativeAt, cometOrbitFields } = await import(`./orbitpositions.js${VERSION_QUERY}`);
const { createLayout, layoutPositions, MODE_COMPRESSED } = await import(`./systemscale.js${VERSION_QUERY}`);
const { formatBodyRadiusKm, formatDistanceKm, formatPeriodYears } = await import(`./unitladder.js${VERSION_QUERY}`);
const { convert } = await import(`./unitregistry.js${VERSION_QUERY}`);
const { createLabelLayer } = await import(`./maplabels.js${VERSION_QUERY}`);

const AU_KM = convert(1, "au", "km");
const SPHERE_SEGMENTS = [28, 20];
// The size a body is taken to be on screen when the host cannot say (the Galaxy Map): the small texture.
const DEFAULT_DIAMETER_PX = 40;
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

// A body's colour: the look record's (MAP.287; a star, planet or moon), a comet's own `color`.
function colorOf(body) {
  return new THREE.Color((body.look && body.look.color) || body.color || "#ffffff");
}

// The shell of a planet's air, as the System Map draws it (systemmap.js).
// MAP.277: the glow of a body's air reaches `look.glowWidth` radii past the limb (this many when the look does not say: the composition of the air is not stored),
// as strongly as `look.glowStrength`; the night side is lit by the light the air scatters round the terminator, `look.twilight` (none for a world with no air).
const AIR_SHELL = 1.14;
const AIR_GLOW_POWER = 1.4;
const NIGHT_FILL_NO_LOOK = 0.12;      // a body with no look record (what the system view draws as a plain sphere) keeps the old fill
const NIGHT_FILL_NO_AIR = 0.04;
const NIGHT_FILL_PER_TWILIGHT = 0.494;
function nightFill(look) {
  return look ? NIGHT_FILL_NO_AIR + NIGHT_FILL_PER_TWILIGHT * (look.twilight || 0) : NIGHT_FILL_NO_LOOK;
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
  const beltAnchors = []; // a belt's symbol point (MAP.321): {ref, kind, name, body, x, y, z}, kept in step by update
  const lights = [];
  const trails = {}; // ref -> its orbit line
  // The painted surfaces (MAP.287): the System Map's painters (solidpaint.js, starpaint.js) one texture per body and tier,
  // painted one per frame, nearest first; until a body's texture is ready it is a sphere in its look colour.
  const textures = makeCache(64, (texture) => texture.dispose());
  const paintQueue = makeQueue();
  const paintable = []; // [{ ref, look }]
  const madeTextures = []; // every texture made here, for dispose

  function hang(parent, object) {
    group.add(object);
    (hangers[parent] = hangers[parent] || []).push(object);
  }

  group.add(new THREE.AmbientLight(0xffffff, 0.35));

  bodies.forEach((item) => {
    const body = item.body;
    const radius = Math.max(layout.radiusOf(item.ref), 1e-9);
    const color = colorOf(body);
    const lump = isLumpy(body.look) ? lumpyGeometry(THREE, body.look, SPHERE_SEGMENTS[0], SPHERE_SEGMENTS[1]) : null;     // MAP.279: a body its weight has not made round
    const mesh = new THREE.Mesh(
      lump ? lump.geometry : new THREE.SphereGeometry(1, SPHERE_SEGMENTS[0], SPHERE_SEGMENTS[1]),
      item.kind === "star"
        ? new THREE.MeshBasicMaterial({ color: color })
        : new THREE.MeshLambertMaterial({ color: color, emissive: color, emissiveIntensity: nightFill(body.look) }));
    mesh.scale.setScalar(radius);
    mesh.userData.ref = item.ref;
    group.add(mesh);
    objects[item.ref] = mesh;
    const look = body.look || null;
    if (look && (look.body === "solid" || look.body === "gas" || look.body === "star")) {
      paintable.push({ ref: item.ref, look: look });
    }
    if (look && look.glowStrength && look.body !== "star" && opts.glow !== false) {
      const shell = look.glowWidth ? 1 + look.glowWidth : AIR_SHELL;
      const air = new THREE.Mesh(
        new THREE.SphereGeometry(1, 20, 14), makeGlowMaterial(THREE, look.glow || "#bcdfff", AIR_GLOW_POWER, look.glowStrength, shell));
      air.scale.setScalar(shell);
      mesh.add(air);
    }
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
      ref: item.ref, key: item.ref, name: item.name, kind: item.kind, parent: item.parent, body: body, x: 0, y: 0, z: 0, radius: radius,
    };
    entries.push(entry);
    byRef[item.ref] = entry;

    if (item.kind === "planet" || item.kind === "moon" || item.kind === "comet") {
      const around = body.orbit.around;
      const path = orbitPath(body.orbit);
      const line = lineFrom(layout, around, path, color, item.kind === "moon" ? 0.45 : 0.6);
      line.userData.ref = item.ref;
      line.userData.base = { color: line.material.color.clone(), opacity: line.material.opacity };
      trails[item.ref] = line;
      hang(around, line);
    }
  });

  // MAP.213: the stars' paths are the orbits of the scene's tree, each child of a node on its own ellipse about the node's
  // centre of mass, the slot-0 child on the other side of it from the slot-1 child, each at the other's share of the
  // separation (MAP.136); a sub-pair's centre follows its own path about its parent node's centre.
  if (scene.stars.length) {
    const colors = {};
    scene.stars.forEach((star) => { colors[star.ref] = colorOf(star); });
    nodePaths(scene).forEach((path) => {
      const trailColor = colors[path.ref] || new THREE.Color("#ffffff");
      const line = lineFrom(layout, path.around, path.points, trailColor, 0.6);
      line.userData.ref = path.ref;
      line.userData.base = { color: line.material.color.clone(), opacity: line.material.opacity };
      trails[path.ref] = line;
      hang(path.around, line);
    });
  }

  // A dot for every body, a few pixels across whatever its size, so a body
  // too small to see at true scale is still found; the sphere covers it once
  // it is bigger than that.
  const dotPositions = new Float32Array(entries.length * 3);
  const dotColors = new Float32Array(entries.length * 3);
  entries.forEach((entry, n) => {
    const c = colorOf(entry.body);
    entry.dotColor = [c.r, c.g, c.b];
    dotColors.set(entry.dotColor, n * 3);
  });
  const dotGeometry = new THREE.BufferGeometry();
  dotGeometry.setAttribute("position", new THREE.BufferAttribute(dotPositions, 3));
  dotGeometry.setAttribute("color", new THREE.BufferAttribute(dotColors, 3));
  const dots = new THREE.Points(dotGeometry, new THREE.PointsMaterial({
    size: 5, sizeAttenuation: false, vertexColors: true, transparent: true, opacity: 0.95, depthWrite: false,
  }));
  dots.frustumCulled = false;
  dots.userData.relativeLayer = "system dots";   // MAP.240: in the system's own units, inside a group at the system
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
    points.userData.relativeLayer = "system belts";
    hang(belt.around, points);
    // MAP.321: where the belt's symbol sits, on the belt's ring half way across it, riding what the belt goes round.
    beltAnchors.push({
      ref: belt.ref, kind: "belt", name: belt.name || "Asteroid belt", body: belt, hung: points,
      local: layout.place(belt.around, [(belt.inner_km + belt.outer_km) / 2, 0, 0]), x: 0, y: 0, z: 0, radius: 0,
    });
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
      if (!entry) continue;    // a tree node's centre of mass: a point the stars hang from, not a body
      objects[ref].position.set(p[0], p[1], p[2]);
      entry.x = p[0];
      entry.y = p[1];
      entry.z = p[2];
    }
    // GEN.182: a comet that has left on an open orbit (or has not arrived) is off its path and is not shown.
    let shown = 0;
    entries.forEach((entry) => {
      entry.hidden = Boolean(at[entry.ref] && at[entry.ref].gone);
      objects[entry.ref].visible = !entry.hidden;
      if (entry.hidden) return;
      dotPositions[shown * 3] = entry.x;
      dotPositions[shown * 3 + 1] = entry.y;
      dotPositions[shown * 3 + 2] = entry.z;
      dotColors.set(entry.dotColor, shown * 3);
      shown += 1;
    });
    dotGeometry.setDrawRange(0, shown);
    dotGeometry.attributes.color.needsUpdate = true;
    dotGeometry.attributes.position.needsUpdate = true;
    lights.forEach((l) => l.light.position.copy(objects[l.ref].position));
    for (const parent of Object.keys(hangers)) {
      const p = parent === "barycenter" ? [0, 0, 0] : placed[parent];
      hangers[parent].forEach((object) => object.position.set(p[0], p[1], p[2]));
    }
    beltAnchors.forEach((anchor) => {
      anchor.x = anchor.hung.position.x + anchor.local[0];
      anchor.y = anchor.hung.position.y + anchor.local[1];
      anchor.z = anchor.hung.position.z + anchor.local[2];
    });
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

  // Paints the textures the bodies call for, one per call: `diameterPx(entry)` is a body's size on screen (the Galaxy
  // Map, which has no such size here, asks for none and every body gets the small tier).
  function paintTick(diameterPx) {
    paintQueue.runOne();
    paintable.forEach((one) => {
      const px = diameterPx ? diameterPx(byRef[one.ref]) : DEFAULT_DIAMETER_PX;
      const tier = tierForDiameter(px);
      const mesh = objects[one.ref];
      if (tier !== "low" && tier !== "high") return;
      const key = one.look.seed + ":" + tier;
      const texture = textures.get(key);
      if (!texture) {
        paintQueue.add(key, -px, () => {
          const made = paintedTexture(THREE, one.look.body === "star" ? paintStar : one.look.body === "gas" ? paintGiant : paintSolid, one.look, TIER_SIZES[tier]);
          madeTextures.push(made);
          textures.set(key, made);
        });
      } else if (mesh.userData.textureKey !== key) {
        // The painted star and the painted solid world already carry their own colour: a tint on top would double it.
        mesh.userData.textureKey = key;
        mesh.material.map = texture;
        mesh.material.color.set("#ffffff");
        mesh.material.needsUpdate = true;
      }
    });
  }

  function dispose() {
    paintQueue.clear();
    madeTextures.forEach((texture) => texture.dispose());
    group.traverse((object) => {
      if (object.geometry) object.geometry.dispose();
      if (object.material) object.material.dispose();
    });
  }

  update(opts.years || 0);

  // MAP.328: every body that has a stored rotation period turns about its stored axis, to where `spinSeconds` of spin time
  // (spin.js `tick`) has taken it; one without a period is left where it was. Called by the 3D system view each frame,
  // not by the Galaxy Map, which draws a system as a point of light.
  function spin(spinSeconds) {
    entries.forEach((entry) => {
      const look = entry.body.look;
      if (look && look.spinH > 0 && !entry.hidden) Spin.turn(THREE, objects[entry.ref], look, spinSeconds, null, Spin.isLocked(look) ? towardParent(entry) : null);
    });
  }

  // UX.244: the direction from a body to what it orbits (a barycentre: the nearest star), for a tidally locked body to keep its face toward.
  function towardParent(entry) {
    const own = objects[entry.ref].position;
    let target = objects[entry.parent];
    if (!target) {
      let best = Infinity;
      entries.forEach((other) => {
        const distance = other.kind === "star" ? own.distanceToSquared(objects[other.ref].position) : Infinity;
        if (distance < best) { best = distance; target = objects[other.ref]; }
      });
    }
    return target ? [target.position.x - own.x, target.position.y - own.y, target.position.z - own.z] : null;
  }

  function trailReach(ref) {
    if (!trails[ref]) return 0;
    const attribute = trails[ref].geometry.attributes.position;
    let farthest = 0;
    for (let n = 0; n < attribute.count; n += 1) {
      farthest = Math.max(farthest, Math.hypot(attribute.getX(n), attribute.getY(n), attribute.getZ(n)));
    }
    return farthest;
  }

  // MAP.319: the radius, in scene units, of the sphere about body `ref` that holds it and everything that goes round
  // it: the outermost orbit line of its moons (their apoapsis) and the moons themselves. Null for a body with no moon.
  function moonSystemReach(ref) {
    const center = byRef[ref];
    if (!center) return null;
    const moons = entries.filter((m) => m.body.orbit && m.body.orbit.around === ref && m.kind === "moon");
    if (!moons.length) return null;
    let far = center.radius;
    moons.forEach((m) => {
      far = Math.max(far, Math.hypot(m.x - center.x, m.y - center.y, m.z - center.z) + m.radius, trailReach(m.ref));
    });
    return far;
  }

  // MAP.319: the radius, in scene units, of the sphere about the scene's origin that holds every body and every orbit
  // line, each orbit measured from where the body it goes round is (the outermost stored body's apoapsis; for a wide
  // binary the second star's orbit about the first).
  function systemReach() {
    let far = 0;
    entries.forEach((e) => {
      far = Math.max(far, Math.hypot(e.x, e.y, e.z) + e.radius);
      const around = e.body.orbit && e.body.orbit.around;
      const parent = around && byRef[around];
      far = Math.max(far, (parent ? Math.hypot(parent.x, parent.y, parent.z) : 0) + trailReach(e.ref));
    });
    return far;
  }

  return { group: group, entries: entries, beltAnchors: beltAnchors, byRef: byRef, objects: objects, update: update, highlight: highlight, paintTick: paintTick, spin: spin,
    moonSystemReach: moonSystemReach, systemReach: systemReach,
    // The look-record texture a body shows now (null: still its flat colour), for tests.
    textureKey: (ref) => objects[ref].userData.textureKey || null,
    // A body's orbit line opacity (null: it has none), for tests.
    trailOpacity: (ref) => (trails[ref] ? trails[ref].material.opacity : null),
    // How far an orbit line reaches from what it goes round, in scene units (null: it has none), for tests.
    trailRadius: (ref) => (trails[ref] ? trailReach(ref) : null),
    dispose: dispose };
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
  // MAP.231: the names are the label engine's (maplabels.js), the same as the Galaxy Map's.
  const labelPoint = new THREE.Vector3();
  const labelLayer = createLabelLayer({
    container: options.labelsEl, budget: 400,
    viewport: () => ({ width: canvasEl.clientWidth || 1, height: canvasEl.clientHeight || 1 }),
    project: (x, y, z) => {
      labelPoint.set(x, y, z).applyMatrix4(camera.matrixWorldInverse);
      if (labelPoint.z >= 0) return null;
      const depth = -labelPoint.z;
      labelPoint.set(x, y, z).project(camera);
      const width = canvasEl.clientWidth || 1;
      const height = canvasEl.clientHeight || 1;
      return { x: (labelPoint.x + 1) / 2 * width, y: (1 - labelPoint.y) / 2 * height, depth };
    },
    reducedMotion: () => typeof matchMedia === "function" && matchMedia("(prefers-reduced-motion: reduce)").matches,
  });

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
    labelLayer.clear();
    if (selected) {
      selected = built.byRef[selected.ref] || null;
      built.highlight(selected ? selected.ref : null);
    }
  }

  // MAP.319: the camera distance at which a sphere of `radius` scene units fills the view.
  function framing(radius) {
    return framingDistance(radius, camera.fov, (canvasEl.clientWidth || 800) / (canvasEl.clientHeight || 600));
  }

  // The whole system fills the view: every body and every orbit line out to the outermost apoapsis.
  function fit() {
    fitDistance = Math.max(framing(Math.max(built.systemReach(), sceneRadius())), 1);
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

  // Moves the center the view turns and zooms about to the selected body, or
  // the middle of the system when nothing is selected, keeping the angle and
  // distance (MAP.138).
  function recenter() {
    follow.ref = null;
    flight = null;
    if (selected && built.byRef[selected.ref]) {
      const entry = built.byRef[selected.ref];
      target.set(entry.x, entry.y, entry.z);
    } else {
      target.set(0, 0, 0);
    }
    placeCamera();
  }

  function flyTo(ref) {
    const entry = built.byRef[ref];
    if (!entry) return;
    const moons = built.moonSystemReach(ref);
    flight = {
      start: performance.now(), duration: 700, ref: ref,
      from: { target: target.clone(), distance: state.distance },
      // A planet with moons fills the view with its moon system (MAP.319), any other body comes close.
      distance: Math.max(moons !== null ? framing(moons) : entry.radius * 6, 1e-6),
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
    if (body.radius_km) fields.push(["Radius", formatBodyRadiusKm(body.radius_km)]);
    if (body.orbit && body.orbit.distance_km) {
      fields.push(["Orbit", formatDistanceKm(body.orbit.distance_km)]);
      if (body.orbit.period_years) fields.push(["Period", formatPeriodYears(body.orbit.period_years)]);
    }
    if (entry.kind === "comet" && body.orbit && body.orbit.kepler) fields.push(...cometOrbitFields(body.orbit));
    return {
      title: entry.name,
      fields: fields,
      view: viewLink(document.getElementById("sysmap-root"), entry.ref, entry.kind),
      nav: navActions(document.getElementById("sysmap-root"), entry.ref),
      links: nearbyLinks(document.getElementById("sysmap-root"), entry.ref),
      buttons: [
        { label: "Go to", onClick: function () { flyTo(entry.ref); } },
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
    points: () => (built ? built.entries.filter((entry) => !entry.hidden) : []),
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
    } else if ((event.key === "c" || event.key === "C") && !event.ctrlKey && !event.metaKey && !event.altKey) {
      recenter();
      event.preventDefault();
    }
  });

  // The nearest bodies first (the selected one before any), each name placed by the engine so none lies over another.
  function layoutLabels() {
    const candidates = built.entries.map((entry) => ({
      id: entry.ref, kind: entry.kind, text: entry.name, x: entry.x, y: entry.y, z: entry.z,
      priority: (entry === selected ? 1e9 : 0) - camera.position.distanceTo(labelPoint.set(entry.x, entry.y, entry.z)),
      selected: entry === selected,
    }));
    labelLayer.update(candidates);
  }

  // A body's size on screen in pixels, for the tier of its texture.
  const eyeToBody = new THREE.Vector3();
  function diameterPx(entry) {
    const distance = Math.max(eyeToBody.set(entry.x, entry.y, entry.z).sub(camera.position).length(), 1e-12);
    return (2 * entry.radius / (distance * 2 * Math.tan(THREE.MathUtils.degToRad(camera.fov) / 2))) * (canvasEl.clientHeight || 1);
  }

  let frame = 0;
  function render(now) {
    if (disposed) return;
    frame = requestAnimationFrame(render);
    if (options.clock) years = options.clock.tick(Date.now());
    built.update(years);
    built.spin(Spin.tick(now));
    built.paintTick(diameterPx);
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
  // MAP.359: the 3D system view is as wide as its column and as tall as the Galaxy Map is for the browser window (95% of it).
  const viewportEl = canvasEl.closest(".sysview-viewport");
  const stopFitting = viewportEl ? watchWindowFit(viewportEl, null, { fill: true }) : null;

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
    recenter: recenter,
    flyTo: flyTo,
    // MAP.319, for the browser tests: how far the camera is, the sphere (scene units) a moon system or (no ref) the whole
    // system fills the view with, and the distance at which a sphere of that radius fills it.
    distance: () => state.distance,
    reachOf: (ref) => (ref ? built.moonSystemReach(ref) : built.systemReach()),
    framing: framing,
    flying: () => flight !== null,
    follow: function (ref) { follow.ref = ref; if (ref) flyTo(ref); },
    select: function (ref) { select(built.byRef[ref] || null); },
    selected: () => selected,
    entries: () => built.entries,
    // What each body is drawn from, for the browser tests: its look record's seed and body family, and the texture it
    // shows now (null: still its flat colour).
    bodies: () => built.entries.map((e) => ({ ref: e.ref, family: e.body.look ? e.body.look.body : null,
      seed: e.body.look ? e.body.look.seed : null, texture: built.textureKey(e.ref) })),
    // MAP.328, for the browser tests: how far each body has turned, its rate (radians a second) and its axis.
    spins: () => built.entries.map((e) => ({ ref: e.ref, angle: Spin.angleAt(e.body.look, Spin.seconds()),
      rate: Spin.rateOf(e.body.look), axis: Spin.axisOf(e.body.look), at: Spin.seconds(), locked: Spin.isLocked(e.body.look), quaternion: built.objects[e.ref].quaternion.toArray() })),
    // MAP.213, for the browser tests: where each body is drawn now, in world units, and how large.
    placed: () => built.entries.map((e) => ({ ref: e.ref, kind: e.kind, x: e.x, y: e.y, z: e.z, radius: e.radius })),
    // The bodies by name, for the screen-reader list.
    list: () => built.entries.map((e) => ({ ref: e.ref, name: e.name, kind: e.kind })),
    camera: camera,
    dispose: function () {
      disposed = true;
      cancelAnimationFrame(frame);
      if (stopWatching) stopWatching();
      if (stopFitting) stopFitting.stop();
      built.dispose();
      renderer.dispose();
    },
  };
}
