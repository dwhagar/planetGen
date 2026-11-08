// html/static/galaxyplane.js
//
// MAP.371: the plane view of the Galaxy Map (Boss, 5:36 PM and 5:48 PM EST on 2026-11-04). With it on, the camera takes a tilt
// preset (galaxystageview.js, `planeTilt`) and a 60 degree arc wedge of the galaxy, a pie slice from the galactic center, gets a
// grid on a plane: rings at round distances (1, 2 and 5 times a power of ten parsecs), the two edges of the wedge, and a label
// on each ring. The plane passes through the focus (the camera's target) and parallel to the galaxy's midplane, or, on the
// one-click choice, is the midplane itself (z = 0). With a star or other object selected the rings are drawn around it, and
// clipped to the wedge like the rest; with none they are arcs about the galactic center. The existing 3D mode is untouched:
// everything outside the wedge is drawn as before and the grid is a few thousand line vertices at most.
//
// MAP.372: stems. A thin line from each of a few objects down (or up) to its shadow on the grid plane, so height reads on the tilted
// map: the selection (the strongest), the objects the map has named (the label budget, MAP.317), the stars about the pointer and the
// brightest stars near the focus, never every star. Only the stems inside the wedge are drawn, one shorter than STEM_MIN_PX on the
// screen is skipped, and a far stem is fainter. They are rebuilt a few times a second, in the same anchor-local coordinates as the
// grid (the MAP.240 class of precision).
//
// MAP.373: the Menu holds a Stems switch and a height stretch (1 to 10 times: a stem is drawn on past the star to a bead at k times its
// height, and the grid says so in its corner note, so the exaggeration stays honest), and a stemmed or selected star carries a label
// of its spherical coordinates, an ordered triple, in the context: (radius, longitude, latitude) from the sector center in an opened
// sector or from the galactic center in the galaxy, and, with a star selected, (distance, bearing, elevation) of every other star from
// the selected one. The same triple is the infobox line of a star.
//
// The lines are built in coordinates local to the plane's anchor (the focus), and the group sits at the anchor, so the
// camera-relative arithmetic stays exact at the scales a few parsecs across (the MAP.240 class: a single-precision position
// 8 kpc out is good to a thousandth of a parsec). No DOM at import; three.js comes in as a parameter.

const VERSION_QUERY = new URL(import.meta.url).search;
const { threeFigures } = await import(`./numberformat.js${VERSION_QUERY}`);

export const WEDGE_DEG = 60;
// Rings aim to be this many to a view's width, and no more than this many are drawn.
const RINGS_ACROSS = 4;
const MAX_RINGS = 40;
const ARC_STEP_DEG = 2;
const REBUILD_SHARE = 0.05;
// MAP.372: the most stems drawn, the share of them the brightest stars near the focus may take, the pointer's reach in pixels, the
// radius about the focus (in views' widths) the bright stars are taken from, the stem a screen shorter than this many pixels is
// skipped, and the milliseconds between two choices of stems.
export const STEM_BUDGET = 60;
const BRIGHTEST_STEMS = 24;
const POINTER_REACH_PX = 90;
const STEM_RADIUS_SPANS = 0.75;
export const STEM_MIN_PX = 2;
const STEM_EVERY_MS = 150;
// MAP.373: the most coordinate labels at once (the selection and the stems the map has named or the pointer is near).
const COORD_LABELS = 10;
export const MAX_STRETCH = 10;

// The round distance (1, 2 or 5 times a power of ten) nearest `value` from above in log terms: the ring spacing.
export function roundStep(value) {
  if (!(value > 0)) return 1;
  const power = Math.pow(10, Math.floor(Math.log10(value)));
  const lead = value / power;
  const pick = lead < 1.5 ? 1 : lead < 3.5 ? 2 : lead < 7.5 ? 5 : 10;
  return pick * power;
}

// "500 pc", "1 kpc", "2.5 kpc", "50 AU" is not needed: the grid lives at parsec scale and above.
export function distanceLabel(pc) {
  if (pc >= 1000) return `${threeFigures(pc / 1000)} kpc`;
  return `${threeFigures(pc)} pc`;
}

// MAP.373: the spherical coordinates of `point` about `origin` in the galaxy frame: the distance (pc), the longitude (degrees from +x,
// 0 to 360, counterclockwise seen from +z) and the latitude (degrees above the xy-plane, -90 to 90).
export function sphericalAbout(point, origin) {
  const dx = point[0] - origin[0];
  const dy = point[1] - origin[1];
  const dz = point[2] - origin[2];
  const r = Math.hypot(dx, dy, dz);
  let lon = (Math.atan2(dy, dx) * 180) / Math.PI;
  if (lon < 0) lon += 360;
  if (lon >= 360) lon -= 360;
  const lat = r > 0 ? (Math.asin(Math.max(-1, Math.min(1, dz / r))) * 180) / Math.PI : 0;
  return { r: r, lonDeg: lon, latDeg: lat };
}

function signed(degrees) {
  return `${degrees < 0 ? "-" : "+"}${Math.abs(degrees).toFixed(1)}\u00b0`;
}

// "7.91 kpc, 82.4\u00b0, +0.4\u00b0": the triple in order (radius, longitude, latitude), with units.
export function tripleText(triple) {
  return `${distanceLabel(triple.r)}, ${triple.lonDeg.toFixed(1)}\u00b0, ${signed(triple.latDeg)}`;
}

function wrap(angle) {
  let a = angle;
  while (a > Math.PI) a -= 2 * Math.PI;
  while (a < -Math.PI) a += 2 * Math.PI;
  return a;
}

export function createPlaneView(options) {
  const THREE = options.THREE;
  const scene = options.scene;
  const camera = options.camera;
  const canvasEl = options.canvasEl;
  const parent = options.parent;
  const accent = options.accentColor || 0x7f8cff;
  const focus = options.focus;           // () => [x, y, z] the camera's target
  const spanNow = options.span;          // () => the width in parsecs the view covers at its focus
  const contextNow = options.context || function () { return { kind: "galaxy", origin: [0, 0, 0] }; };    // () => {kind, origin}
  const stemSources = options.stemSources || function () { return { named: [], stars: [] }; };    // () => {named, stars}, plain {x, y, z} lists
  const half = (WEDGE_DEG / 2) * Math.PI / 180;

  const group = new THREE.Group();
  group.visible = false;
  scene.add(group);
  const rings = new THREE.LineSegments(new THREE.BufferGeometry(), new THREE.LineBasicMaterial({ color: accent, transparent: true, opacity: 0.4, depthWrite: false }));
  const edges = new THREE.LineSegments(new THREE.BufferGeometry(), new THREE.LineBasicMaterial({ color: accent, transparent: true, opacity: 0.75, depthWrite: false }));
  rings.renderOrder = 2;
  edges.renderOrder = 2;
  rings.frustumCulled = false;
  edges.frustumCulled = false;
  group.add(rings);
  group.add(edges);
  const stems = new THREE.LineSegments(new THREE.BufferGeometry(), new THREE.LineBasicMaterial({ vertexColors: true, transparent: true, opacity: 0.8, depthWrite: false }));
  stems.renderOrder = 3;
  stems.frustumCulled = false;
  group.add(stems);
  const extensions = new THREE.LineSegments(new THREE.BufferGeometry(), new THREE.LineBasicMaterial({ vertexColors: true, transparent: true, opacity: 0.6, depthWrite: false }));
  extensions.renderOrder = 3;
  extensions.frustumCulled = false;
  group.add(extensions);
  const beads = new THREE.Points(new THREE.BufferGeometry(), new THREE.PointsMaterial({ size: 5, sizeAttenuation: false, vertexColors: true, transparent: true, depthWrite: false }));
  beads.renderOrder = 4;
  beads.frustumCulled = false;
  group.add(beads);
  const stemColor = new THREE.Color(accent);
  const stemSelectedColor = new THREE.Color(0xffffff);
  let pointerAt = null;
  if (canvasEl.addEventListener) {
    canvasEl.addEventListener("pointermove", function (event) { pointerAt = { x: event.clientX, y: event.clientY }; });
    canvasEl.addEventListener("pointerleave", function () { pointerAt = null; });
  }

  const labelBox = document.createElement("div");
  labelBox.className = "map-labels map-plane-labels";
  labelBox.setAttribute("aria-hidden", "true");
  if (parent) parent.appendChild(labelBox);

  const state = {
    on: false, midplane: false, centerAngle: -Math.PI / 2, selected: null,
    built: null, step: 0, ringCount: 0, labels: [], planeZ: 0, around: "galactic center",
    stemsOn: true, stretch: 1, coordLabels: [], labelPool: new Map(),
    stemCount: 0, stemShortestPx: null, stemMs: 0, stemRuns: 0, stemChosenAt: -1e9, stemKinds: { selected: 0, named: 0, pointer: 0, bright: 0 },
  };
  const point = new THREE.Vector3();

  function setWedgeToward(target) {
    if (Math.hypot(target[0], target[1]) > 1e-6) state.centerAngle = Math.atan2(target[1], target[0]);
  }

  function planeHeight(target) {
    if (state.midplane) return 0;
    return state.selected ? state.selected[2] : target[2];
  }

  function clearStems() {
    [stems, extensions, beads].forEach(function (layer) {
      layer.geometry.dispose();
      layer.geometry = new THREE.BufferGeometry();
    });
    state.coordLabels.forEach(function (label) { label.el.remove(); });
    state.coordLabels = [];
    state.labelPool = new Map();
    state.stemCount = 0;
    state.stemKinds = { selected: 0, named: 0, pointer: 0, bright: 0 };
  }

  // The height a point is drawn at with the stretch: k times its height above the plane (k = 1 is true).
  function stretchedZ(z, planeZ) {
    return planeZ + (z - planeZ) * state.stretch;
  }

  // The note on the grid that says what the heights are drawn at, so the exaggeration is never hidden.
  const note = document.createElement("div");
  note.className = "map-label map-plane-note";
  note.dataset.kind = "plane-note";
  function noteText() {
    return state.stretch > 1
      ? `Heights drawn ${state.stretch} times taller than true (the bead tops each stem)`
      : "Heights at true scale";
  }

  // The label of one object's coordinates in the context: absolute (radius, longitude, latitude) about the context's center, or, for
  // every object but the selected one when something is selected, (distance, bearing, elevation) from the selected one.
  function coordinateText(at, isSelected) {
    if (state.selected && !isSelected) {
      return `${tripleText(sphericalAbout(at, state.selected))} from the selection`;
    }
    const context = contextNow();
    return tripleText(sphericalAbout(at, context.origin));
  }

  // The screen length, in pixels, of the line from a point straight down to the plane.
  const upper = new THREE.Vector3();
  const lower = new THREE.Vector3();
  function stemPixels(x, y, z, planeZ) {
    const width = canvasEl.clientWidth || 1;
    const height = canvasEl.clientHeight || 1;
    upper.set(x, y, z).project(camera);
    lower.set(x, y, planeZ).project(camera);
    if (upper.z > 1 || lower.z > 1 || upper.z < -1 || lower.z < -1) return 0;
    return Math.hypot((upper.x - lower.x) * width / 2, (upper.y - lower.y) * height / 2);
  }

  // Chooses the stems and builds their lines: {selected, named, pointer, bright} counts are kept for the tests. Called at most every
  // STEM_EVERY_MS, and at once when the plane or the selection changes.
  function buildStems(target, now) {
    const began = performance.now();
    const span = Math.max(1e-9, spanNow());
    const planeZ = state.planeZ;
    const width = canvasEl.clientWidth || 1;
    const height = canvasEl.clientHeight || 1;
    const rect = canvasEl.getBoundingClientRect ? canvasEl.getBoundingClientRect() : { left: 0, top: 0 };
    const sources = stemSources();
    const chosen = [];
    const taken = new Set();
    const kinds = { selected: 0, named: 0, pointer: 0, bright: 0 };
    let shortest = Infinity;
    function offer(item, kind, strong) {
      if (!item || chosen.length >= STEM_BUDGET) return;
      const key = `${item.x},${item.y},${item.z}`;
      if (taken.has(key) || !insideWedge(item.x, item.y)) return;
      const pixels = stemPixels(item.x, item.y, stretchedZ(item.z, planeZ), planeZ);
      if (pixels < STEM_MIN_PX) return;
      shortest = Math.min(shortest, pixels);
      taken.add(key);
      kinds[kind] += 1;
      chosen.push({ x: item.x, y: item.y, z: item.z, strong: strong, kind: kind });
    }
    if (state.selected) offer({ x: state.selected[0], y: state.selected[1], z: state.selected[2] }, "selected", true);
    (sources.named || []).forEach(function (item) { offer(item, "named", false); });
    const stars = sources.stars || [];
    if (pointerAt) {
      const px = pointerAt.x - rect.left;
      const py = pointerAt.y - rect.top;
      const near = [];
      for (let i = stars.length - 1; i >= 0 && near.length < 4 * STEM_BUDGET; i -= 1) {
        point.set(stars[i].x, stars[i].y, stars[i].z).project(camera);
        if (point.z < -1 || point.z > 1) continue;
        const sx = ((point.x + 1) / 2) * width;
        const sy = ((1 - point.y) / 2) * height;
        if (Math.hypot(sx - px, sy - py) <= POINTER_REACH_PX) near.push(stars[i]);
      }
      near.forEach(function (star) { offer(star, "pointer", false); });
    }
    // The brightest are last in the list the map keeps (brightest last, so they win a tie on the screen).
    const reach = span * STEM_RADIUS_SPANS;
    let bright = 0;
    for (let i = stars.length - 1; i >= 0 && bright < BRIGHTEST_STEMS && chosen.length < STEM_BUDGET; i -= 1) {
      const star = stars[i];
      if (Math.hypot(star.x - target[0], star.y - target[1], star.z - target[2]) > reach) continue;
      const before = chosen.length;
      offer(star, "bright", false);
      if (chosen.length > before) bright += 1;
    }
    // With the stems switched off nothing is drawn but the selection's label.
    if (!state.stemsOn) {
      chosen.length = 0;
      kinds.selected = kinds.named = kinds.pointer = kinds.bright = 0;
    }
    const positions = [];
    const colors = [];
    const extra = [];
    const extraColors = [];
    const beadAt = [];
    const beadColors = [];
    const labelled = [];
    const anchor = state.built.anchor;
    const eye = camera.position;
    chosen.forEach(function (item) {
      // Fainter with distance from the eye, as far as three spans; the selection stays full strength.
      const away = Math.hypot(item.x - eye.x, item.y - eye.y, item.z - eye.z);
      const fade = item.strong ? 1 : Math.max(0.25, 1 - away / (3 * span));
      const color = item.strong ? stemSelectedColor : stemColor;
      const level = item.strong ? 1 : 0.55 * fade;
      positions.push(item.x - state.built.anchor[0], item.y - state.built.anchor[1], item.z - state.built.anchor[2],
        item.x - state.built.anchor[0], item.y - state.built.anchor[1], planeZ - state.built.anchor[2]);
      for (let k = 0; k < 2; k += 1) colors.push(color.r * level, color.g * level, color.b * level);
      // Stretched, the line runs on from the star to a bead at k times its height, in the same color.
      const top = stretchedZ(item.z, planeZ);
      if (state.stretch > 1) {
        extra.push(item.x - anchor[0], item.y - anchor[1], item.z - anchor[2], item.x - anchor[0], item.y - anchor[1], top - anchor[2]);
        for (let k = 0; k < 2; k += 1) extraColors.push(color.r * level, color.g * level, color.b * level);
        beadAt.push(item.x - anchor[0], item.y - anchor[1], top - anchor[2]);
        beadColors.push(color.r * level, color.g * level, color.b * level);
      }
      if (item.kind !== "bright" && labelled.length < COORD_LABELS) labelled.push({ item: item, top: top });
    });
    // The selection is labelled whether or not it has a stem (the plane may pass through it, or the stems be off).
    if (state.selected && !labelled.some(function (entry) { return entry.item.strong; })) {
      const at = state.selected;
      labelled.unshift({ item: { x: at[0], y: at[1], z: at[2], strong: true, kind: "selected" }, top: stretchedZ(at[2], planeZ) });
    }
    function fill(layer, verts, cols) {
      layer.geometry.dispose();
      layer.geometry = new THREE.BufferGeometry();
      layer.geometry.setAttribute("position", new THREE.Float32BufferAttribute(verts, 3));
      layer.geometry.setAttribute("color", new THREE.Float32BufferAttribute(cols, 3));
    }
    fill(stems, positions, colors);
    fill(extensions, extra, extraColors);
    fill(beads, beadAt, beadColors);
    // The coordinate labels, rebuilt with the stems: the selection first, then the names and the pointer's stars.
    // An element is kept for as long as its star is labelled (keyed by its place), so a rebuild never removes and re-adds a label (which
    // blinked it for a frame); only the stars that drop out lose theirs.
    const kept = new Map();
    state.coordLabels = labelled.map(function (entry) {
      const at = [entry.item.x, entry.item.y, entry.item.z];
      const key = `${at[0]},${at[1]},${at[2]}`;
      const old = state.labelPool.get(key);
      const el = old ? old.el : document.createElement("div");
      if (!old) {
        el.className = "map-label map-plane-coords";
        el.dataset.kind = "plane-coords";
        el.style.opacity = "0";
        labelBox.appendChild(el);
      }
      if (entry.item.strong) el.dataset.selected = "";
      else delete el.dataset.selected;
      const text = coordinateText(at, entry.item.strong);
      if (el.textContent !== text) el.textContent = text;
      const label = { el: el, text: text, world: [entry.item.x, entry.item.y, entry.top], star: at, selected: !!entry.item.strong };
      kept.set(key, label);
      return label;
    });
    state.labelPool.forEach(function (label, key) { if (!kept.has(key)) label.el.remove(); });
    state.labelPool = kept;
    note.textContent = noteText();
    state.stemCount = chosen.length;
    state.stemKinds = kinds;
    state.stemShortestPx = shortest;
    state.stemChosenAt = now;
    state.stemRuns += 1;
    state.stemMs += (performance.now() - began - state.stemMs) / Math.min(state.stemRuns, 20);
  }

  function insideWedge(x, y) {
    return Math.abs(wrap(Math.atan2(y, x) - state.centerAngle)) <= half + 1e-9;
  }

  function clear() {
    clearStems();
    rings.geometry.dispose();
    edges.geometry.dispose();
    rings.geometry = new THREE.BufferGeometry();
    edges.geometry = new THREE.BufferGeometry();
    labelBox.textContent = "";
    state.labels = [];
    state.ringCount = 0;
  }

  function pushLabel(text, x, y, z) {
    const el = document.createElement("div");
    el.className = "map-label";
    el.dataset.kind = "plane";
    el.textContent = text;
    labelBox.appendChild(el);
    state.labels.push({ text: text, world: [x, y, z], el: el });
  }

  function build(target) {
    clear();
    const span = Math.max(1e-9, spanNow());
    const step = roundStep(span / RINGS_ACROSS);
    const z = planeHeight(target);
    const anchor = [target[0], target[1], z];
    group.position.set(anchor[0], anchor[1], anchor[2]);
    state.step = step;
    state.planeZ = z;
    const ringVerts = [];
    const edgeVerts = [];
    const local = (x, y) => [x - anchor[0], y - anchor[1], 0];
    const inWedge = (x, y) => Math.abs(wrap(Math.atan2(y, x) - state.centerAngle)) <= half + 1e-9;
    const reach = span * 1.5;
    const around = state.selected ? [state.selected[0], state.selected[1]] : [0, 0];
    state.around = state.selected ? "the selection" : "galactic center";
    const arcStep = (ARC_STEP_DEG * Math.PI) / 180;

    if (!state.selected) {
      // Arcs about the galactic center across the wedge, from just inside the view to just outside.
      const radius = Math.hypot(target[0], target[1]);
      const first = Math.max(1, Math.ceil((radius - reach) / step));
      const last = Math.floor((radius + reach) / step);
      let drawn = 0;
      for (let k = first; k <= last && drawn < MAX_RINGS; k += 1, drawn += 1) {
        const r = k * step;
        let previous = null;
        const count = Math.max(2, Math.ceil((2 * half) / arcStep));
        for (let i = 0; i <= count; i += 1) {
          const angle = state.centerAngle - half + (2 * half * i) / count;
          const p = local(r * Math.cos(angle), r * Math.sin(angle));
          if (previous) ringVerts.push(...previous, ...p);
          previous = p;
        }
        pushLabel(distanceLabel(r), r * Math.cos(state.centerAngle), r * Math.sin(state.centerAngle), z);
        state.ringCount += 1;
      }
      const near = Math.max(0, radius - reach);
      const far = radius + reach;
      for (const sign of [-1, 1]) {
        const angle = state.centerAngle + sign * half;
        edgeVerts.push(...local(near * Math.cos(angle), near * Math.sin(angle)), ...local(far * Math.cos(angle), far * Math.sin(angle)));
      }
    } else {
      // Rings about the selection, drawn only where they lie inside the wedge.
      const count = Math.ceil((2 * Math.PI) / arcStep);
      const most = Math.min(MAX_RINGS, Math.max(1, Math.ceil(reach / step)));
      for (let k = 1; k <= most; k += 1) {
        const r = k * step;
        let labelled = false;
        let previous = null;
        for (let i = 0; i <= count; i += 1) {
          const angle = (2 * Math.PI * i) / count;
          const x = around[0] + r * Math.cos(angle);
          const y = around[1] + r * Math.sin(angle);
          const ok = inWedge(x, y);
          const p = local(x, y);
          if (previous && previous.ok && ok) ringVerts.push(...previous.p, ...p);
          if (ok && !labelled && previous && previous.ok && i > count / 8) {
            pushLabel(distanceLabel(r), x, y, z);
            labelled = true;
          }
          previous = { p: p, ok: ok };
        }
        state.ringCount += 1;
      }
      const far = Math.hypot(around[0], around[1]) + reach;
      const near = Math.max(0, Math.hypot(around[0], around[1]) - reach);
      for (const sign of [-1, 1]) {
        const angle = state.centerAngle + sign * half;
        edgeVerts.push(...local(near * Math.cos(angle), near * Math.sin(angle)), ...local(far * Math.cos(angle), far * Math.sin(angle)));
      }
    }
    rings.geometry.setAttribute("position", new THREE.Float32BufferAttribute(ringVerts, 3));
    edges.geometry.setAttribute("position", new THREE.Float32BufferAttribute(edgeVerts, 3));
    state.segments = ringVerts.length / 6 + edgeVerts.length / 6;
    state.built = { target: target.slice(), span: span, midplane: state.midplane, selected: state.selected ? state.selected.slice() : null, angle: state.centerAngle, anchor: anchor };
    state.stemChosenAt = -1e9;
  }

  function stale(target) {
    const built = state.built;
    if (!built) return true;
    const span = Math.max(1e-9, spanNow());
    if (roundStep(span / RINGS_ACROSS) !== state.step) return true;
    if (built.midplane !== state.midplane || built.angle !== state.centerAngle) return true;
    const was = built.selected;
    const now = state.selected;
    if ((was === null) !== (now === null)) return true;
    if (was && (was[0] !== now[0] || was[1] !== now[1] || was[2] !== now[2])) return true;
    if (Math.hypot(target[0] - built.target[0], target[1] - built.target[1]) > REBUILD_SHARE * span) return true;
    return !state.midplane && !state.selected && Math.abs(target[2] - built.target[2]) > REBUILD_SHARE * span;
  }

  function placeLabels() {
    const width = canvasEl.clientWidth || 1;
    const height = canvasEl.clientHeight || 1;
    state.labels.concat(state.coordLabels).forEach(function (label) {
      point.set(label.world[0], label.world[1], label.world[2]).project(camera);
      const visible = point.z > -1 && point.z < 1 && Math.abs(point.x) <= 1 && Math.abs(point.y) <= 1;
      label.el.style.opacity = visible ? "0.9" : "0";
      if (visible) label.el.style.transform = `translate(${((point.x + 1) / 2) * width}px, ${((1 - point.y) / 2) * height}px)`;
    });
    note.style.transform = `translate(8px, ${height - 24}px)`;
    note.style.opacity = "0.9";
  }

  function update() {
    if (!state.on) return;
    const target = focus();
    if (!target) return;
    // The wedge stays where it was while the focus stays inside it, and swings round to the focus when it leaves.
    if (Math.abs(wrap(Math.atan2(target[1], target[0]) - state.centerAngle)) > half && Math.hypot(target[0], target[1]) > 1e-6) {
      setWedgeToward(target);
    }
    if (stale(target)) build(target);
    const now = performance.now();
    if (!note.isConnected) labelBox.appendChild(note);
    note.textContent = noteText();
    if (now - state.stemChosenAt >= STEM_EVERY_MS) buildStems(target, now);
    placeLabels();    // after the stems, so a label is placed in the frame it is made
  }

  function set(on) {
    state.on = !!on;
    group.visible = state.on;
    labelBox.hidden = !state.on;
    if (state.on) {
      const target = focus();
      if (target) setWedgeToward(target);
      state.built = null;
    } else clear();
  }

  // The widest angle, in degrees, a ring's vertex lies from the middle of the wedge, seen from the galactic center (read by the browser tests).
  function widestOffsetDeg() {
    const position = rings.geometry.getAttribute("position");
    if (!position) return 0;
    let widest = 0;
    for (let i = 0; i < position.count; i += 1) {
      const x = position.getX(i) + group.position.x;
      const y = position.getY(i) + group.position.y;
      widest = Math.max(widest, Math.abs(wrap(Math.atan2(y, x) - state.centerAngle)));
    }
    return (widest * 180) / Math.PI;
  }

  canvasEl.galaxyPlane = function () {
    return {
      maxOffsetDeg: widestOffsetDeg(),
      stemsOn: state.stemsOn, stretch: state.stretch, note: noteText(),
      coords: state.coordLabels.map(function (label) { return { text: label.text, star: label.star.slice(), selected: label.selected, shown: label.el.style.opacity === "0.9" }; }),
      context: contextNow().kind,
      beads: beads.geometry.getAttribute("position") ? beads.geometry.getAttribute("position").count : 0,
      stems: state.stemCount, stemShortestPx: state.stemShortestPx === Infinity ? null : state.stemShortestPx, stemKinds: Object.assign({}, state.stemKinds), stemMs: state.stemMs, stemBudget: STEM_BUDGET,
      stemPositions: stems.geometry.getAttribute("position") ? Array.from(stems.geometry.getAttribute("position").array) : [],
      on: state.on, midplane: state.midplane, planeZ: state.planeZ, step: state.step, rings: state.ringCount, segments: state.segments || 0,
      wedgeDeg: WEDGE_DEG, centerDeg: (state.centerAngle * 180) / Math.PI, around: state.around, anchor: group.position.toArray(),
      labels: state.labels.map(function (label) { return { text: label.text, world: label.world.slice(), shown: label.el.style.opacity === "0.9" }; }),
      visible: group.visible,
    };
  };

  return {
    set: set,
    isOn: function () { return state.on; },
    // MAP.376: the wedge chosen by its middle line (radians from +x); the grid is rebuilt about it.
    setWedgeAngle: function (angle) { state.centerAngle = wrap(angle); state.stemChosenAt = -1e9; },
    wedgeAngle: function () { return state.centerAngle; },
    setMidplane: function (on) { state.midplane = !!on; },
    isMidplane: function () { return state.midplane; },
    select: function (at) { state.selected = at ? [at[0], at[1], at[2]] : null; state.stemChosenAt = -1e9; },
    setStems: function (on) { state.stemsOn = !!on; state.stemChosenAt = -1e9; },
    setStretch: function (value) {
      state.stretch = Math.max(1, Math.min(MAX_STRETCH, Math.round(Number(value) || 1)));
      state.stemChosenAt = -1e9;
    },
    stretch: function () { return state.stretch; },
    // The selection's coordinates in the context as the infobox line has them, or null with the plane view off.
    describe: function (at) { return state.on ? coordinateText([at[0], at[1], at[2]], true) : null; },
    update: update,
  };
}
