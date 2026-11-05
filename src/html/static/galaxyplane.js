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
// The lines are built in coordinates local to the plane's anchor (the focus), and the group sits at the anchor, so the
// camera-relative arithmetic stays exact at the scales a few parsecs across (the MAP.240 class: a single-precision position
// 8 kpc out is good to a thousandth of a parsec). No DOM at import; three.js comes in as a parameter.

export const WEDGE_DEG = 60;
// Rings aim to be this many to a view's width, and no more than this many are drawn.
const RINGS_ACROSS = 4;
const MAX_RINGS = 40;
const ARC_STEP_DEG = 2;
const REBUILD_SHARE = 0.05;

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
  if (pc >= 1000) return `${Number((pc / 1000).toPrecision(3))} kpc`;
  if (pc >= 1) return `${Number(pc.toPrecision(3))} pc`;
  return `${Number(pc.toPrecision(3))} pc`;
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

  const labelBox = document.createElement("div");
  labelBox.className = "map-labels map-plane-labels";
  labelBox.setAttribute("aria-hidden", "true");
  if (parent) parent.appendChild(labelBox);

  const state = {
    on: false, midplane: false, centreAngle: -Math.PI / 2, selected: null,
    built: null, step: 0, ringCount: 0, labels: [], planeZ: 0, around: "galactic center",
  };
  const point = new THREE.Vector3();

  function setWedgeToward(target) {
    if (Math.hypot(target[0], target[1]) > 1e-6) state.centreAngle = Math.atan2(target[1], target[0]);
  }

  function planeHeight(target) {
    if (state.midplane) return 0;
    return state.selected ? state.selected[2] : target[2];
  }

  function clear() {
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
    const inWedge = (x, y) => Math.abs(wrap(Math.atan2(y, x) - state.centreAngle)) <= half + 1e-9;
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
          const angle = state.centreAngle - half + (2 * half * i) / count;
          const p = local(r * Math.cos(angle), r * Math.sin(angle));
          if (previous) ringVerts.push(...previous, ...p);
          previous = p;
        }
        pushLabel(distanceLabel(r), r * Math.cos(state.centreAngle), r * Math.sin(state.centreAngle), z);
        state.ringCount += 1;
      }
      const near = Math.max(0, radius - reach);
      const far = radius + reach;
      for (const sign of [-1, 1]) {
        const angle = state.centreAngle + sign * half;
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
        const angle = state.centreAngle + sign * half;
        edgeVerts.push(...local(near * Math.cos(angle), near * Math.sin(angle)), ...local(far * Math.cos(angle), far * Math.sin(angle)));
      }
    }
    rings.geometry.setAttribute("position", new THREE.Float32BufferAttribute(ringVerts, 3));
    edges.geometry.setAttribute("position", new THREE.Float32BufferAttribute(edgeVerts, 3));
    state.segments = ringVerts.length / 6 + edgeVerts.length / 6;
    state.built = { target: target.slice(), span: span, midplane: state.midplane, selected: state.selected ? state.selected.slice() : null, angle: state.centreAngle };
  }

  function stale(target) {
    const built = state.built;
    if (!built) return true;
    const span = Math.max(1e-9, spanNow());
    if (roundStep(span / RINGS_ACROSS) !== state.step) return true;
    if (built.midplane !== state.midplane || built.angle !== state.centreAngle) return true;
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
    state.labels.forEach(function (label) {
      point.set(label.world[0], label.world[1], label.world[2]).project(camera);
      const visible = point.z > -1 && point.z < 1 && Math.abs(point.x) <= 1 && Math.abs(point.y) <= 1;
      label.el.style.opacity = visible ? "0.9" : "0";
      if (visible) label.el.style.transform = `translate(${((point.x + 1) / 2) * width}px, ${((1 - point.y) / 2) * height}px)`;
    });
  }

  function update() {
    if (!state.on) return;
    const target = focus();
    if (!target) return;
    // The wedge stays where it was while the focus stays inside it, and swings round to the focus when it leaves.
    if (Math.abs(wrap(Math.atan2(target[1], target[0]) - state.centreAngle)) > half && Math.hypot(target[0], target[1]) > 1e-6) {
      setWedgeToward(target);
    }
    if (stale(target)) build(target);
    placeLabels();
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
      widest = Math.max(widest, Math.abs(wrap(Math.atan2(y, x) - state.centreAngle)));
    }
    return (widest * 180) / Math.PI;
  }

  canvasEl.galaxyPlane = function () {
    return {
      maxOffsetDeg: widestOffsetDeg(),
      on: state.on, midplane: state.midplane, planeZ: state.planeZ, step: state.step, rings: state.ringCount, segments: state.segments || 0,
      wedgeDeg: WEDGE_DEG, centreDeg: (state.centreAngle * 180) / Math.PI, around: state.around, anchor: group.position.toArray(),
      labels: state.labels.map(function (label) { return { text: label.text, world: label.world.slice(), shown: label.el.style.opacity === "0.9" }; }),
      visible: group.visible,
    };
  };

  return {
    set: set,
    isOn: function () { return state.on; },
    setMidplane: function (on) { state.midplane = !!on; },
    isMidplane: function () { return state.midplane; },
    select: function (at) { state.selected = at ? [at[0], at[1], at[2]] : null; },
    update: update,
  };
}
