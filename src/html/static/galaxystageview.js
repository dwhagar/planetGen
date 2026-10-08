// html/static/galaxystageview.js
//
// The Galaxy Map's drill-down (docs/design/galaxy-drilldown-navigation.md,
// sections 4, 5, 8.1 and 10-11), drawn and driven: the whole galaxy in
// 3D with an arc of it picked on the map (MAP.85), then a slab of the
// arc (picked with the buttons beside the map or on the map), then a
// segment of the slab (one block), then a slab and a segment inside that
// block, ... down to a sector (MAP.56). galaxystages.js has the rules;
// this file has the scene, the camera moves, the breadcrumb, the slab
// buttons and their leader lines, the tooltip, keys, touch, the stage
// URLs and the map's own Back and Forward.
// galaxymap3d.js creates it (createStageView), hands it the pointer and
// key events and calls step() every frame.
//
// host (from galaxymap3d.js):
// - THREE, scene, camera, canvasEl;
// - edgePc, shape (or null), galaxyRadius, reducedMotion, accentColor;
// - blockScene: galaxyblocks' scene (buildCells);
// - makeBlockMesh(part, translucent): a mesh with the map's block shader;
// - setCamera({target: [x, y, z], dist, quaternion}): moves the camera
//   to `dist` from `target` along the quaternion's back (+z) axis, turned
//   by it
//   (and the tiles, scale bar and wedge lines with it);
// - setWedgeClip(clip): only the wedge in view shown, {r0, r1, a0, a1,
//   z0, z1, cells (its blocks' bounds)} (null: the whole galaxy, with its
//   bearing labels);
// - fetchStage(query): a Promise of GET /galaxy/stage's JSON for a
//   container ("?at=m.ring.wedge.slab", or "" for the galaxy);
// - showBlockInfo(info), showPlacedInfo(entry), showCellInfo(cell),
//   showHint(text): the info panel (showBlockInfo's info: {block, total,
//   generated, hint, enter, generate});
// - picker: the map's picker (mappick.js; a new one when missing), whose
//   other layers (phenomena, clouds) take a click or hover inside a
//   container before the blocks: their tooltip(entry), hover(entry) and
//   select(entry) are called; the drill-down adds its choices to it;
// - clearSelection(): a block was picked, so a phenomenon's ring goes;
// - canGenerate, courseSectors, sectorUrl(id), locate(name);
// - els: {crumbs, slabs, tooltip, notice, address, matches, controls}
//   (any may be missing).

const VERSION_QUERY = new URL(import.meta.url).search;

const S = await import(`./galaxystages.js${VERSION_QUERY}`);
const MC = await import(`./mapcontrol.js${VERSION_QUERY}`);
const { worldUnitsPerPixel } = await import(`./mapcore.js${VERSION_QUERY}`);
const B = await import(`./bookmarks.js${VERSION_QUERY}`);
const { createPicker, createTooltip } = await import(`./mappick.js${VERSION_QUERY}`);

// The other choices while one is hovered (MAP.18).
const OTHER_FADE = 0.25;
// A new stage's blocks fade in over the last part of a flight.
const FADE_IN_MS = 200;
const DRAG_CLICK_PX = 6;
// Each zoom step flies the camera to a preset for what it shows (MAP.97,
// Boss 2026-10-02: "when it zooms to a block it moves as isometric, wen
// it zooms to a slab it moves to top-down and the other direction as
// well. The only exception is the galaxy strtas out top-down so we can
// see the spiral arms."): the whole galaxy and a slab straight down, a
// block (several slabs: an arc, an entered block, the cube of sectors)
// at the isometric slant, so its slabs show side by side and can be
// picked on the map as well as with the slab buttons. Layers and blocks
// touch: no space between them (Boss, 2026-10-01).
const CUBE_MAX_SECTORS = 27;
export const ISO_TILT = Math.atan(Math.SQRT2);
// Every view can be turned, moved and zoomed (Boss, 2026-10-01; the
// whole galaxy too since MAP.85): drag turns it, right-drag or
// Shift-drag moves it, the wheel or a pinch zooms. It turns any way by
// any amount, through edge-on and under the plane, trackball style
// (MAP.96: "rotate the contents of the galaxy map any direction any
// amount"); zoom runs from MIN_ZOOM to MAX_ZOOM times the stage's own fit
// (the whole galaxy: GALAXY_MIN_ZOOM, about twice as close, out to its
// fit and no further, MAP.58), and the view's middle can't wander more
// than PAN_REACH fits away.
const ROTATE_PER_PX = (0.4 * Math.PI) / 180;
// Shift and an arrow key turn the view this far.
const KEY_TURN = (5 * Math.PI) / 180;
export const MIN_ZOOM = 1 / 8;
export const MAX_ZOOM = 2.5;
export const PAN_REACH = 1.5;
export const GALAXY_MIN_ZOOM = 0.5;
const WHEEL_ZOOM_PER_PX = 0.0025;
// The zoom policies (mapcontrol.js): a short range around the fit, and
// on the whole galaxy only closer than its fit.
const FREE_VIEW_ZOOM = MC.zoomPolicy(MC.ZOOM_RANGE, MIN_ZOOM, MAX_ZOOM);
const GALAXY_ZOOM = MC.zoomPolicy(MC.ZOOM_RANGE, GALAXY_MIN_ZOOM, 1);
// The arc under the pointer is outlined in full; its neighbors' outlines
// are this faint (MAP.85).
const NEIGHBOR_OUTLINE_OPACITY = 0.35;
// The lines between slabs while one is picked (MAP.77).
const SLAB_LINE_OPACITY = 0.45;
const TWO_PI = 2 * Math.PI;

// The free view's limits: the camera's distance from MIN_ZOOM to MAX_ZOOM times the stage's
// fit ...
export function clampZoom(dist, fitDist) {
  return MC.clampDistance(FREE_VIEW_ZOOM, dist, fitDist);
}

// ... and the view's middle no further than `reach` from the fit's
// (`target` is moved back in place).
export function clampPan(target, fitTarget, reach) {
  const off = [target[0] - fitTarget[0], target[1] - fitTarget[1], target[2] - fitTarget[2]];
  const far = Math.hypot(off[0], off[1], off[2]);
  if (far > reach) {
    for (let k = 0; k < 3; k++) target[k] = fitTarget[k] + (off[k] * reach) / far;
  }
  return target;
}

export function createStageView(host) {
  const THREE = host.THREE;
  const camera = host.camera;
  const canvasEl = host.canvasEl;
  const edgePc = host.edgePc;
  const els = host.els || {};
  const accent = new THREE.Color(host.accentColor || "#4f5fe8");

  let active = false;
  let outline = null;
  const dataCache = new Map();
  let stage = { at: null, picks: [] };
  let resolved = null;
  let selectedSector = null;
  let display = null;
  let leaving = [];
  let animation = null;
  // {option (index into display.options), sticky}
  let hover = null;
  let generatedOnly = false;
  let view = null;
  let goToken = 0;
  // The map's own history (Back and Forward buttons): this entry's index,
  // and the furthest one known.
  let mapIndex = 0;
  let maxIndex = 0;

  // --- The outline and the stage API's counts ----------------------------

  function getOutline() {
    if (!outline) {
      outline = S.galaxyOutline(edgePc, host.shape, host.galaxyRadius);
      outline.shapeless = !host.shape;
    }
    return outline;
  }

  function resolve(s) {
    return S.settleStage(s, getOutline(), edgePc);
  }

  function containerKey(at) {
    return at ? S.formatDrillKey(at) : "galaxy";
  }

  // The stage API's answer for container `at`, as {generated: Map(key ->
  // count), sectors: Map("ring/slot/layer" -> sector)}; null while it
  // loads (and on an error, which leaves counts at 0).
  function dataFor(at) {
    const cached = dataCache.get(containerKey(at));
    return cached && cached.ready ? cached.value : null;
  }

  function loadData(at) {
    const key = containerKey(at);
    let entry = dataCache.get(key);
    if (!entry) {
      entry = { ready: false, value: null };
      entry.promise = host.fetchStage(at ? S.stageQuery({ at: at, picks: [] }) : "").then(function (payload) {
        const generated = new Map();
        const looks = new Map();
        (payload.children || []).forEach(function (child) {
          generated.set(child.ring + "/" + child.wedge + "/" + child.slab, child.generated);
          looks.set(child.ring + "/" + child.wedge + "/" + child.slab, child.look || null);
        });
        const sectors = new Map();
        (payload.sectors || []).forEach(function (sector) {
          sectors.set(sector.ring + "/" + sector.slot + "/" + sector.layer, sector);
        });
        entry.value = { generated: generated, looks: looks, sectors: sectors, sectorGenerated: new Map() };
        payload.sectors && payload.sectors.forEach(function (sector) {
          entry.value.sectorGenerated.set(sector.ring + "/" + sector.slot + "/" + sector.layer, 1);
        });
        // A thin block's view is its sectors (galaxystages.thinSectors):
        // their own counts come from each level-3 child holding any.
        if (!at || !S.thinSectors(at, getOutline(), edgePc)) {
          entry.ready = true;
          return entry.value;
        }
        const filled = (payload.children || []).filter(function (child) { return child.generated > 0; });
        return Promise.all(filled.map(function (child) {
          const query = S.stageQuery({ at: { m: 3, ring: child.ring, wedge: child.wedge, slab: child.slab }, picks: [] });
          return host.fetchStage(query).then(function (inner) {
            (inner.sectors || []).forEach(function (sector) {
              const key = sector.ring + "/" + sector.slot + "/" + sector.layer;
              sectors.set(key, sector);
              entry.value.sectorGenerated.set(key, 1);
            });
          }, function () { /* that child's sectors read as not generated */ });
        })).then(function () {
          entry.ready = true;
          return entry.value;
        });
      }).then(null, function () {
        entry.value = { generated: new Map(), looks: new Map(), sectors: new Map(), sectorGenerated: new Map(), failed: true };
        entry.ready = true;
        dataCache.delete(key);
        return entry.value;
      });
      dataCache.set(key, entry);
    }
    return entry.promise;
  }

  // Forgets the generated counts (the galaxy changed); the shown stage
  // reloads its own.
  function invalidate() {
    dataCache.clear();
    if (active && display) {
      loadData(stage.at).then(function () {
        if (display && !animation) rebuildDisplay();
      });
    }
  }

  function generatedOf(block, data) {
    if (!data) return 0;
    const key = block.ring + "/" + block.wedge + "/" + block.slab;
    if (block.m === 1) return data.sectorGenerated.get(key) || 0;
    return data.generated.get(key) || 0;
  }

  // MAP.86: what a block's generated sectors hold, for its color and
  // opacity ({share, color, colored} from the stage API), or null.
  function lookOf(block, data) {
    if (!data) return null;
    const key = block.ring + "/" + block.wedge + "/" + block.slab;
    if (block.m === 1) {
      const sector = data.sectors.get(key);
      return sector ? sector.look || null : null;
    }
    return data.looks.get(key) || null;
  }

  function sumOf(blocks, data) {
    let generated = 0;
    let total = 0;
    blocks.forEach(function (block) {
      generated += generatedOf(block, data);
      total += block.total || 0;
    });
    return { generated: generated, total: total };
  }

  // The choices a stage offers: [{pick, blocks, a0, a1}], or, once there
  // is nothing left to pick, the view's blocks as one. Without an outline
  // (no shape yet), only blocks holding generated sectors.
  function choicesOf(r, data) {
    // A slab pick has one choice per slab, lowest first; the cube too: a
    // whole slab is lit and picked, never one of its sectors (MAP.91).
    let options = r.kind ? r.options : [{ pick: null, blocks: r.view.blocks, a0: r.view.a0, a1: r.view.a1 }];
    if (getOutline().shapeless) {
      options = options.map(function (o) {
        return Object.assign({}, o, { blocks: o.blocks.filter(function (b) { return generatedOf(b, data) > 0; }) });
      }).filter(function (o) { return o.blocks.length; });
    }
    return options;
  }

  function slabsIn(blocks) {
    return Array.from(new Set(blocks.map(function (b) { return b.slab; }))).sort(function (p, q) { return p - q; });
  }

  function isSectorView(r) {
    return !!(r && r.view && r.view.blocks.length && r.view.blocks[0].m === 1);
  }

  // A view of sectors across several layers, few enough to show as a
  // cube and pick one by one.
  function isCube(r) {
    return !!(r && r.kind === "layer" && isSectorView(r) && r.view.blocks.length <= CUBE_MAX_SECTORS);
  }

  // Whether the view can be turned and moved: every stage, since the
  // whole galaxy turns too (MAP.85).
  function isFree(r) {
    return !!(r && r.stage);
  }

  // The view's zoom policy: a short range round the fit; the whole galaxy
  // only zooms in from its own.
  function zoomPolicyFor(r) {
    return r && r.view && isWholeGalaxy(r) ? GALAXY_ZOOM : FREE_VIEW_ZOOM;
  }

  function isWholeGalaxy(r) {
    return !r.stage.at && r.view.a1 - r.view.a0 >= TWO_PI - 1e-9;
  }

  // The bounds a set of blocks covers in the plane, bearings counted
  // either way from `mid` (the middle of the view's bearings, so a block
  // reaching a little before the view's first bearing doesn't wrap all the
  // way round): {r0, r1, t0, t1, z0, z1}.
  function spanOf(blocks, mid) {
    const span = { r0: Infinity, r1: 0, t0: Infinity, t1: -Infinity, z0: Infinity, z1: -Infinity };
    blocks.forEach(function (block) {
      const b = block.bounds;
      const t0 = mid + wrapAngle(b.t0 - mid);
      span.r0 = Math.min(span.r0, b.r0);
      span.r1 = Math.max(span.r1, b.r1);
      span.t0 = Math.min(span.t0, t0);
      span.t1 = Math.max(span.t1, t0 + (b.t1 - b.t0));
      span.z0 = Math.min(span.z0, b.z0);
      span.z1 = Math.max(span.z1, b.z1);
    });
    return span;
  }

  // --- Drawing a stage ---------------------------------------------------------

  // One stage's meshes: a group per choice (so the hovered one stays and
  // the others fade), with each cell's choice behind the meshes for
  // picking.
  function buildDisplay(r) {
    const data = dataFor(r.stage.at);
    const options = choicesOf(r, data);
    const root = new THREE.Group();
    const groups = [];
    const eye = camera.position.toArray();
    options.forEach(function (option, index) {
      const cells = option.blocks.map(function (block) {
        const b = block.bounds;
        return {
          ring: block.ring, seg: block.wedge, slab: block.slab,
          r0: b.r0, r1: b.r1, t0: b.t0, t1: b.t1, z0: b.z0, z1: b.z1,
          filled: generatedOf(block, data), total: block.total, look: lookOf(block, data), block: block,
          option: index,
        };
      });
      const dim = generatedOnly ? function (cell) { return !(cell.filled > 0); } : null;
      const built = host.blockScene.buildCells(cells, eye, dim);
      const group = new THREE.Group();
      const meshes = [];
      // The whole galaxy shows no sector or block lines (MAP.85); while a
      // slab is picked the blocks show none either, only the lines between
      // slabs below (MAP.77); on one slab, the lines between its blocks,
      // the segments picked next.
      const gridEdges = isWholeGalaxy(r) || r.kind === "layer" ? 0 : 1;
      [built.solid, built.glass].forEach(function (part, n) {
        if (!part.vertexCount) return;
        const mesh = host.makeBlockMesh(part, n === 1);
        if (mesh.material.uniforms && mesh.material.uniforms.gridEdges) mesh.material.uniforms.gridEdges.value = gridEdges;
        mesh.userData.cells = built.cells[n];
        mesh.renderOrder = n;
        group.add(mesh);
        meshes.push(mesh);
      });
      root.add(group);
      groups.push({ option: option, meshes: meshes, fade: 1 });
    });
    const lines = r.kind === "layer" ? slabLines(options) : null;
    if (lines) root.add(lines);
    host.scene.add(root);
    root.updateMatrixWorld(true);
    return { resolved: r, root: root, options: options, groups: groups, fade: 1, data: data, lines: lines };
  }

  // The boundaries between slabs (MAP.77): each slab's outline along its
  // blocks' own sides, top and bottom, with its corners joined, and no
  // lines between the blocks inside it.
  function slabLines(options) {
    const points = [];
    options.forEach(function (option) {
      const mid = (option.a0 + option.a1) / 2;
      const span = spanOf(option.blocks, mid);
      Array.prototype.push.apply(points, edgePoints(S.outlineEdges(option.blocks, mid), span.z1, span));
    });
    if (!points.length) return null;
    const line = new THREE.LineSegments(
      new THREE.BufferGeometry().setFromPoints(points),
      new THREE.LineBasicMaterial({ color: accent, transparent: true, opacity: SLAB_LINE_OPACITY, depthWrite: false }),
    );
    line.renderOrder = 3;
    line.frustumCulled = false;
    return line;
  }

  function disposeDisplay(d) {
    if (!d) return;
    host.scene.remove(d.root);
    if (d.lines) {
      d.lines.geometry.dispose();
      d.lines.material.dispose();
    }
    d.groups.forEach(function (group) {
      group.meshes.forEach(function (mesh) {
        mesh.geometry.dispose();
        mesh.material.dispose();
      });
    });
  }

  function setDisplayFade(d, fade) {
    d.fade = fade;
    if (d.lines) d.lines.material.opacity = SLAB_LINE_OPACITY * fade;
    d.groups.forEach(function (group) {
      group.meshes.forEach(function (mesh) {
        const value = fade * group.fade;
        mesh.material.uniforms.fade.value = value;
        mesh.material.transparent = mesh.renderOrder === 1 || value < 1;
        mesh.material.depthWrite = mesh.renderOrder === 0 && value >= 1;
      });
    });
  }

  function rebuildDisplay() {
    const old = display;
    display = buildDisplay(resolved);
    setDisplayFade(display, 1);
    disposeDisplay(old);
    applyHover();
    renderStrip();
  }

  // --- Framing -----------------------------------------------------------------

  function fovHalf() {
    const vertical = THREE.MathUtils.degToRad(camera.fov) / 2;
    const aspect = camera.aspect || 1;
    return Math.min(vertical, Math.atan(Math.tan(vertical) * aspect));
  }

  // Points round blocks' bounds, for fitting the view to them: each
  // block's corners, and along its outer arc every 22.5 degrees, top and
  // bottom.
  function boundsPoints(blocks) {
    const points = [];
    blocks.forEach(function (block) {
      const b = block.bounds;
      const steps = Math.max(1, Math.ceil((b.t1 - b.t0) / (Math.PI / 8)));
      [b.z0, b.z1].forEach(function (z) {
        points.push(b.r0 * Math.cos(b.t0), b.r0 * Math.sin(b.t0), z, b.r0 * Math.cos(b.t1), b.r0 * Math.sin(b.t1), z);
        for (let k = 0; k <= steps; k++) {
          const t = b.t0 + ((b.t1 - b.t0) * k) / steps;
          points.push(b.r1 * Math.cos(t), b.r1 * Math.sin(t), z);
        }
      });
    });
    return points;
  }

  // The points a stage's view is fitted round, worked out once per stage.
  const fitPointsCache = new WeakMap();
  function fitPointsOf(r) {
    let points = fitPointsCache.get(r);
    if (!points) {
      points = boundsPoints(r.view.blocks);
      fitPointsCache.set(r, points);
    }
    return points;
  }

  const fitBasis = { m: null, x: null, y: null, z: null };

  // The camera distance from `target` along turn `quat`'s back axis at which every
  // point (a flat [x, y, z, ...] list) shows on the map with FIT_MARGIN to
  // spare, across the map's actual width and height (MAP.53, MAP.78): a
  // wide map fits a long arc by its width, a tall one by its height, and a
  // bigger map shows the same fit larger.
  function fitDistance(points, target, quat) {
    if (!fitBasis.m) {
      fitBasis.m = new THREE.Matrix4();
      fitBasis.x = new THREE.Vector3();
      fitBasis.y = new THREE.Vector3();
      fitBasis.z = new THREE.Vector3();
    }
    fitBasis.m.makeRotationFromQuaternion(quat);
    fitBasis.m.extractBasis(fitBasis.x, fitBasis.y, fitBasis.z);
    const X = fitBasis.x;
    const Y = fitBasis.y;
    const Z = fitBasis.z;
    const tanV = Math.tan(THREE.MathUtils.degToRad(camera.fov) / 2);
    const aspect = canvasEl.clientWidth && canvasEl.clientHeight ? canvasEl.clientWidth / canvasEl.clientHeight : camera.aspect || 1;
    const tanH = tanV * aspect;
    const margin = S.FIT_MARGIN;
    let dist = 0;
    for (let k = 0; k + 2 < points.length; k += 3) {
      const vx = points[k] - target[0];
      const vy = points[k + 1] - target[1];
      const vz = points[k + 2] - target[2];
      const x = vx * X.x + vy * X.y + vz * X.z;
      const y = vx * Y.x + vy * Y.y + vz * Y.z;
      const z = vx * Z.x + vy * Z.y + vz * Z.z;
      dist = Math.max(dist, z + (margin * Math.abs(x)) / tanH, z + (margin * Math.abs(y)) / tanV);
    }
    return dist > 0 ? dist : host.galaxyRadius * 2.4;
  }

  // fitDistance with the target moved (across the screen, three rounds)
  // so the points sit in the middle of the map: a slanted view sees the
  // near side of a block bigger than the far side, so the middle of the
  // blocks isn't the middle of their picture. {target, dist}.
  function centeredFit(points, target, quat) {
    let at = target.slice();
    let dist = fitDistance(points, at, quat);
    const tanV = Math.tan(THREE.MathUtils.degToRad(camera.fov) / 2);
    const aspect = canvasEl.clientWidth && canvasEl.clientHeight ? canvasEl.clientWidth / canvasEl.clientHeight : camera.aspect || 1;
    const tanH = tanV * aspect;
    for (let round = 0; round < 3; round++) {
      const X = fitBasis.x;
      const Y = fitBasis.y;
      const Z = fitBasis.z;
      let x0 = Infinity;
      let x1 = -Infinity;
      let y0 = Infinity;
      let y1 = -Infinity;
      for (let k = 0; k + 2 < points.length; k += 3) {
        const vx = points[k] - at[0];
        const vy = points[k + 1] - at[1];
        const vz = points[k + 2] - at[2];
        const depth = dist - (vx * Z.x + vy * Z.y + vz * Z.z);
        if (depth <= 0) continue;
        const sx = (vx * X.x + vy * X.y + vz * X.z) / (depth * tanH);
        const sy = (vx * Y.x + vy * Y.y + vz * Y.z) / (depth * tanV);
        x0 = Math.min(x0, sx);
        x1 = Math.max(x1, sx);
        y0 = Math.min(y0, sy);
        y1 = Math.max(y1, sy);
      }
      if (!(x1 >= x0)) break;
      const cx = ((x0 + x1) / 2) * dist * tanH;
      const cy = ((y0 + y1) / 2) * dist * tanV;
      if (Math.abs(cx) + Math.abs(cy) < 1e-6 * dist) break;
      at = [at[0] + X.x * cx + Y.x * cy, at[1] + X.y * cx + Y.y * cy, at[2] + X.z * cx + Y.z * cy];
      dist = fitDistance(points, at, quat);
    }
    return { target: at, dist: dist };
  }

  // The camera's turn looking from bearing `theta` (radians from +x) and
  // tilt `phi` (from straight down) at its target, with the screen's up
  // pointing away from the camera across the plane (galactic north up
  // for theta -90 degrees); straight down needs no special case.
  function presetQuat(theta, phi) {
    const back = new THREE.Vector3(Math.sin(phi) * Math.cos(theta), Math.sin(phi) * Math.sin(theta), Math.cos(phi));
    const upward = new THREE.Vector3(-Math.cos(phi) * Math.cos(theta), -Math.cos(phi) * Math.sin(theta), Math.sin(phi));
    const right = new THREE.Vector3().crossVectors(upward, back);
    return new THREE.Quaternion().setFromRotationMatrix(new THREE.Matrix4().makeBasis(right, upward, back));
  }

  // The preset tilt for a stage (MAP.97): straight down for the whole
  // galaxy and for a slab (one layer left: its segments, or its sectors),
  // the isometric slant for a block of several slabs.
  function presetTilt(r) {
    if (isWholeGalaxy(r)) return 0;
    return r.kind === "layer" ? ISO_TILT : 0;
  }

  // The camera for a stage, with the view's middle bearing pointing up
  // the screen (galactic north up for the whole galaxy), at the stage's
  // preset tilt, turned about the middle of what it shows and fitted
  // round its blocks: {target, dist, quat}.
  function cameraFor(r) {
    const blocks = r.view.blocks;
    if (!blocks.length) {
      return { target: [0, 0, 0], dist: host.galaxyRadius * 2.4, quat: presetQuat(-Math.PI / 2, 0) };
    }
    const fp = S.footprint(blocks);
    let z0 = Infinity;
    let z1 = -Infinity;
    blocks.forEach(function (block) {
      z0 = Math.min(z0, block.bounds.z0);
      z1 = Math.max(z1, block.bounds.z1);
    });
    const theta = isWholeGalaxy(r) ? -Math.PI / 2 : (r.view.a0 + r.view.a1) / 2 + Math.PI;
    const quat = presetQuat(theta, presetTilt(r));
    const fit = centeredFit(fitPointsOf(r), [fp.center[0], fp.center[1], (z0 + z1) / 2], quat);
    return { target: fit.target, dist: fit.dist, quat: quat };
  }

  // The view on arrival at camera `to`, keeping `to` as its fit. The
  // target is a copy: a pan moves view.target in place, and the pan's
  // reach is measured from the fit's. `zoom` is the user's zoom, the
  // distance over the fit's at the view's own turn.
  function settledView(to) {
    return Object.assign({}, to, { target: to.target.slice(), quat: to.quat.clone(), fit: to, zoom: 1 });
  }

  // The view as a flight's start: {target, dist, quat}, copied.
  function viewNow() {
    return { target: view.target.slice(), dist: view.dist, quat: view.quat.clone() };
  }

  // Turns the view about its target, trackball style (MAP.96): `yaw`
  // about the screen's up axis, `pitch` about its right axis, radians
  // (a drag right or down by that much); no limit and no flip at the
  // poles, since the turn is kept as a quaternion.
  const turnAxis = { x: null, y: null, q: null };
  function turnView(yaw, pitch) {
    if (!turnAxis.q) {
      turnAxis.x = new THREE.Vector3(1, 0, 0);
      turnAxis.y = new THREE.Vector3(0, 1, 0);
      turnAxis.q = new THREE.Quaternion();
    }
    view.quat.multiply(turnAxis.q.setFromAxisAngle(turnAxis.y, -yaw));
    view.quat.multiply(turnAxis.q.setFromAxisAngle(turnAxis.x, -pitch));
    view.quat.normalize();
  }

  // Where the stage's blocks fall on the map now, in canvas pixels:
  // {left, top, right, bottom, width, height} (the map's own size). Read
  // by the browser tests through the canvas.
  canvasEl.galaxyFrame = function () {
    if (!resolved || !resolved.view) return null;
    camera.updateMatrixWorld();
    const points = fitPointsOf(resolved);
    const v = new THREE.Vector3();
    const w = canvasEl.clientWidth;
    const h = canvasEl.clientHeight;
    const box = { left: Infinity, top: Infinity, right: -Infinity, bottom: -Infinity, width: w, height: h };
    for (let k = 0; k + 2 < points.length; k += 3) {
      v.set(points[k], points[k + 1], points[k + 2]).project(camera);
      const x = ((v.x + 1) / 2) * w;
      const y = ((1 - v.y) / 2) * h;
      box.left = Math.min(box.left, x);
      box.right = Math.max(box.right, x);
      box.top = Math.min(box.top, y);
      box.bottom = Math.max(box.bottom, y);
    }
    return box;
  };

  // Which lines the stage draws (MAP.77): {blockEdges} (the blocks' own
  // edges on or off) and {slabLines} (how many line pieces trace the
  // boundaries between slabs). Read by the browser tests.
  canvasEl.galaxyLines = function () {
    if (!display) return null;
    let blockEdges = 0;
    display.groups.forEach(function (group) {
      group.meshes.forEach(function (mesh) {
        if (mesh.material.uniforms && mesh.material.uniforms.gridEdges) blockEdges = Math.max(blockEdges, mesh.material.uniforms.gridEdges.value);
      });
    });
    const slabLines = display.lines ? display.lines.geometry.getAttribute("position").count / 2 : 0;
    return { blockEdges: blockEdges, slabLines: slabLines, kind: display.resolved.kind };
  };

  // Each slab button's slab outline on the map, as line pieces in client
  // pixels cut to the map ({slab: [[[x, y], [x, y]], ...]}), and the
  // strip's mode (MAP.98, MAP.99). Read by the browser tests.
  canvasEl.galaxySlabOutlines = function () {
    if (!strip || !view) return null;
    const rect = canvasEl.getBoundingClientRect();
    const out = { mode: strip.mode, slabs: {} };
    strip.rows.forEach(function (r) { out.slabs[r.option.pick.lo] = outlineOnScreen(r.outline, rect); });
    return out;
  };

  // The camera's tilt from straight down, in degrees (0 top-down, 90
  // edge-on, 180 from under the plane), and whether a flight is running.
  // Read by the browser tests.
  canvasEl.galaxyCamera = function () {
    if (!view) return null;
    const back = new THREE.Vector3(0, 0, 1).applyQuaternion(view.quat);
    return { tilt: THREE.MathUtils.radToDeg(Math.acos(Math.max(-1, Math.min(1, back.z)))), flying: !!animation };
  };

  // Keeps the view fitted round the stage as it turns or the map changes
  // size (MAP.78): the distance is the user's zoom times the fit at the
  // view's own turn.
  function reframe() {
    if (!view || !view.fit || !resolved || !resolved.view || !resolved.view.blocks.length) return;
    view.dist = view.zoom * fitDistance(fitPointsOf(resolved), view.fit.target, view.quat);
  }

  function applyView() {
    host.setCamera({ target: view.target.slice(), dist: view.dist, quaternion: view.quat.toArray() });
    drawLeaders();
  }

  function wrapAngle(a) {
    return Math.atan2(Math.sin(a), Math.cos(a));
  }

  // Animates the camera from `from` to `to` along van Wijk and Nuij's
  // path in the plane, easing the height and turn (the shortest way
  // between the two turns) alongside; calls
  // progress(t) each frame (t 0..1) and done() at the end.
  function flyCamera(from, to, progress, done) {
    const tanHalf = Math.tan(THREE.MathUtils.degToRad(camera.fov) / 2);
    const path = S.flightPath(from.target, 2 * from.dist * tanHalf, to.target, 2 * to.dist * tanHalf);
    animation = {
      startedAt: performance.now(),
      duration: host.reducedMotion ? 0 : S.flightMs(path.S),
      frame: function (t) {
        const e = S.easeInOut(t);
        const p = path.at(path.S * t);
        view = {
          target: [p.center[0], p.center[1], from.target[2] + (to.target[2] - from.target[2]) * e],
          dist: p.w / (2 * tanHalf),
          quat: new THREE.Quaternion().slerpQuaternions(from.quat, to.quat, e),
        };
        applyView();
        if (progress) progress(t);
      },
      done: done,
    };
  }

  // --- Moving between stages -------------------------------------------------

  function setWedgeClip(r) {
    if (!host.setWedgeClip) return;
    if (!r.view.blocks.length || isWholeGalaxy(r)) {
      host.setWedgeClip(null);
      return;
    }
    const span = spanOf(r.view.blocks, (r.view.a0 + r.view.a1) / 2);
    host.setWedgeClip({
      r0: span.r0, r1: span.r1, a0: span.t0, a1: span.t1, z0: span.z0, z1: span.z1,
      cells: r.view.blocks.map(function (block) { return block.bounds; }),
    });
  }

  // A stage's query with the NAV pick kept (host.pickQuery, "?pick=...").
  function withPick(query) {
    if (!host.pickQuery) return query;
    return query ? query + "&" + host.pickQuery.slice(1) : host.pickQuery;
  }

  // Goes to stage `next` (carried on through any choice of one). push:
  // record it in the browser history (default yes).
  function go(next, options) {
    options = options || {};
    if (animation) finishAnimation();
    const r = resolve(next);
    if (r.problem) {
      notice(r.problem);
      return;
    }
    notice("");
    if (!options.keepSector) selectedSector = null;
    if (r.sector && !selectedSector) {
      selectedSector = { ring: r.sector.ring, layer: r.sector.slab, slot: r.sector.wedge };
    }
    const token = ++goToken;
    const query = withPick(options.query || S.stageQuery(r.stage));
    if (options.push !== false && (!S.sameStage(stage, r.stage) || query !== location.search)) {
      mapIndex += 1;
      maxIndex = mapIndex;
      history.pushState({ galaxyStage: true, mapIndex: mapIndex, maxIndex: maxIndex }, "", location.pathname + query + location.hash);
    }
    stage = r.stage;
    resolved = r;
    renderCrumbs();
    renderStrip();
    renderTravel();
    hover = null;
    showTooltip("", 0, 0);
    loadData(r.stage.at).then(function () {
      if (!active || token !== goToken) return;
      if (animation) finishAnimation();
      flyTo(r);
    });
  }

  function flyTo(r) {
    const from = display && view ? viewNow() : null;
    const to = cameraFor(r);
    const old = display;
    const incoming = buildDisplay(r);
    setWedgeClip(r);
    setDisplayFade(incoming, from ? 0 : 1);
    display = incoming;
    if (!from) {
      view = settledView(to);
      applyView();
      afterArrival();
      return;
    }
    if (old) leaving.push(old);
    flyCamera(from, to, function (t) {
      const ms = animation ? animation.duration : 0;
      const fadeStart = ms > 0 ? Math.max(0, 1 - FADE_IN_MS / ms) : 0;
      const f = fadeStart >= 1 ? 1 : Math.max(0, (t - fadeStart) / (1 - fadeStart));
      setDisplayFade(incoming, f);
      if (old) setDisplayFade(old, 1 - f);
    }, function () {
      setDisplayFade(display, 1);
      view = settledView(to);
      applyView();
      afterArrival();
    });
  }

  function afterArrival() {
    leaving.forEach(function (d) { if (d && d !== display) disposeDisplay(d); });
    leaving = [];
    hover = null;
    renderStrip();
    const index = selectedSector ? sectorOption(selectedSector) : -1;
    if (index >= 0) {
      setHover({ option: index, sticky: true });
      showSectorInfo(display.options[index].blocks[0]);
    } else {
      applyHover();
      showStageInfo();
    }
  }

  function finishAnimation() {
    const a = animation;
    animation = null;
    a.frame(1);
    a.done();
  }

  // Every frame: the running move, if any.
  function step(now) {
    if (!active || !animation) return;
    const a = animation;
    const t = a.duration > 0 ? Math.min(1, (now - a.startedAt) / a.duration) : 1;
    a.frame(t);
    if (t >= 1) {
      animation = null;
      a.done();
    }
  }

  function up() {
    if (selectedSector && !resolved.sector) {
      selectedSector = null;
      go(stage);
      return;
    }
    const parent = S.parentStage(stage, getOutline(), edgePc);
    if (parent) go(parent);
  }

  function home() {
    go({ at: null, picks: [] });
  }

  // The map's Back and Forward: the browser's own history, kept to this
  // map's stages.
  function travel(direction) {
    if (direction < 0 && mapIndex > 0) history.back();
    else if (direction > 0 && mapIndex < maxIndex) history.forward();
  }

  function renderTravel() {
    if (!els.controls) return;
    const back = els.controls.querySelector('[data-action="back"]');
    const forward = els.controls.querySelector('[data-action="forward"]');
    const upButton = els.controls.querySelector('[data-action="up"]');
    if (back) back.disabled = mapIndex <= 0;
    if (forward) forward.disabled = mapIndex >= maxIndex;
    if (upButton) upButton.disabled = !stage.at && !stage.picks.length && !selectedSector;
    const resetButton = els.controls.querySelector('[data-action="reset-view"]');
    if (resetButton) resetButton.disabled = !isFree(resolved);
  }

  // --- Picking and hover -----------------------------------------------------

  // The stage's choices are the picker's last layer (mappick.js): a
  // phenomenon the map puts in front takes the pointer first.
  const picker = host.picker || createPicker(camera, canvasEl);
  const tooltip = createTooltip(els.tooltip);
  const choiceLayer = picker.addLayer({
    name: "choices", priority: 10,
    enabled: function () { return !!display && !!resolved; },
    meshes: function () {
      const meshes = [];
      display.groups.forEach(function (group) { Array.prototype.push.apply(meshes, group.meshes); });
      return meshes;
    },
    entryOf: function (hit) {
      const cell = hit.object.userData.cells[hit.object.userData.part.owners[hit.face.a]];
      return cell ? cell.option : null;
    },
  });

  // Whether choice `index` can be picked: not while it holds nothing
  // generated with "Generated only" on.
  function pickable(index) {
    if (!display || !display.options[index]) return false;
    if (!generatedOnly) return true;
    const data = display.data;
    return display.options[index].blocks.some(function (b) { return generatedOf(b, data) > 0; });
  }

  // The choice under a screen point (an index into display.options), or
  // -1; while the pick is a slab, the slab of the block under it. Every
  // choice lights on hover, with "Generated only" on (always on while
  // picking a NAV end) too: only taking one needs something generated
  // (pickable, NAV.31).
  function optionAt(clientX, clientY) {
    const found = picker.pick(clientX, clientY, function (layer) { return layer === choiceLayer; });
    return found ? found.entry : -1;
  }

  let outlineGroup = null;

  function clearOutline() {
    if (outlineGroup) {
      host.scene.remove(outlineGroup);
      outlineGroup.children.forEach(function (line) {
        line.geometry.dispose();
        line.material.dispose();
      });
      outlineGroup = null;
    }
  }

  // Line pieces [x, y, z] pairs tracing `edges` (galaxystages.outlineEdges)
  // at height z; with `walls` ({z0, z1}), at both heights and up the side
  // from z0 to z1 at each corner where a radial edge turns along a circle
  // (an arc runs the disk's whole height).
  function edgePoints(edges, z, walls) {
    const points = [];
    const at = function (r, t, h) { return new THREE.Vector3(r * Math.cos(t), r * Math.sin(t), h); };
    const heights = walls ? [walls.z0, walls.z1] : [z];
    heights.forEach(function (h) {
      edges.radial.forEach(function (e) { points.push(at(e.r0, e.t, h), at(e.r1, e.t, h)); });
      edges.circles.forEach(function (e) {
        const steps = Math.max(1, Math.ceil((e.t1 - e.t0) / (Math.PI / 90)));
        for (let k = 0; k < steps; k++) {
          points.push(at(e.r, e.t0 + ((e.t1 - e.t0) * k) / steps, h), at(e.r, e.t0 + ((e.t1 - e.t0) * (k + 1)) / steps, h));
        }
      });
    });
    if (walls) {
      // A corner is a radial edge's end no other radial edge carries on
      // from (a straight side running across a ring boundary has none).
      const ends = new Map();
      edges.radial.forEach(function (e) {
        [e.r0, e.r1].forEach(function (r) {
          const key = Math.round(r * 1e3) + "/" + Math.round(Math.cos(e.t) * 1e9) + "/" + Math.round(Math.sin(e.t) * 1e9);
          const end = ends.get(key) || { r: r, t: e.t, n: 0 };
          end.n += 1;
          ends.set(key, end);
        });
      });
      ends.forEach(function (end) {
        if (end.n === 1 && end.r > 0) points.push(at(end.r, end.t, walls.z0), at(end.r, end.t, walls.z1));
      });
    }
    return points;
  }

  function addOutline(points, opacity) {
    if (!outlineGroup) {
      outlineGroup = new THREE.Group();
      outlineGroup.renderOrder = 6;
      host.scene.add(outlineGroup);
    }
    const line = new THREE.LineSegments(
      new THREE.BufferGeometry().setFromPoints(points),
      new THREE.LineBasicMaterial({ color: accent, depthTest: false, transparent: true, opacity: opacity }),
    );
    line.renderOrder = 6;
    line.frustumCulled = false;
    outlineGroup.add(line);
  }

  // An accent outline around a choice's area along its blocks' own sides
  // (MAP.18, MAP.52): on top of it, or for an arc of the galaxy round its
  // whole height, with its neighboring arcs' outlines faint (MAP.85).
  function outlineOption(index) {
    clearOutline();
    const option = display.options[index];
    const mid = (option.a0 + option.a1) / 2;
    const span = spanOf(option.blocks, mid);
    const isArc = option.pick && option.pick.kind === "arc";
    // An arc or a slab is outlined round its whole height (MAP.91).
    const walls = isArc || (option.pick && option.pick.kind === "layer");
    addOutline(edgePoints(S.outlineEdges(option.blocks, mid), span.z1, walls ? span : null), 1);
    if (!isArc) return;
    neighborArcs(option.pick.n).forEach(function (n) {
      const other = display.options.find(function (o) { return o.pick && o.pick.kind === "arc" && o.pick.n === n; });
      if (!other) return;
      const otherMid = (other.a0 + other.a1) / 2;
      const otherSpan = spanOf(other.blocks, otherMid);
      addOutline(edgePoints(S.outlineEdges(other.blocks, otherMid), otherSpan.z1, null), NEIGHBOR_OUTLINE_OPACITY);
    });
  }

  // The arcs beside arc n: either side in its band, and in and out.
  function neighborArcs(n) {
    const per = S.ARCS_PER_TURN;
    const band = Math.floor(n / per);
    const bin = n % per;
    const out = [band * per + ((bin + 1) % per), band * per + ((bin + per - 1) % per)];
    if (band > 0) out.push(n - per);
    if (band < S.ARC_BANDS - 1) out.push(n + per);
    return out;
  }

  function setHover(next) {
    hover = next;
    applyHover();
  }

  // The hovered choice stays as it is and the others fade (a layer
  // hovered in the strip too); a choice on the map gets an outline, a
  // slab round all its blocks (MAP.91).
  function applyHover() {
    if (!display) return;
    const index = hover && hover.option != null ? hover.option : -1;
    const layer = hover && hover.layer ? hover.layer : null;
    display.groups.forEach(function (group, n) {
      const slab = group.option.blocks[0].slab;
      const lit = layer ? slab >= layer.lo && slab <= layer.hi : index < 0 || n === index;
      group.fade = lit ? 1 : OTHER_FADE;
    });
    if (!animation) setDisplayFade(display, 1);
    const option = index >= 0 ? display.options[index] : null;
    if (option) {
      outlineOption(index);
    } else {
      clearOutline();
    }
    let row = layer;
    if (!row && option && resolved.kind === "layer" && option.pick) row = { lo: option.pick.lo, hi: option.pick.hi };
    markStripRow(row);
  }

  function showTooltip(text, clientX, clientY) {
    tooltip.show(text, clientX, clientY);
  }

  // The tooltip beside a choice picked by keyboard: at its middle on
  // screen.
  function tooltipAtOption(index) {
    const option = display.options[index];
    const fp = S.footprint(option.blocks);
    tooltip.showAt(optionText(index), camera, canvasEl, fp.center[0], fp.center[1], option.blocks[0].bounds.z1);
  }

  function blockText(block, data) {
    const b = block.bounds;
    const bearing = S.bearingRange(b);
    const generated = generatedOf(block, data);
    const parts = [(block.m === 1 ? "Sector " : "Block ") + S.blockLabel(block),
      "bearing " + bearing[0].toFixed(1) + "°–" + bearing[1].toFixed(1) + "°",
      Math.round(b.r0) + "–" + Math.round(b.r1) + " pc from the core"];
    if (block.m === 1) parts.push(generated > 0 ? "generated" : "not generated");
    else parts.push(S.formatCount(generated) + (block.total != null ? " of " + S.formatCount(block.total) : "") + " sectors generated");
    return parts.join(", ");
  }

  function layerText(pick) {
    const first = S.slabLayers(stage.at, pick.lo).first;
    const last = S.slabLayers(stage.at, pick.hi).last;
    const z0 = (first - 0.5) * edgePc;
    const z1 = (last + 0.5) * edgePc;
    const where = z1 <= 0 ? Math.round(-z1) + "–" + Math.round(-z0) + " pc below the plane"
      : z0 >= 0 ? Math.round(z0) + "–" + Math.round(z1) + " pc above the plane"
        : "across the plane";
    const span = first === last ? "layer " + first : "layers " + first + " to " + last;
    return S.pickLabel(pick, stage.at, resolved.view) + ": " + span + ", " + where;
  }

  function optionText(index) {
    return choiceText(index) + (pickable(index) ? "" : " (nothing generated here to pick)");
  }

  function choiceText(index) {
    const option = display.options[index];
    const data = display.data;
    if (option.blocks.length === 1) return blockText(option.blocks[0], data);
    const sum = sumOf(option.blocks, data);
    const counts = S.formatCount(sum.generated) + (getOutline().shapeless ? "" : " of " + S.formatCount(sum.total)) + " sectors generated";
    if (option.pick && option.pick.kind === "layer") return layerText(option.pick) + ", " + counts;
    const span = spanOf(option.blocks, (option.a0 + option.a1) / 2);
    return (option.pick ? S.pickLabel(option.pick, stage.at, option) : "Here") + ", "
      + Math.round(span.r0) + "–" + Math.round(span.r1) + " pc from the core, " + counts;
  }

  // The panel on arrival: the stage's own block (its hint under it), or
  // just the hint at the galaxy.
  function showStageInfo() {
    const hint = hintFor(resolved);
    if (!stage.at) {
      host.showHint(hint);
      return;
    }
    const data = dataFor(stage.at);
    let generated = 0;
    if (data) data.generated.forEach(function (n) { generated += n; });
    const total = getOutline().shapeless ? null : S.drillBlockTotal(stage.at, getOutline());
    host.showBlockInfo({
      block: stage.at, total: total,
      generated: data ? generated : null, hint: hint, enter: null,
      generate: generateOffer(total, data ? generated : null),
    });
  }

  // In a 3-block, for an admin: {block, layer, wholeBlock} for the
  // "Generate this layer" and "Generate this block" buttons
  // (generatebuttons.js's blockGenerateButtons), `layer` set while one
  // layer of sectors is shown and not all generated. Null when there is
  // nothing to offer.
  function generateOffer(total, generated) {
    if (!host.canGenerate || !stage.at || stage.at.m !== 3) return null;
    const blockDone = total != null && generated != null && total > 0 && generated >= total;
    const blocks = resolved.view.blocks;
    let layer = isSectorView(resolved) && blocks.every(function (b) { return b.slab === blocks[0].slab; }) ? blocks[0].slab : null;
    if (layer != null && !getOutline().shapeless) {
      const sum = sumOf(blocks, dataFor(stage.at));
      if (sum.total > 0 && sum.generated >= sum.total) layer = null;
    }
    if (blockDone && layer == null) return null;
    return { block: S.formatDrillKey(stage.at), layer: layer, wholeBlock: !blockDone };
  }

  function sectorEntry(block, sector) {
    const c = block.bounds;
    const r = (c.r0 + c.r1) / 2;
    const t = (c.t0 + c.t1) / 2;
    const z = (c.z0 + c.z1) / 2;
    return {
      id: sector.id, name: sector.name, system_count: sector.system_count,
      ring_index: block.ring, layer_index: block.slab, ring_slot_index: block.wedge,
      designation: S.sectorDesignation(block.ring, block.slab, block.wedge),
      galactic_radius_pc: Math.hypot(r * Math.cos(t), r * Math.sin(t), z),
    };
  }

  function sectorRecord(block) {
    const data = dataFor(stage.at);
    return data ? data.sectors.get(block.ring + "/" + block.wedge + "/" + block.slab) : null;
  }

  function showSectorInfo(block) {
    const sector = sectorRecord(block);
    if (sector) host.showPlacedInfo(sectorEntry(block, sector));
    else host.showCellInfo({ m: 1, bounds: block.bounds, address: { ring: block.ring, layer: block.slab, slot: block.wedge }, filled: 0 });
  }

  // A stage's own hint in the info panel.
  function hintFor(r) {
    const base = baseHint(r);
    if (!isFree(r)) return base;
    return (base ? base + " " : "") + "Drag (or Shift and the arrow keys) to turn the view about its middle, right-drag (or "
      + "Shift-drag) to move it, scroll or pinch to zoom; Reset view brings it back.";
  }

  function baseHint(r) {
    if (!r || !r.kind) return "";
    if (isCube(r)) {
      return "Click a layer of sectors on the map, or pick one with the buttons beside the map, to see just that layer.";
    }
    if (r.kind === "layer") {
      return "Click a " + S.slabNoun(r.stage.at).toLowerCase() + " (a layer of the disk) on the map, or pick one with the buttons beside the map; each button's line points at its slab.";
    }
    if (r.kind === "arc") return "Click an arc of the galaxy (a piece of the disk, top to bottom) to look at it more closely.";
    if (isSectorView(r)) {
      return "Click a sector to open it; a sector that isn't generated yet shows where it is"
        + (host.canGenerate ? " and how to generate it." : ".");
    }
    return "Click a block of the slab to zoom into it.";
  }

  // --- Acting on a choice ----------------------------------------------------

  function act(index) {
    if (!pickable(index)) return;
    const option = display.options[index];
    // One sector is opened, unless it is all of a slab to pick (MAP.91).
    const slabPick = option.pick && option.pick.kind === "layer";
    if (!slabPick && option.blocks.length === 1 && option.blocks[0].m === 1) {
      const block = option.blocks[0];
      const sector = sectorRecord(block);
      if (sector && host.sectorUrl(sector.id)) {
        window.location.assign(host.sectorUrl(sector.id));
        return;
      }
      selectedSector = { ring: block.ring, layer: block.slab, slot: block.wedge };
      setHover({ option: index, sticky: true });
      showSectorInfo(block);
      renderCrumbs();
      renderTravel();
      return;
    }
    if (!option.pick) return;
    go({ at: stage.at, picks: stage.picks.concat([option.pick]) });
  }

  function sectorOption(sector) {
    if (!display) return -1;
    return display.options.findIndex(function (option) {
      return option.blocks.length === 1 && option.blocks[0].m === 1 && option.blocks[0].ring === sector.ring
        && option.blocks[0].wedge === sector.slot && option.blocks[0].slab === sector.layer;
    });
  }

  // --- Input -------------------------------------------------------------------

  let touchPending = -1;
  // The zoom a pinch started from.
  let pinchZoom = 1;

  // Turns the view (drag) or moves it in the screen's plane (pan).
  function drag(dx, dy, pan) {
    if (!pan) {
      turnView(dx * ROTATE_PER_PX, dy * ROTATE_PER_PX);
      reframe();
    } else {
      const perPx = worldUnitsPerPixel(camera, view.dist, canvasEl.clientHeight);
      const fit = view.fit;
      const reach = PAN_REACH * fit.dist * Math.tan(fovHalf());
      MC.panInScreenPlane(THREE, camera, view.target, dx, dy, perPx);
      clampPan(view.target, fit.target, reach);
    }
    applyView();
  }

  function zoomBy(factor) {
    if (!view || animation || !MC.canZoom(zoomPolicyFor(resolved))) return;
    view.zoom = MC.clampDistance(zoomPolicyFor(resolved), view.zoom * factor, 1);
    reframe();
    applyView();
  }

  // The wheel zooms within the view's zoom policy (a locked view leaves
  // it to scroll the page). True when it was used.
  function onWheel(event) {
    if (!MC.canZoom(zoomPolicyFor(resolved)) || !view) return false;
    const deltaPx = MC.wheelPixels(event, canvasEl.clientHeight, 200);
    if (deltaPx) zoomBy(Math.exp(deltaPx * WHEEL_ZOOM_PER_PX));
    return true;
  }

  // Back to the stage's own view, after turning or moving it.
  function resetView() {
    if (!resolved || animation || !display) return;
    const from = viewNow();
    const to = cameraFor(resolved);
    flyCamera(from, to, null, function () {
      view = settledView(to);
      applyView();
    });
  }

  // The other layer (a phenomenon) under the pointer, while there is one.
  let hoveredOther = null;

  function hoverOther(found) {
    const layer = found && found.layer !== choiceLayer ? found.layer : null;
    if (hoveredOther && hoveredOther !== layer && hoveredOther.hover) hoveredOther.hover(null);
    hoveredOther = layer;
    if (layer && layer.hover) layer.hover(found.entry);
  }

  function hoverAt(event) {
    if (animation || event.pointerType === "touch") return;
    const found = picker.pick(event.clientX, event.clientY);
    hoverOther(found);
    if (found && found.layer !== choiceLayer) {
      if (hover && !hover.sticky) setHover(null);
      showTooltip(found.layer.tooltip ? found.layer.tooltip(found.entry) : "", event.clientX, event.clientY);
      return;
    }
    const index = found ? found.entry : -1;
    if (index < 0) {
      if (hover && !hover.sticky) setHover(null);
      showTooltip("", 0, 0);
      return;
    }
    if (!hover || hover.option !== index) setHover({ option: index });
    showTooltip(optionText(index), event.clientX, event.clientY);
  }

  function clickAt(event, type) {
    if (animation) return;
    // Inside a container, a phenomenon under the click is shown rather
    // than the block picked. Stars never take the click (MAP.101): in a
    // dense sector they would hide it.
    const found = picker.pick(event.clientX, event.clientY);
    if (found && found.layer !== choiceLayer) {
      showTooltip("", 0, 0);
      if (found.layer.select) found.layer.select(found.entry);
      return;
    }
    if (host.clearSelection) host.clearSelection();
    const index = found ? found.entry : -1;
    if (index < 0) return;
    if (type === "touch" && touchPending !== index) {
      // First tap highlights, a second on the same thing acts.
      touchPending = index;
      setHover({ option: index, sticky: true });
      showTooltip(optionText(index), event.clientX, event.clientY);
      return;
    }
    touchPending = -1;
    showTooltip("", 0, 0);
    act(index);
  }

  // Drag turns the view, right-drag or Shift-drag moves it, two fingers
  // zoom it, where the view is free; a press that hardly moves picks.
  const pointerControl = MC.createPointerControl(canvasEl, {
    dragClickPx: DRAG_CLICK_PX,
    buttons: [0, 2],
    isPan: function (event) { return event.button === 2 || event.shiftKey; },
    canDrag: function () { return isFree(resolved) && view && !animation; },
    onDragStart: function () { showTooltip("", 0, 0); },
    onDrag: function (dx, dy, pan) { if (view && !animation) drag(dx, dy, pan); },
    pinch: {
      canStart: function () { return MC.canZoom(zoomPolicyFor(resolved)) && !!view; },
      start: function () { pinchZoom = view.zoom; },
      move: function (ratio) {
        if (!view || animation) return;
        view.zoom = MC.clampDistance(zoomPolicyFor(resolved), pinchZoom * ratio, 1);
        reframe();
        applyView();
      },
    },
    onClick: clickAt,
    onHover: hoverAt,
    onLeave: onPointerLeave,
  });

  // Right-drag moves the view, so the map has no context menu where the
  // view is free.
  canvasEl.addEventListener("contextmenu", function (event) {
    if (isFree(resolved)) event.preventDefault();
  });

  function onPointerLeave() {
    showTooltip("", 0, 0);
    hoverOther(null);
    if (hover && !hover.sticky) setHover(null);
  }

  // Arrow keys move among the choices: Left and Right through them in
  // order, Up and Down to the next layer up or down, or (among the
  // segments of a slab) to the nearest one a ring further out or in;
  // Enter takes it.
  function onKey(event) {
    const key = event.key;
    if (key === "Escape" || key === "Backspace") {
      event.preventDefault();
      up();
      return;
    }
    if (key === "Home") {
      event.preventDefault();
      home();
      return;
    }
    if (["ArrowLeft", "ArrowRight", "ArrowUp", "ArrowDown", "Enter"].indexOf(key) < 0) return;
    event.preventDefault();
    // Shift and an arrow turn the view about the middle of what it shows
    // (MAP.53).
    if (event.shiftKey && key !== "Enter") {
      if (view && !animation && isFree(resolved)) {
        const turns = { ArrowLeft: [-KEY_TURN, 0], ArrowRight: [KEY_TURN, 0], ArrowUp: [0, KEY_TURN], ArrowDown: [0, -KEY_TURN] };
        turnView(turns[key][0], turns[key][1]);
        reframe();
        applyView();
      }
      return;
    }
    if (animation || !display || !display.options.length) return;
    const options = display.options;
    const current = hover && hover.option != null ? hover.option : -1;
    if (key === "Enter") {
      if (current >= 0) act(current);
      else if (hover && hover.layer) go({ at: stage.at, picks: stage.picks.concat([hover.layer]) });
      return;
    }
    let next = -1;
    if (current < 0) {
      next = 0;
    } else if (resolved.kind !== "segment" || key === "ArrowLeft" || key === "ArrowRight") {
      const forward = key === "ArrowRight" || key === "ArrowUp";
      next = (current + (forward ? 1 : -1) + options.length) % options.length;
    } else {
      next = segmentBeside(options, current, key === "ArrowUp" ? 1 : -1);
    }
    if (next < 0) return;
    setHover({ option: next, sticky: true });
    tooltipAtOption(next);
  }

  // The segment a ring further out (step 1) or in (-1) from segment
  // `current` whose middle bearing is nearest its own, or -1.
  function segmentBeside(options, current, step) {
    const here = options[current].blocks[0];
    const mid = function (b) { return (b.bounds.t0 + b.bounds.t1) / 2; };
    let best = -1;
    let bestTurn = Infinity;
    options.forEach(function (option, n) {
      const block = option.blocks[0];
      if (block.ring !== here.ring + step) return;
      const turn = Math.abs(wrapAngle(mid(block) - mid(here)));
      if (turn < bestTurn) {
        bestTurn = turn;
        best = n;
      }
    });
    return best;
  }

  // --- Breadcrumb, slab buttons, notice --------------------------------------

  function notice(text) {
    if (!els.notice) return;
    els.notice.textContent = text || "";
    els.notice.hidden = !text;
  }

  // The breadcrumb's ☆ (MAP.23, design doc section 8.2): saves the
  // selected sector, else the stage's own URL, in bookmarks.js.
  let bookmarkButton = null;
  let refreshBookmark = null;

  function bookmarkEntry() {
    if (!resolved) return null;
    if (selectedSector) {
      const designation = S.sectorDesignation(selectedSector.ring, selectedSector.layer, selectedSector.slot);
      const data = dataFor(stage.at);
      const sector = data ? data.sectors.get(selectedSector.ring + "/" + selectedSector.slot + "/" + selectedSector.layer) : null;
      // The sector page itself, without a pick mode's query.
      const page = sector && host.sectorUrl(sector.id) ? host.sectorUrl(sector.id).split("?")[0] : null;
      return {
        kind: "sector", value: designation,
        name: sector ? sector.name : "Sector " + S.blockLabel({ m: 1, ring: selectedSector.ring, wedge: selectedSector.slot, slab: selectedSector.layer }),
        url: page || location.pathname + "?sector=" + encodeURIComponent(designation),
        sectorId: sector ? sector.id : null,
      };
    }
    const labels = S.crumbs(stage, getOutline(), edgePc).map(function (crumb) { return crumb.label; });
    return { kind: "stage", value: location.pathname + S.stageQuery(stage), name: labels.join(" › ") };
  }

  // The breadcrumb (MAP.93): one line at any width. When the steps don't
  // all fit, it shows the first, a "…" button (a menu of the steps it
  // hides) and as many of the last ones as fit, then the current one,
  // measured again whenever the line's width changes. At phone width the
  // line gives way to the round Steps button between Back and Forward
  // (MAP.94, els.steps), a menu of every step.
  let crumbItems = [];
  let crumbWidth = -1;
  // The line's "…" menu, while it has one.
  let moreMenu = null;

  function crumbSteps() {
    const items = S.crumbs(stage, getOutline(), edgePc);
    if (selectedSector && !(resolved && resolved.sector)) {
      items.push({ label: "Sector " + S.blockLabel({ m: 1, ring: selectedSector.ring, wedge: selectedSector.slot, slab: selectedSector.layer }), last: true, sector: true });
      items[items.length - 2].last = false;
    }
    return items;
  }

  // One step: the current one as text, the others as buttons going there.
  function crumbNode(crumb, className) {
    if (crumb.last) {
      const here = document.createElement("span");
      here.textContent = crumb.label;
      here.setAttribute("aria-current", "location");
      return here;
    }
    const button = document.createElement("button");
    button.type = "button";
    button.className = className;
    button.textContent = crumb.label;
    button.addEventListener("click", function () { go(crumb.stage); });
    return button;
  }

  // A <details> menu of `steps` (a summary showing `text`, labelled
  // `label`); a step taken closes it.
  function stepsMenu(steps, text, label, className) {
    const menu = document.createElement("details");
    menu.className = className;
    const summary = document.createElement("summary");
    summary.textContent = text;
    summary.setAttribute("aria-label", label);
    summary.title = label;
    menu.appendChild(summary);
    const list = document.createElement("ul");
    steps.forEach(function (crumb) {
      const item = document.createElement("li");
      const node = crumbNode(crumb, "");
      if (!crumb.last) node.addEventListener("click", function () { menu.open = false; });
      item.appendChild(node);
      list.appendChild(item);
    });
    menu.appendChild(list);
    return menu;
  }

  // The line's steps, those between the first and the last `keep` folded
  // into "…" (none folded while keep covers them all).
  function fillCrumbs(list, keep) {
    list.textContent = "";
    const items = crumbItems;
    const folded = keep >= items.length - 1 ? [] : items.slice(1, items.length - keep);
    const shown = folded.length ? [items[0], null].concat(items.slice(items.length - keep)) : items;
    moreMenu = null;
    shown.forEach(function (crumb) {
      const item = document.createElement("li");
      if (crumb) {
        item.appendChild(crumbNode(crumb, "galaxy-crumb"));
      } else {
        item.className = "galaxy-crumb-more";
        moreMenu = stepsMenu(folded, "…", folded.length === 1 ? "1 more step" : folded.length + " more steps", "galaxy-crumb-menu");
        item.appendChild(moreMenu);
      }
      list.appendChild(item);
    });
  }

  function layoutCrumbs() {
    const nav = els.crumbs;
    const list = nav && nav.querySelector("ol");
    if (!list) return;
    crumbWidth = nav.clientWidth || 0;
    fillCrumbs(list, crumbItems.length);
    // Hidden, or no layout (a phone, where the line gives way): nothing
    // to measure.
    if (!list.clientWidth) return;
    // Measured at full length: the current step is cut short only when
    // even the shortest line doesn't fit.
    list.classList.add("galaxy-crumbs-measuring");
    for (let keep = crumbItems.length - 2; keep >= 1 && list.scrollWidth > list.clientWidth + 1; keep--) {
      fillCrumbs(list, keep);
    }
    list.classList.remove("galaxy-crumbs-measuring");
  }

  function renderSteps() {
    const box = els.steps;
    if (!box) return;
    const panel = box.querySelector("[data-steps-panel]");
    if (!panel) return;
    panel.textContent = "";
    const list = document.createElement("ol");
    crumbItems.forEach(function (crumb) {
      const item = document.createElement("li");
      const node = crumbNode(crumb, "");
      if (!crumb.last) node.addEventListener("click", function () { box.open = false; });
      item.appendChild(node);
      list.appendChild(item);
    });
    panel.appendChild(list);
  }

  function renderCrumbs() {
    const nav = els.crumbs;
    if (!nav) return;
    crumbItems = crumbSteps();
    let list = nav.querySelector("ol");
    if (!list) {
      nav.textContent = "";
      list = document.createElement("ol");
      nav.appendChild(list);
    }
    if (!bookmarkButton) {
      bookmarkButton = document.createElement("button");
      bookmarkButton.type = "button";
      bookmarkButton.className = "galaxy-bookmark";
      refreshBookmark = B.toggleButton(bookmarkButton, bookmarkEntry, true);
    }
    nav.appendChild(bookmarkButton);
    layoutCrumbs();
    renderSteps();
    refreshBookmark();
  }

  // The line is measured again when its width changes (a resize, the
  // panel shown).
  if (els.crumbs && typeof ResizeObserver === "function") {
    new ResizeObserver(function () {
      if (els.crumbs.clientWidth !== crumbWidth) layoutCrumbs();
    }).observe(els.crumbs);
  }
  // A "…" or Steps menu closes on a click elsewhere or Escape.
  function openStepMenus() {
    return [moreMenu, els.steps].filter(function (menu) { return menu && menu.open; });
  }
  document.addEventListener("click", function (event) {
    const path = event.composedPath ? event.composedPath() : [];
    openStepMenus().forEach(function (menu) {
      if (path.indexOf(menu) < 0) menu.open = false;
    });
  });
  document.addEventListener("keydown", function (event) {
    if (event.key !== "Escape") return;
    openStepMenus().forEach(function (menu) {
      menu.open = false;
      const summary = menu.querySelector("summary");
      if (summary && menu.contains(event.target)) summary.focus();
    });
  });

  // The slab buttons and their leader lines (MAP.54, layout MAP.76):
  // while the next pick is a slab, one button per slab beside the map
  // (in one column below it on a phone), each with a line from it to its
  // slab on the map. The lines are drawn on an SVG over the map's row and
  // redrawn whenever the view moves (applyView) or the layout changes, so
  // they always point at their slabs; the buttons are always in slab
  // number order, running the same way as the slabs on screen, so the
  // lines don't cross (MAP.110). The
  // lines are faint, and the slab hovered on the map or with its button
  // (or the focused button) lights its own; clicking a button picks its
  // slab. Each line ends on its slab's outline at the point nearest the
  // button (MAP.98); a slab off the map ends its line at the map's edge
  // with an arrow. Each button reads on one line, "#4 Unknown" or "#6 ≈
  // 2.43% charted" (MAP.100). A column taller than the map splits into
  // two, one each side of it, then the buttons shrink to their numbers,
  // and when even that won't fit they give way to picking on the map
  // (MAP.99). Otherwise the box says which layers the view holds.
  const SVG_NS = "http://www.w3.org/2000/svg";
  // {rows: [{option, item, button, line, outline, side}], list, other,
  // side, note, svg, mode, phone, fitFor}, while there are buttons.
  let strip = null;
  let leaderSvg = null;

  function slabChoices() {
    if (!resolved || resolved.kind !== "layer" || !display || display.resolved !== resolved) return null;
    return display.options;
  }

  // The SVG the lines are drawn on, over the map's row (created once).
  function leaderLayer() {
    const row = els.slabs && els.slabs.parentElement;
    if (!row) return null;
    if (!leaderSvg) {
      leaderSvg = document.createElementNS(SVG_NS, "svg");
      leaderSvg.setAttribute("class", "galaxy-slab-leaders");
      leaderSvg.setAttribute("aria-hidden", "true");
      const defs = document.createElementNS(SVG_NS, "defs");
      const marker = document.createElementNS(SVG_NS, "marker");
      marker.setAttribute("id", "galaxy-slab-arrow");
      marker.setAttribute("viewBox", "0 0 10 10");
      marker.setAttribute("refX", "9");
      marker.setAttribute("refY", "5");
      marker.setAttribute("markerWidth", "7");
      marker.setAttribute("markerHeight", "7");
      marker.setAttribute("orient", "auto-start-reverse");
      const tip = document.createElementNS(SVG_NS, "path");
      tip.setAttribute("d", "M0,0 L10,5 L0,10 z");
      tip.setAttribute("fill", "currentColor");
      marker.appendChild(tip);
      defs.appendChild(marker);
      leaderSvg.appendChild(defs);
      row.appendChild(leaderSvg);
    }
    return leaderSvg;
  }

  function clearLeaders() {
    if (!leaderSvg) return;
    Array.from(leaderSvg.querySelectorAll("polyline")).forEach(function (line) { line.remove(); });
    leaderSvg.hidden = true;
  }

  function renderStrip() {
    const box = els.slabs;
    strip = null;
    clearLeaders();
    clearSideBox();
    if (!box || !resolved || !resolved.view) return;
    box.textContent = "";
    const heading = document.createElement("h3");
    heading.id = "galaxymap3d-slabs-heading";
    // A thin block's view is sectors, so its "slabs" are sector layers.
    const noun = isSectorView(resolved) ? "Layer" : S.slabNoun(stage.at);
    heading.textContent = noun + "s";
    box.appendChild(heading);
    const choices = slabChoices();
    box.classList.toggle("is-picking", !!(choices && choices.length));
    if (!choices || !choices.length) {
      const slabs = slabsIn(resolved.view.blocks);
      const note = document.createElement("p");
      note.className = "galaxy-slab-note";
      if (slabs.length) {
        const lo = slabs[0];
        const hi = slabs[slabs.length - 1];
        note.textContent = "Showing " + (lo === hi ? noun.toLowerCase() + " " + lo : noun.toLowerCase() + "s " + lo + " to " + hi)
          + (resolved.kind === "segment" ? "; pick " + (isSectorView(resolved) ? "a sector" : "a block") + " of it on the map." : ".");
      }
      box.appendChild(note);
      return;
    }
    const data = display.data;
    const svg = leaderLayer();
    const list = document.createElement("ol");
    list.className = "galaxy-slab-buttons";
    box.appendChild(list);
    // Said in place of the buttons when none fit (MAP.99).
    const note = document.createElement("p");
    note.className = "galaxy-slab-note galaxy-slab-none-note";
    note.textContent = "Too many " + noun.toLowerCase() + "s for buttons here; pick one on the map.";
    note.hidden = true;
    box.appendChild(note);
    const side = sideBox();
    let other = null;
    if (side) {
      other = document.createElement("ol");
      other.className = "galaxy-slab-buttons";
      side.appendChild(other);
    }
    const rows = choices.slice().reverse().map(function (option) {
      const pick = option.pick;
      const sum = sumOf(option.blocks, data);
      const takeable = !(generatedOnly && !(sum.generated > 0));
      const item = document.createElement("li");
      const button = document.createElement("button");
      button.type = "button";
      button.className = "galaxy-slab-button";
      button.dataset.slab = String(pick.lo);
      // One line (MAP.100): "#4" and how much is charted, the second part
      // left out when the buttons are small (MAP.99).
      const label = S.slabButtonLabel(pick, sum.generated, getOutline().shapeless ? 0 : sum.total);
      const number = S.slabNumber(pick);
      const name = document.createElement("span");
      name.className = "galaxy-slab-name";
      name.textContent = number;
      const count = document.createElement("span");
      count.className = "galaxy-slab-count";
      count.textContent = label.slice(number.length + 1);
      button.appendChild(name);
      button.appendChild(count);
      button.title = S.pickLabel(pick, stage.at, resolved.view) + ": " + label.slice(number.length + 1);
      button.setAttribute("aria-label", layerText(pick) + ", " + S.formatInt(sum.generated)
        + (getOutline().shapeless ? "" : " of " + S.formatInt(sum.total)) + " sectors generated"
        + (takeable ? "" : " (nothing generated here to pick)"));
      if (!takeable) button.setAttribute("aria-disabled", "true");
      const light = function () { if (!animation) setHover({ layer: pick, sticky: true }); };
      const unlight = function () {
        const inStrip = list.contains(document.activeElement) || (other && other.contains(document.activeElement));
        if (hover && hover.layer && hover.layer.lo === pick.lo && !inStrip) setHover(null);
      };
      button.addEventListener("pointerenter", light);
      button.addEventListener("focus", light);
      button.addEventListener("pointerleave", unlight);
      button.addEventListener("blur", function () { setTimeout(unlight, 0); });
      button.addEventListener("click", function () {
        if (!takeable || animation) return;
        go({ at: stage.at, picks: stage.picks.concat([pick]) });
      });
      item.appendChild(button);
      list.appendChild(item);
      let line = null;
      if (svg) {
        line = document.createElementNS(SVG_NS, "polyline");
        line.setAttribute("class", "galaxy-slab-leader");
        line.dataset.slab = String(pick.lo);
        svg.appendChild(line);
      }
      // The slab's outline on the map (as its hover draws it), where its
      // line ends (MAP.98).
      const mid = (option.a0 + option.a1) / 2;
      const span = spanOf(option.blocks, mid);
      const outline = edgePoints(S.outlineEdges(option.blocks, mid), span.z1, span);
      return { option: option, item: item, button: button, line: line, outline: outline, side: "right" };
    });
    strip = { rows: rows, list: list, other: other, side: side, note: note, svg: svg, mode: null, fitFor: null };
    fitStrip(true);
    applyHover();
  }

  // --- Fitting the buttons to the map (MAP.99) ---------------------------------

  // The box on the map's other side (left) that takes half the buttons
  // when one column is taller than the map; made once, beside els.slabs.
  let sideEl = null;

  function sideBox() {
    const main = els.slabs && els.slabs.parentElement;
    if (!main) return null;
    if (!sideEl) {
      sideEl = document.createElement("div");
      sideEl.className = "galaxy-slabs galaxy-slabs-side";
      sideEl.id = "galaxymap3d-slabs-side";
      sideEl.setAttribute("role", "group");
      sideEl.setAttribute("aria-labelledby", "galaxymap3d-slabs-heading");
      sideEl.hidden = true;
      main.insertBefore(sideEl, main.firstChild);
    }
    sideEl.textContent = "";
    return sideEl;
  }

  function clearSideBox() {
    if (sideEl) {
      sideEl.textContent = "";
      sideEl.hidden = true;
    }
    const main = els.slabs && els.slabs.parentElement;
    if (main) main.classList.remove("is-slab-split", "is-slab-small", "is-slab-none");
  }

  function isPhone() {
    return typeof window.matchMedia === "function" && window.matchMedia("(max-width: 599px)").matches;
  }

  // Puts the strip in `mode`: "one" (one column beside the map, below it
  // on a phone), "split" (two columns, one each side), "small" (the
  // number alone: two columns beside the map, one below it on a phone) or
  // "none" (no buttons; the box says to pick on the map).
  function setStripMode(mode, phone) {
    const main = els.slabs.parentElement;
    const two = !phone && (mode === "split" || mode === "small");
    main.classList.toggle("is-slab-split", two);
    main.classList.toggle("is-slab-small", mode === "small");
    main.classList.toggle("is-slab-none", mode === "none");
    strip.list.hidden = mode === "none";
    strip.note.hidden = mode !== "none";
    if (strip.side) strip.side.hidden = !two;
    strip.mode = mode;
    strip.phone = phone;
    orderStrip(true);
  }

  // Whether the buttons fit: each column no taller than the map (below
  // the map on a phone, no taller than the window left under it).
  function stripFits(phone) {
    const rect = canvasEl.getBoundingClientRect();
    if (!rect.height) return true;
    const budget = phone ? Math.max(rect.height, window.innerHeight - rect.height - 48) : rect.height;
    const heading = els.slabs.querySelector("h3");
    const top = heading ? heading.getBoundingClientRect().height + 6 : 0;
    return [strip.list, strip.other].every(function (list) {
      if (!list || !list.children.length || list.parentElement.hidden) return true;
      const first = list.children[0].getBoundingClientRect();
      const last = list.children[list.children.length - 1].getBoundingClientRect();
      return top + last.bottom - first.top <= budget + 1;
    });
  }

  // Picks the first mode whose buttons fit, then draws the lines. Only
  // when the strip is new or the window changed size, so the map
  // resizing as the columns change doesn't refit it again.
  function fitStrip(force) {
    if (!strip) return;
    const phone = isPhone();
    const key = window.innerWidth + "x" + window.innerHeight;
    if (force || strip.fitFor !== key) {
      strip.fitFor = key;
      const modes = phone ? ["one", "small", "none"] : ["one", "split", "small", "none"];
      for (let n = 0; n < modes.length; n++) {
        setStripMode(modes[n], phone);
        if (modes[n] === "none" || stripFits(phone)) break;
      }
    }
    drawLeaders();
  }

  // Where a slab's line ends, in client pixels, when its outline gives
  // none: the middle of its blocks, moved toward the camera to the side of
  // the slab facing it, at the slab's mid height. {x, y, off}: off when
  // that point is not on the map (then x, y are on its edge). It also
  // says which way the buttons' number order runs (MAP.110).
  function slabAnchor(blocks, rect) {
    const fp = S.footprint(blocks);
    let z0 = Infinity;
    let z1 = -Infinity;
    blocks.forEach(function (block) {
      z0 = Math.min(z0, block.bounds.z0);
      z1 = Math.max(z1, block.bounds.z1);
    });
    const dx = camera.position.x - fp.center[0];
    const dy = camera.position.y - fp.center[1];
    const far = Math.hypot(dx, dy);
    const reach = far > 1e-9 ? Math.min(0.5 * fp.radius, far) / far : 0;
    const point = new THREE.Vector3(fp.center[0] + dx * reach, fp.center[1] + dy * reach, (z0 + z1) / 2);
    const ahead = point.clone().applyMatrix4(camera.matrixWorldInverse).z < 0;
    point.project(camera);
    let x = rect.left + ((point.x + 1) / 2) * rect.width;
    let y = rect.top + ((1 - point.y) / 2) * rect.height;
    if (!ahead) {
      // Behind the camera: toward the map's bottom edge.
      x = rect.left + rect.width / 2;
      y = rect.bottom;
    }
    const inset = 3;
    const cx = Math.max(rect.left + inset, Math.min(rect.right - inset, x));
    const cy = Math.max(rect.top + inset, Math.min(rect.bottom - inset, y));
    return { x: cx, y: cy, off: !ahead || cx !== x || cy !== y };
  }

  // The slab's outline on screen as line pieces [[x, y], [x, y]] in client
  // pixels, each cut to the map (pieces with an end behind the camera
  // left out).
  function outlineOnScreen(points, rect) {
    const inset = 3;
    const box = { x0: rect.left + inset, y0: rect.top + inset, x1: rect.right - inset, y1: rect.bottom - inset };
    const screen = points.map(function (p) {
      if (p.clone().applyMatrix4(camera.matrixWorldInverse).z >= 0) return null;
      const q = p.clone().project(camera);
      return [rect.left + ((q.x + 1) / 2) * rect.width, rect.top + ((1 - q.y) / 2) * rect.height];
    });
    const pieces = [];
    for (let n = 0; n + 1 < screen.length; n += 2) {
      const a = screen[n];
      const b = screen[n + 1];
      if (!a || !b) continue;
      const cut = clipPiece(a, b, box);
      if (cut) pieces.push(cut);
    }
    return pieces;
  }

  // Piece a-b cut to box {x0, y0, x1, y1} (Liang-Barsky), or null.
  function clipPiece(a, b, box) {
    const dx = b[0] - a[0];
    const dy = b[1] - a[1];
    let t0 = 0;
    let t1 = 1;
    const sides = [[-dx, a[0] - box.x0], [dx, box.x1 - a[0]], [-dy, a[1] - box.y0], [dy, box.y1 - a[1]]];
    for (let n = 0; n < sides.length; n++) {
      const p = sides[n][0];
      const q = sides[n][1];
      if (p === 0) {
        if (q < 0) return null;
        continue;
      }
      const t = q / p;
      if (p < 0) t0 = Math.max(t0, t);
      else t1 = Math.min(t1, t);
      if (t0 > t1) return null;
    }
    return [[a[0] + dx * t0, a[1] + dy * t0], [a[0] + dx * t1, a[1] + dy * t1]];
  }

  // Where a slab's line ends (MAP.98): on its outline as drawn on the map,
  // at the point nearest `from` (the line's last bend), or with `lane` (an
  // x), nearest that upright lane (a phone's lines come across from it).
  // With `minY`, no higher than that where the outline allows (MAP.110).
  // Falls back to slabAnchor when none of the outline is on the map.
  function slabEnd(r, from, lane, rect, minY) {
    const best = S.nearestOnPieces(outlineOnScreen(r.outline, rect), from, lane, minY);
    if (!best) return slabAnchor(r.option.blocks, rect);
    return { x: best.x, y: best.y, off: false };
  }

  // Orders the buttons by slab number (MAP.110): the highest first when
  // the stack's top slab shows above its bottom one on screen, else the
  // lowest first, so the lines run the same way as the slabs (on a phone
  // the other way round, so the lanes nest). With two columns deals them
  // out in turn, right first, so each column runs the map's whole height.
  // Leaves them be while one has focus, unless `force` (the columns
  // changed).
  function orderStrip(force) {
    if (!strip || !view) return;
    const rect = canvasEl.getBoundingClientRect();
    if (!rect.width || !rect.height) return;
    const byNumber = S.slabOrder(strip.rows, true, function (r) { return r.option.pick; });
    const top = slabAnchor(byNumber[0].option.blocks, rect).y;
    const bottom = slabAnchor(byNumber[byNumber.length - 1].option.blocks, rect).y;
    const order = (top <= bottom) !== !!strip.phone ? byNumber : byNumber.reverse();
    const two = strip.side && !strip.side.hidden;
    const focused = strip.list.contains(document.activeElement) || (strip.other && strip.other.contains(document.activeElement));
    if (focused && !force) return;
    order.forEach(function (r, n) {
      r.side = two && n % 2 === 1 ? "left" : "right";
      const list = r.side === "left" ? strip.other : strip.list;
      if (list.children[list.children.length - 1] !== r.item) list.appendChild(r.item);
    });
  }

  // Redraws the lines from each button to its slab, reordering the
  // buttons first if their slabs' order on screen changed.
  function drawLeaders() {
    if (!strip || !strip.svg || !view || !active) return;
    const svg = strip.svg;
    const row = svg.parentElement;
    const base = row.getBoundingClientRect();
    const rect = canvasEl.getBoundingClientRect();
    if (!rect.width || !rect.height || !base.width || strip.mode === "none") {
      svg.hidden = true;
      return;
    }
    svg.hidden = false;
    svg.setAttribute("width", String(base.width));
    svg.setAttribute("height", String(base.height));
    svg.setAttribute("viewBox", "0 0 " + base.width + " " + base.height);
    orderStrip();
    const rows = strip.rows;
    // Below the map (a phone) the lines run up lanes along the map's right
    // edge, so they nest: the bottom button's line takes the outer lane
    // to the highest slab, and the column runs from the lowest slab down.
    const below = rows[0].button.getBoundingClientRect().top >= rect.bottom - 1;
    const items = Array.from(strip.list.children);
    // Beside the map each line ends where its slab's outline comes
    // nearest its button, but down each column no higher than the line
    // above it ends, where the outline allows, so the lines keep the
    // buttons' number order and don't cross (MAP.110).
    const ends = below ? null : new Map();
    if (ends) {
      [strip.list, strip.other].forEach(function (list) {
        if (!list) return;
        let floor = -Infinity;
        Array.from(list.children).forEach(function (item) {
          const r = rows.find(function (row) { return row.item === item; });
          const b = r.button.getBoundingClientRect();
          const from = [r.side === "left" ? b.right : b.left, b.top + b.height / 2];
          const end = slabEnd(r, from, null, rect, floor);
          ends.set(r, { from: from, end: end });
          floor = Math.max(floor, end.y);
        });
      });
    }
    rows.forEach(function (r) {
      if (!r.line) return;
      const b = r.button.getBoundingClientRect();
      const points = [];
      let end;
      if (below) {
        // Out of the button's right end to its own lane, up the lane to
        // the slab's height, then across to the slab.
        const lane = rect.right - 8 - 7 * (items.length - 1 - items.indexOf(r.item));
        end = slabEnd(r, [lane, b.top + b.height / 2], lane, rect);
        points.push([b.right, b.top + b.height / 2], [lane, b.top + b.height / 2], [lane, end.y]);
      } else {
        end = ends.get(r).end;
        points.push(ends.get(r).from);
      }
      points.push([end.x, end.y]);
      r.line.setAttribute("points", points.map(function (p) {
        return (p[0] - base.left).toFixed(1) + "," + (p[1] - base.top).toFixed(1);
      }).join(" "));
      if (end.off) r.line.setAttribute("marker-end", "url(#galaxy-slab-arrow)");
      else r.line.removeAttribute("marker-end");
      r.line.classList.toggle("is-off", end.off);
    });
  }

  // Lights the button and line of the slab holding `layer` ({lo, hi}, or
  // null for none).
  function markStripRow(layer) {
    if (!strip) return;
    strip.rows.forEach(function (r) {
      const lit = !!layer && r.option.pick.lo <= layer.lo && r.option.pick.hi >= layer.hi;
      r.button.classList.toggle("is-lit", lit);
      if (r.line) r.line.classList.toggle("is-lit", lit);
    });
  }

  // A resized map (the window resized or turned) refits the view to its
  // new size (MAP.53), and the lines follow the layout: the map or the
  // buttons resized, or the page reflowed.
  function relayout() {
    if (active && view && !animation) {
      reframe();
      applyView();
    }
    fitStrip(false);
  }
  if (typeof ResizeObserver === "function") {
    const observer = new ResizeObserver(relayout);
    observer.observe(canvasEl);
    if (els.slabs) observer.observe(els.slabs);
  }
  window.addEventListener("resize", relayout);

  // --- The address bar (section 9.3) -----------------------------------------

  // Flies to the stage showing sector {ring, layer, slot} among its
  // layer's neighbours (MAP.26) and selects it.
  function locate(sector) {
    const next = S.sectorStage(sector.ring, sector.layer, sector.slot, getOutline(), edgePc);
    if (!next) {
      notice("Sector " + S.blockLabel({ m: 1, ring: sector.ring, wedge: sector.slot, slab: sector.layer })
        + " is outside the galaxy.");
      return false;
    }
    selectedSector = { ring: sector.ring, layer: sector.layer, slot: sector.slot };
    go(next, { keepSector: true, query: "?sector=" + S.sectorDesignation(sector.ring, sector.layer, sector.slot) });
    return true;
  }

  function clearMatches() {
    if (els.matches) {
      els.matches.textContent = "";
      els.matches.hidden = true;
    }
  }

  // The name lookup's matches, as buttons that fly to each one's sector.
  function showMatches(matches) {
    const box = els.matches;
    if (!box) return;
    box.textContent = "";
    const list = document.createElement("ul");
    matches.forEach(function (match) {
      const item = document.createElement("li");
      const button = document.createElement("button");
      button.type = "button";
      button.textContent = match.kind === "system"
        ? match.name + " (system in " + (match.sector_name || "sector " + match.sector_id) + ")"
        : match.name + " (sector)";
      button.addEventListener("click", function () {
        clearMatches();
        locate(match);
      });
      item.appendChild(button);
      list.appendChild(item);
    });
    box.appendChild(list);
    box.hidden = false;
  }

  function onAddress(event) {
    event.preventDefault();
    const input = els.address.querySelector("input");
    clearMatches();
    const asked = S.parseAddress(input.value, edgePc);
    if (asked.problem) {
      notice(asked.problem);
      return;
    }
    if (asked.sector) {
      locate(asked.sector);
      return;
    }
    notice("Looking up " + asked.name + "…");
    host.locate(asked.name).then(function (matches) {
      if (!matches.length) {
        notice("Nothing is named like " + asked.name + ".");
        return;
      }
      const exact = matches.filter(function (m) { return (m.name || "").toLowerCase() === asked.name.toLowerCase(); });
      const sectors = new Set((exact.length ? exact : matches).map(function (m) { return m.sector_id; }));
      if (matches.length === 1 || (exact.length && sectors.size === 1)) {
        notice("");
        locate(exact[0] || matches[0]);
        return;
      }
      notice(matches.length + " names match; pick one.");
      showMatches(matches);
    }, function () {
      notice("The lookup failed. Please try again shortly.");
    });
  }

  if (els.address) els.address.addEventListener("submit", onAddress);

  // --- Turning on, and the URL -----------------------------------------------

  function setActive(on) {
    if (on === active) return;
    active = on;
    [els.crumbs, els.slabs, els.address].forEach(function (el) { if (el) el.hidden = !on; });
    if (!on) {
      clearMatches();
      if (animation) finishAnimation();
      disposeDisplay(display);
      display = null;
      view = null;
      clearOutline();
      strip = null;
      clearLeaders();
      clearSideBox();
      showTooltip("", 0, 0);
      notice("");
      return;
    }
    const r = resolved || resolve(stage);
    loadData(r.stage.at).then(function () {
      if (!active) return;
      display = null;
      flyTo(r);
    });
  }

  // The stage a URL's query asks for: {stage, sector, problem}. ?sector=
  // opens the stage showing that sector; with a NAV course and no stage
  // asked for, the smallest stage showing the whole course (section 9.4).
  function stageFromLocation(useCourse) {
    const parsed = S.parseStageQuery(location.search);
    if (parsed.problem) return { stage: { at: null, picks: [] }, sector: null, problem: parsed.problem };
    if (parsed.sector) {
      const s = parsed.sector;
      const next = S.sectorStage(s.ring, s.layer, s.slot, getOutline(), edgePc);
      if (!next) {
        return {
          stage: { at: null, picks: [] }, sector: null,
          problem: "Sector " + S.blockLabel({ m: 1, ring: s.ring, wedge: s.slot, slab: s.layer }) + " is outside the galaxy.",
        };
      }
      return { stage: next, sector: s, problem: null };
    }
    let next = parsed.stage;
    if (useCourse && !next.at && !next.picks.length && host.courseSectors && host.courseSectors.length) {
      next = S.courseStage(host.courseSectors, getOutline(), edgePc);
    }
    return { stage: next, sector: null, problem: null };
  }

  // Opens the stage the URL asks for (or the galaxy, with a notice).
  function openFromLocation() {
    const state = history.state && history.state.galaxyStage ? history.state : null;
    mapIndex = state && state.mapIndex != null ? state.mapIndex : 0;
    maxIndex = Math.max(mapIndex, state && state.maxIndex != null ? state.maxIndex : mapIndex);
    history.replaceState({ galaxyStage: true, mapIndex: mapIndex, maxIndex: maxIndex }, "");
    const asked = stageFromLocation(true);
    let r = resolve(asked.stage);
    let problem = asked.problem || r.problem;
    if (r.problem) r = resolve({ at: null, picks: [] });
    selectedSector = asked.sector || (r.sector ? { ring: r.sector.ring, layer: r.sector.slab, slot: r.sector.wedge } : null);
    stage = r.stage;
    resolved = r;
    renderCrumbs();
    renderStrip();
    renderTravel();
    notice(problem || "");
  }

  function onPopState(event) {
    if (!active) return;
    const state = event.state && event.state.galaxyStage ? event.state : null;
    mapIndex = state && state.mapIndex != null ? state.mapIndex : 0;
    maxIndex = Math.max(maxIndex, mapIndex);
    // Remembered on this entry too, so a reload here still knows how far
    // Forward goes.
    if (state) history.replaceState({ galaxyStage: true, mapIndex: mapIndex, maxIndex: maxIndex }, "");
    const asked = stageFromLocation(false);
    if (asked.problem) {
      selectedSector = null;
      go({ at: null, picks: [] }, { push: false });
      notice(asked.problem);
      return;
    }
    selectedSector = asked.sector;
    go(asked.stage, { push: false, keepSector: true });
  }
  window.addEventListener("popstate", onPopState);

  function setGeneratedOnly(on) {
    generatedOnly = on;
    if (display && !animation) rebuildDisplay();
  }

  return {
    setActive: setActive,
    isActive: function () { return active; },
    openFromLocation: openFromLocation,
    step: step,
    invalidate: invalidate,
    onPointerDown: pointerControl.down,
    onPointerMove: pointerControl.move,
    onPointerUp: pointerControl.up,
    onPointerLeave: pointerControl.leave,
    onKey: onKey,
    onWheel: onWheel,
    resetView: resetView,
    home: home,
    up: up,
    travel: travel,
    setGeneratedOnly: setGeneratedOnly,
    stage: function () { return stage; },
    // Inside a block or an arc (not the whole galaxy), where phenomena
    // take the pointer before the blocks.
    inContainer: function () { return !!stage.at; },
    go: go,
    locate: locate,
  };
}
