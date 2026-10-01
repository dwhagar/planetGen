// html/static/galaxystageview.js
//
// The Galaxy Map's drill-down (docs/design/galaxy-drilldown-navigation.md,
// sections 4, 5 and 8.1), drawn and driven: eight stages from the galaxy
// to a sector, alternating a 3D view of a container's children (pick a
// slab) and a top-down view of one slab (pick a block). galaxystages.js
// has the rules; this file has the scene, the camera moves, the
// breadcrumb, the slab strip, the tooltip, keys, touch and the stage URLs.
// galaxymap3d.js creates it (createStageView) and hands it the pointer,
// wheel and key events while the map is in stage mode, and calls step()
// every frame.
//
// host (from galaxymap3d.js):
// - THREE, scene, camera, canvasEl;
// - edgePc, shape (or null), galaxyRadius, reducedMotion;
// - blockScene: galaxyblocks' scene (buildCells);
// - makeBlockMesh(part, translucent): a mesh with the map's block shader;
// - setCamera({target: [x, y, z], dist, theta, phi}): moves the camera
//   (and the tiles, scale bar and wedge lines with it);
// - fetchStage(query): a Promise of GET /galaxy/stage's JSON for a query
//   ("?at=m.ring.wedge.slab", or "" for the galaxy);
// - setBlockSize(m): the scale readout's "1 block = m sectors";
// - showBlockInfo(info), showPlacedInfo(entry), showCellInfo(cell),
//   showHint(text): the info panel (showBlockInfo's info: {block, total,
//   generated, hint, enter} -- enter() flies into it, null for the block
//   the stage shows);
// - canGenerate: whether the visitor gets Generate buttons;
// - sectorUrl(id): a generated sector's page;
// - locate(name): a Promise of GET /galaxy/locate's matches for a name
//   (the address bar);
// - els: {crumbs, slabs, tooltip, notice, address, matches} (any may be
//   missing).

const VERSION_QUERY = new URL(import.meta.url).search;

const S = await import(`./galaxystages.js${VERSION_QUERY}`);

// The 3D stages' tilt from straight down, and the top-down views' (not
// quite 0: the camera keeps galactic north as its up vector, which needs
// the view direction off vertical by a hair).
const TOP_DOWN_PHI = 1e-3;
// Other slabs, while one is hovered.
const OTHER_SLAB_FADE = 0.35;
// The 3D stages' rotation limits and drag speed.
const MIN_TILT = (8 * Math.PI) / 180;
const MAX_TILT = (82 * Math.PI) / 180;
const ROTATE_PER_PX = (0.4 * Math.PI) / 180;
// Wheel zoom inside a stage: 3D views from fitting the container to 1.5
// times closer; stage 2 (the galaxy's slab) up to 4 times.
const ZOOM_3D = 1.5;
const ZOOM_GALAXY_SLAB = 4;
const WHEEL_ZOOM_PER_PX = 0.0025;
// A new stage's blocks fade in over the last part of a flight.
const FADE_IN_MS = 200;
const DRAG_CLICK_PX = 4;

export function createStageView(host) {
  const THREE = host.THREE;
  const camera = host.camera;
  const canvasEl = host.canvasEl;
  const edgePc = host.edgePc;
  const els = host.els || {};
  const accent = new THREE.Color(host.accentColor || "#4f5fe8");

  let active = false;
  let outline = null;
  const childCache = new Map();
  const dataCache = new Map();
  let stage = { at: null, slab: null };
  let selectedSector = null;
  let display = null;
  let leaving = [];
  let animation = null;
  let hover = null;
  let generatedOnly = false;
  // The camera as the stage left it, and the fit it zooms within.
  let view = null;
  // Bumped by every go(), so a stage whose data arrives after a newer
  // go() is dropped.
  let goToken = 0;

  // --- The outline and the children --------------------------------------

  function getOutline() {
    if (!outline) {
      outline = S.galaxyOutline(edgePc, host.shape, host.galaxyRadius);
      outline.shapeless = !host.shape;
    }
    return outline;
  }

  function containerKey(at) {
    return at ? [at.m, at.ring, at.wedge, at.slab].join(".") : "galaxy";
  }

  function childrenOf(at) {
    const key = containerKey(at);
    let groups = childCache.get(key);
    if (!groups) {
      groups = S.stageChildren(at, getOutline(), edgePc);
      childCache.set(key, groups);
    }
    return groups;
  }

  // The stage API's answer for container `at`, as {generated: Map(key ->
  // count), sectors: Map("ring/slot/layer" -> sector)}; null while it
  // loads (and on an error, which leaves counts at 0).
  function dataFor(at) {
    const key = containerKey(at);
    const cached = dataCache.get(key);
    if (cached && cached.ready) return cached.value;
    return null;
  }

  function loadData(at) {
    const key = containerKey(at);
    let entry = dataCache.get(key);
    if (!entry) {
      entry = { ready: false, value: null };
      entry.promise = host.fetchStage(at ? S.stageQuery({ at: at, slab: null }) : "").then(function (payload) {
        const generated = new Map();
        (payload.children || []).forEach(function (child) {
          generated.set(child.ring + "/" + child.wedge + "/" + child.slab, child.generated);
        });
        const sectors = new Map();
        (payload.sectors || []).forEach(function (sector) {
          sectors.set(sector.ring + "/" + sector.slot + "/" + sector.layer, sector);
        });
        entry.value = { generated: generated, sectors: sectors };
        entry.ready = true;
        return entry.value;
      }, function () {
        entry.value = { generated: new Map(), sectors: new Map(), failed: true };
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
    return data.generated.get(block.ring + "/" + block.wedge + "/" + block.slab) || 0;
  }

  // The groups a stage shows: every slab in 3D, the one slab top-down.
  // Without an outline (no shape yet), only blocks holding generated
  // sectors.
  function stageGroups(s, data) {
    let groups = childrenOf(s.at);
    if (s.slab != null) groups = groups.filter(function (g) { return g.slab === s.slab; });
    if (getOutline().shapeless) {
      groups = groups.map(function (g) {
        return { slab: g.slab, blocks: g.blocks.filter(function (b) { return generatedOf(b, data) > 0; }) };
      }).filter(function (g) { return g.blocks.length; });
    }
    return groups;
  }

  // --- Drawing a stage ---------------------------------------------------------

  // One stage's meshes: a group per slab (so a slab can rise and the
  // others fade), with the blocks behind each mesh for picking.
  function buildDisplay(s) {
    const data = dataFor(s.at);
    const groups = stageGroups(s, data);
    const root = new THREE.Group();
    const slabs = [];
    const eye = camera.position.toArray();
    groups.forEach(function (g) {
      const cells = g.blocks.map(function (block) {
        const b = block.bounds;
        return {
          ring: block.ring, seg: block.wedge, slab: block.slab,
          r0: b.r0, r1: b.r1, t0: b.t0, t1: b.t1, z0: b.z0, z1: b.z1,
          filled: generatedOf(block, data), total: block.total, block: block,
        };
      });
      const dim = generatedOnly ? function (cell) { return !(cell.filled > 0); } : null;
      const built = host.blockScene.buildCells(cells, eye, dim);
      const group = new THREE.Group();
      const meshes = [];
      [built.solid, built.glass].forEach(function (part, n) {
        if (!part.vertexCount) return;
        const mesh = host.makeBlockMesh(part, n === 1);
        mesh.userData.cells = built.cells[n];
        mesh.renderOrder = n;
        group.add(mesh);
        meshes.push(mesh);
      });
      root.add(group);
      slabs.push({ slab: g.slab, blocks: g.blocks, group: group, meshes: meshes, fade: 1 });
    });
    host.scene.add(root);
    return { stage: s, root: root, slabs: slabs, fade: 1, data: data };
  }

  function disposeDisplay(d) {
    if (!d) return;
    host.scene.remove(d.root);
    d.slabs.forEach(function (slab) {
      slab.meshes.forEach(function (mesh) {
        mesh.geometry.dispose();
        mesh.material.dispose();
      });
    });
  }

  function setDisplayFade(d, fade) {
    d.fade = fade;
    d.slabs.forEach(function (slab) {
      slab.meshes.forEach(function (mesh, n) {
        const value = fade * slab.fade;
        mesh.material.uniforms.fade.value = value;
        mesh.material.transparent = n === 1 || mesh.userData.cells === undefined || value < 1 || mesh.renderOrder === 1;
        mesh.material.depthWrite = mesh.renderOrder === 0 && value >= 1;
      });
    });
  }

  function rebuildDisplay() {
    const old = display;
    display = buildDisplay(stage);
    if (old) {
      display.slabs.forEach(function (slab) {
        const before = old.slabs.find(function (o) { return o.slab === slab.slab; });
        if (before) slab.fade = before.fade;
      });
    }
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

  function zRange(groups) {
    let z0 = Infinity;
    let z1 = -Infinity;
    groups.forEach(function (g) {
      g.blocks.forEach(function (block) {
        z0 = Math.min(z0, block.bounds.z0);
        z1 = Math.max(z1, block.bounds.z1);
      });
    });
    return [z0, z1];
  }

  function containerTheta(at) {
    if (!at) return null;
    const b = S.drillBlockBounds(at, edgePc);
    return (b.t0 + b.t1) / 2 + Math.PI;
  }

  // The camera a stage starts from: {target, dist, theta, phi}.
  function cameraFor(s, groups) {
    groups = groups || stageGroups(s, dataFor(s.at));
    const blocks = [];
    groups.forEach(function (g) { Array.prototype.push.apply(blocks, g.blocks); });
    if (!blocks.length) {
      return { target: [0, 0, 0], dist: host.galaxyRadius * 2.4, theta: -Math.PI * 32 / 180, phi: (S.STAGE_TILT_DEG * Math.PI) / 180 };
    }
    const fp = S.footprint(blocks);
    const z = zRange(groups);
    const half = fovHalf();
    if (s.slab != null) {
      const theta = s.at ? containerTheta(s.at) : -Math.PI / 2;
      return {
        target: [fp.center[0], fp.center[1], (z[0] + z[1]) / 2],
        dist: (S.FIT_MARGIN * fp.radius) / Math.tan(half) + (z[1] - z[0]) / 2,
        theta: theta, phi: TOP_DOWN_PHI,
      };
    }
    const radius = Math.hypot(fp.radius, (z[1] - z[0]) / 2);
    const theta = s.at ? containerTheta(s.at) : (view && !view.stage.at ? view.theta : (-32 * Math.PI) / 180);
    return {
      target: [fp.center[0], fp.center[1], (z[0] + z[1]) / 2],
      dist: (S.FIT_MARGIN * radius) / Math.sin(half),
      theta: theta, phi: (S.STAGE_TILT_DEG * Math.PI) / 180,
    };
  }

  function applyView() {
    host.setCamera({ target: view.target.slice(), dist: view.dist, theta: view.theta, phi: view.phi });
  }

  // --- Moving between stages -------------------------------------------------

  function childM(at) {
    return at ? at.m / (at.m === 3 ? 3 : 9) : 243;
  }

  function wrapAngle(a) {
    return Math.atan2(Math.sin(a), Math.cos(a));
  }

  // Animates the camera from `from` to `to` along van Wijk and Nuij's
  // path in the plane, easing the height, tilt and turn alongside; calls
  // progress(t) each frame (t 0..1) and done() at the end.
  function flyCamera(from, to, ms, progress, done) {
    const tanHalf = Math.tan(THREE.MathUtils.degToRad(camera.fov) / 2);
    const path = S.flightPath(from.target, 2 * from.dist * tanHalf, to.target, 2 * to.dist * tanHalf);
    const duration = ms != null ? ms : S.flightMs(path.S);
    const turn = wrapAngle(to.theta - from.theta);
    animation = {
      startedAt: performance.now(),
      duration: host.reducedMotion ? 0 : duration,
      frame: function (t) {
        const e = S.easeInOut(t);
        const p = path.at(path.S * t);
        view = {
          stage: view ? view.stage : stage,
          target: [p.center[0], p.center[1], from.target[2] + (to.target[2] - from.target[2]) * e],
          dist: p.w / (2 * tanHalf),
          theta: from.theta + turn * e,
          phi: from.phi + (to.phi - from.phi) * e,
          fit: to.fit,
        };
        applyView();
        if (progress) progress(t);
      },
      done: done,
    };
  }

  // Goes to stage `next`. A 3D stage's own slab is pulled out (it rises
  // while the others fade, the camera turns to look straight down); any
  // other move flies there. push: record it in the browser history.
  function go(next, options) {
    options = options || {};
    if (animation) finishAnimation();
    const problem = S.validStage(next, getOutline(), edgePc);
    if (problem) {
      notice(problem);
      return;
    }
    notice("");
    if (!options.keepSector) selectedSector = null;
    const previous = stage;
    const token = ++goToken;
    const query = options.query || S.stageQuery(next);
    if (options.push !== false && (!S.sameStage(previous, next) || query !== location.search)) {
      history.pushState({ galaxyStage: true }, "", location.pathname + query + location.hash);
    }
    loadData(next.at).then(function () {
      if (!active || token !== goToken) return;
      if (animation) finishAnimation();
      if (next.slab != null && S.sameBlock(previous.at, next.at) && previous.slab == null && display
          && S.sameStage(display.stage, previous)) {
        pullOut(next);
      } else {
        flyTo(next);
      }
    });
    stage = next;
    renderCrumbs();
    renderStrip();
    host.setBlockSize(childM(next.at));
  }

  function pullOut(next) {
    const slabInfo = display.slabs.find(function (s) { return s.slab === next.slab; });
    const from = currentView();
    const to = cameraFor(next);
    const thickness = slabInfo ? slabInfo.blocks[0].bounds.z1 - slabInfo.blocks[0].bounds.z0 : 0;
    const old = display;
    to.fit = to.dist;
    flyCamera(from, to, S.PULL_OUT_MS, function (t) {
      const e = S.easeInOut(t);
      old.slabs.forEach(function (slab) {
        if (slab.slab === next.slab) {
          slab.group.position.z = 0.6 * thickness * e;
          slab.fade = 1;
        } else {
          slab.fade = 1 - e;
        }
      });
      setDisplayFade(old, 1);
    }, function () {
      view.stage = next;
      display = null;
      disposeDisplay(old);
      display = buildDisplay(next);
      setDisplayFade(display, 1);
      view = Object.assign({}, to, { stage: next, fit: to.dist });
      applyView();
      afterArrival();
    });
  }

  function flyTo(next) {
    const from = display ? currentView() : null;
    const to = cameraFor(next);
    to.fit = to.dist;
    const old = display;
    const incoming = buildDisplay(next);
    setDisplayFade(incoming, from ? 0 : 1);
    if (!from) {
      display = incoming;
      view = Object.assign({}, to, { stage: next });
      applyView();
      afterArrival();
      return;
    }
    flyCamera(from, to, null, function (t) {
      const ms = animation ? animation.duration : 0;
      const fadeStart = ms > 0 ? Math.max(0, 1 - FADE_IN_MS / ms) : 0;
      const f = fadeStart >= 1 ? 1 : Math.max(0, (t - fadeStart) / (1 - fadeStart));
      setDisplayFade(incoming, f);
      if (old) setDisplayFade(old, 1 - f);
    }, function () {
      display = incoming;
      setDisplayFade(display, 1);
      view = Object.assign({}, to, { stage: next });
      applyView();
      afterArrival();
    });
    display = incoming;
    if (old) leaving.push(old);
  }

  function afterArrival() {
    leaving.forEach(function (d) { if (d && d !== display) disposeDisplay(d); });
    leaving = [];
    hover = null;
    applyHover();
    renderStrip();
    if (selectedSector) {
      const block = sectorBlock(selectedSector);
      if (block) setHover({ block: block, sticky: true });
    } else {
      showStageInfo();
    }
  }

  function finishAnimation() {
    const a = animation;
    animation = null;
    a.frame(1);
    a.done();
  }

  function currentView() {
    if (view) return { target: view.target.slice(), dist: view.dist, theta: view.theta, phi: view.phi };
    return cameraFor(stage);
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
    const parent = S.parentStage(stage);
    if (parent) go(parent);
  }

  // --- Picking and hover -----------------------------------------------------

  const raycaster = new THREE.Raycaster();

  // The block under a screen point: {block, slab} or null. Dimmed blocks
  // ("Generated only") can't be picked.
  function blockAt(clientX, clientY) {
    if (!display) return null;
    const rect = canvasEl.getBoundingClientRect();
    if (!rect.width || !rect.height) return null;
    const ndc = new THREE.Vector2(((clientX - rect.left) / rect.width) * 2 - 1, -((clientY - rect.top) / rect.height) * 2 + 1);
    raycaster.setFromCamera(ndc, camera);
    const meshes = [];
    display.slabs.forEach(function (slab) {
      if (slab.fade > 0.05) Array.prototype.push.apply(meshes, slab.meshes);
    });
    const hits = raycaster.intersectObjects(meshes, false);
    for (let n = 0; n < hits.length; n++) {
      const hit = hits[n];
      const owners = hit.object.userData.part.owners;
      const cell = hit.object.userData.cells[owners[hit.face.a]];
      if (!cell) continue;
      if (generatedOnly && !(cell.filled > 0)) continue;
      return { block: cell.block, slab: cell.slab, cell: cell };
    }
    return null;
  }

  let outlineLine = null;

  function clearOutline() {
    if (outlineLine) {
      host.scene.remove(outlineLine);
      outlineLine.geometry.dispose();
      outlineLine.material.dispose();
      outlineLine = null;
    }
  }

  // An accent outline around a block's top face (top-down stages).
  function outlineBlock(block) {
    clearOutline();
    const b = block.bounds;
    const steps = Math.max(2, Math.ceil((b.t1 - b.t0) / (Math.PI / 90)));
    const points = [];
    const z = b.z1;
    for (let k = 0; k <= steps; k++) {
      const t = b.t0 + ((b.t1 - b.t0) * k) / steps;
      points.push(new THREE.Vector3(b.r1 * Math.cos(t), b.r1 * Math.sin(t), z));
    }
    for (let k = steps; k >= 0; k--) {
      const t = b.t0 + ((b.t1 - b.t0) * k) / steps;
      points.push(new THREE.Vector3(b.r0 * Math.cos(t), b.r0 * Math.sin(t), z));
    }
    outlineLine = new THREE.LineLoop(
      new THREE.BufferGeometry().setFromPoints(points),
      new THREE.LineBasicMaterial({ color: accent, depthTest: false, transparent: true }),
    );
    outlineLine.renderOrder = 6;
    host.scene.add(outlineLine);
  }

  function setHover(next) {
    hover = next;
    applyHover();
  }

  function applyHover() {
    if (!display) return;
    const topDown = stage.slab != null;
    display.slabs.forEach(function (slab) {
      if (!topDown && !animation) {
        slab.fade = hover && hover.slab != null && hover.slab !== slab.slab ? OTHER_SLAB_FADE : 1;
      }
    });
    if (!animation) setDisplayFade(display, 1);
    if (topDown && hover && hover.block) outlineBlock(hover.block);
    else clearOutline();
    markStripRow(hover ? hover.slab : null);
    if (hover && hover.sticky) showInfo(hover);
  }

  function slabSummary(slab) {
    const group = display && display.slabs.find(function (s) { return s.slab === slab; });
    const data = dataFor(stage.at);
    let generated = 0;
    let total = 0;
    (group ? group.blocks : []).forEach(function (block) {
      generated += generatedOf(block, data);
      total += block.total || 0;
    });
    return { generated: generated, total: total };
  }

  function slabText(at, slab) {
    const layers = S.slabLayers(at, slab);
    const z0 = (layers.first - 0.5) * edgePc;
    const z1 = (layers.last + 0.5) * edgePc;
    const where = z1 <= 0 ? Math.round(-z1) + "–" + Math.round(-z0) + " pc below the plane"
      : z0 >= 0 ? Math.round(z0) + "–" + Math.round(z1) + " pc above the plane"
        : "across the plane";
    const span = layers.first === layers.last ? "layer " + layers.first : "layers " + layers.first + " to " + layers.last;
    return S.slabNoun(at) + " " + slab + ": " + span + ", " + where;
  }

  function blockText(block, data) {
    const b = block.bounds;
    const bearing = S.bearingRange(b);
    const generated = generatedOf(block, data);
    const parts = [(block.m === 1 ? "Sector " : "Block ") + S.blockLabel(block),
      "bearing " + bearing[0].toFixed(1) + "°–" + bearing[1].toFixed(1) + "°",
      Math.round(b.r0) + "–" + Math.round(b.r1) + " pc from the core"];
    if (block.m === 1) parts.push(generated > 0 ? "generated" : "not generated");
    else parts.push(S.formatInt(generated) + (block.total != null ? " of " + S.formatInt(block.total) : "") + " sectors generated");
    return parts.join(", ");
  }

  function showTooltip(text, clientX, clientY) {
    const tip = els.tooltip;
    if (!tip) return;
    if (!text) {
      tip.hidden = true;
      return;
    }
    const rect = tip.parentElement.getBoundingClientRect();
    tip.textContent = text;
    tip.hidden = false;
    const x = Math.min(clientX - rect.left + 12, rect.width - tip.offsetWidth - 4);
    const y = Math.min(clientY - rect.top + 12, rect.height - tip.offsetHeight - 4);
    tip.style.left = Math.max(4, x) + "px";
    tip.style.top = Math.max(4, y) + "px";
  }

  function hoverText(found) {
    if (stage.slab == null) {
      const summary = slabSummary(found.slab);
      return slabText(stage.at, found.slab) + ", " + S.formatInt(summary.generated)
        + (getOutline().shapeless ? "" : " of " + S.formatInt(summary.total)) + " sectors generated";
    }
    return blockText(found.block, dataFor(stage.at));
  }

  // The info panel for the selected block (top-down stages; a 3D stage's
  // slabs have the strip).
  function showInfo(h) {
    const data = dataFor(stage.at);
    if (stage.slab == null || !h.block) return;
    const block = h.block;
    if (block.m === 1) {
      const sector = data && data.sectors.get(block.ring + "/" + block.wedge + "/" + block.slab);
      if (sector) {
        host.showPlacedInfo(sectorEntry(block, sector));
      } else {
        host.showCellInfo({ m: 1, bounds: block.bounds, address: { ring: block.ring, layer: block.slab, slot: block.wedge }, filled: 0 });
      }
      return;
    }
    host.showBlockInfo({
      block: block, generated: generatedOf(block, data), total: block.total,
      enter: function () { go({ at: S.stripBlock(block), slab: null }); },
    });
  }

  // The panel on arrival: the stage's own block (its hint under it), or
  // just the hint at the galaxy.
  function showStageInfo() {
    const hint = hintFor(stage);
    if (!stage.at) {
      host.showHint(hint);
      return;
    }
    // Its generated sectors are its children's.
    const data = dataFor(stage.at);
    let generated = 0;
    if (data) data.generated.forEach(function (n) { generated += n; });
    host.showBlockInfo({
      block: stage.at, total: getOutline().shapeless ? null : S.drillBlockTotal(stage.at, getOutline()),
      generated: data ? generated : null, hint: hint, enter: null,
    });
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

  // A stage's own hint in the info panel.
  function hintFor(s) {
    if (s.slab == null) {
      return "Hover a " + S.slabNoun(s.at).toLowerCase() + " to see it, click it (or pick it in the list) to pull it out"
        + " and look at it from above. Drag to turn the view.";
    }
    if (s.at && s.at.m === 3) {
      return "Click a sector to open it; a sector that isn't generated yet shows where it is"
        + (host.canGenerate ? " and how to generate it." : ".");
    }
    return "Click a block to fly into it.";
  }

  // --- Acting on a block or slab -------------------------------------------

  function act(found) {
    if (stage.slab == null) {
      go({ at: stage.at, slab: found.slab });
      return;
    }
    const block = found.block;
    if (block.m === 1) {
      const data = dataFor(stage.at);
      const sector = data && data.sectors.get(block.ring + "/" + block.wedge + "/" + block.slab);
      if (sector && host.sectorUrl(sector.id)) {
        window.location.assign(host.sectorUrl(sector.id));
        return;
      }
      selectedSector = { ring: block.ring, layer: block.slab, slot: block.wedge };
      setHover({ block: block, slab: block.slab, sticky: true });
      renderCrumbs();
      return;
    }
    go({ at: S.stripBlock(block), slab: null });
  }

  function sectorBlock(sector) {
    if (!display) return null;
    for (const slab of display.slabs) {
      const found = slab.blocks.find(function (b) {
        return b.ring === sector.ring && b.wedge === sector.slot && b.slab === sector.layer;
      });
      if (found) return found;
    }
    return null;
  }

  // --- Input -------------------------------------------------------------------

  let pointer = null;
  let touchPending = null;

  function onPointerDown(event) {
    if (event.button !== 0) return;
    pointer = { x: event.clientX, y: event.clientY, moved: 0, id: event.pointerId, type: event.pointerType };
    try { canvasEl.setPointerCapture(event.pointerId); } catch (err) { /* not essential */ }
  }

  function onPointerMove(event) {
    if (pointer && event.pointerId === pointer.id) {
      const dx = event.clientX - pointer.x;
      const dy = event.clientY - pointer.y;
      pointer.moved += Math.abs(dx) + Math.abs(dy);
      pointer.x = event.clientX;
      pointer.y = event.clientY;
      if (pointer.moved > DRAG_CLICK_PX && !animation && view) {
        drag(dx, dy);
        showTooltip("", 0, 0);
      }
      return;
    }
    if (animation || event.pointerType === "touch") return;
    const found = blockAt(event.clientX, event.clientY);
    if (!found) {
      if (hover && !hover.sticky) setHover(null);
      showTooltip("", 0, 0);
      return;
    }
    if (!hover || hover.sticky || hover.slab !== found.slab || hover.block !== found.block) {
      setHover({ block: found.block, slab: found.slab });
    }
    showTooltip(hoverText(found), event.clientX, event.clientY);
  }

  function onPointerUp(event) {
    if (!pointer || event.pointerId !== pointer.id) return;
    const wasClick = pointer.moved <= DRAG_CLICK_PX;
    const type = pointer.type;
    pointer = null;
    try { canvasEl.releasePointerCapture(event.pointerId); } catch (err) { /* already released */ }
    if (!wasClick || animation) return;
    const found = blockAt(event.clientX, event.clientY);
    if (!found) return;
    if (type === "touch") {
      // First tap highlights, a second on the same thing acts.
      const key = stage.slab == null ? "slab:" + found.slab : "block:" + S.blockLabel(found.block);
      if (touchPending !== key) {
        touchPending = key;
        setHover({ block: found.block, slab: found.slab, sticky: true });
        showTooltip(hoverText(found), event.clientX, event.clientY);
        return;
      }
      touchPending = null;
    }
    showTooltip("", 0, 0);
    act(found);
  }

  function onPointerLeave() {
    showTooltip("", 0, 0);
    if (hover && !hover.sticky) setHover(null);
  }

  // Drag: turns a 3D stage; pans the galaxy's slab (stage 2) when zoomed
  // in; does nothing on the other top-down stages, which fit the canvas.
  function drag(dx, dy) {
    if (stage.slab == null) {
      view.theta -= dx * ROTATE_PER_PX;
      view.phi = Math.max(MIN_TILT, Math.min(MAX_TILT, view.phi - dy * ROTATE_PER_PX));
      applyView();
      return;
    }
    if (stage.at) return;
    const heightPx = canvasEl.clientHeight || 1;
    const perPx = (2 * view.dist * Math.tan(THREE.MathUtils.degToRad(camera.fov) / 2)) / heightPx;
    const right = new THREE.Vector3().setFromMatrixColumn(camera.matrixWorld, 0);
    const upward = new THREE.Vector3().setFromMatrixColumn(camera.matrixWorld, 1);
    view.target[0] -= (right.x * dx - upward.x * dy) * perPx;
    view.target[1] -= (right.y * dx - upward.y * dy) * perPx;
    const reach = host.galaxyRadius;
    const r = Math.hypot(view.target[0], view.target[1]);
    if (r > reach) {
      view.target[0] *= reach / r;
      view.target[1] *= reach / r;
    }
    applyView();
  }

  function zoomBy(factor) {
    if (!view || animation) return;
    const closest = view.fit / (stage.slab == null ? ZOOM_3D : stage.at ? 1 : ZOOM_GALAXY_SLAB);
    view.dist = Math.max(closest, Math.min(view.fit, view.dist * factor));
    applyView();
  }

  function onWheel(event) {
    let deltaPx = event.deltaY;
    if (event.deltaMode === 1) deltaPx *= 33;
    else if (event.deltaMode === 2) deltaPx *= canvasEl.clientHeight || 400;
    deltaPx = Math.max(-200, Math.min(200, deltaPx));
    if (deltaPx) zoomBy(Math.exp(deltaPx * WHEEL_ZOOM_PER_PX));
  }

  // The blocks of the top-down slab in ring order, for the arrow keys.
  function slabBlocks() {
    if (!display || !display.slabs.length) return [];
    return display.slabs[0].blocks.filter(function (b) {
      return !generatedOnly || generatedOf(b, dataFor(stage.at)) > 0;
    });
  }

  function onKey(event) {
    const key = event.key;
    if (key === "Escape" || key === "Backspace") {
      event.preventDefault();
      up();
      return;
    }
    if (key === "Home") {
      event.preventDefault();
      go({ at: null, slab: null });
      return;
    }
    if (["ArrowLeft", "ArrowRight", "ArrowUp", "ArrowDown", "Enter"].indexOf(key) < 0) return;
    event.preventDefault();
    if (animation || !display) return;
    if (stage.slab == null) {
      const slabs = display.slabs.map(function (s) { return s.slab; });
      if (key === "Enter") {
        if (hover && hover.slab != null) go({ at: stage.at, slab: hover.slab });
        return;
      }
      if (key === "ArrowUp" || key === "ArrowDown") {
        const n = hover && hover.slab != null ? slabs.indexOf(hover.slab) : -1;
        const next = n < 0 ? (key === "ArrowUp" ? slabs.length - 1 : 0) : Math.max(0, Math.min(slabs.length - 1, n + (key === "ArrowUp" ? 1 : -1)));
        setHover({ slab: slabs[next], sticky: true });
        return;
      }
      // Left/right: the neighbouring block, or turn the galaxy.
      if (stage.at) {
        const siblings = S.crumbSiblings(stage, getOutline(), edgePc);
        const n = siblings.findIndex(function (s) { return S.sameStage(s.stage, stage); });
        const next = siblings[(n + (key === "ArrowRight" ? 1 : -1) + siblings.length) % siblings.length];
        if (next) go(next.stage);
      } else {
        view.theta += key === "ArrowLeft" ? 0.1 : -0.1;
        applyView();
      }
      return;
    }
    const blocks = slabBlocks();
    if (!blocks.length) return;
    let current = hover && hover.block ? hover.block : null;
    if (key === "Enter") {
      if (current) act({ block: current, slab: current.slab });
      return;
    }
    if (!current) {
      setHover({ block: blocks[0], slab: blocks[0].slab, sticky: true });
      return;
    }
    let next = null;
    if (key === "ArrowLeft" || key === "ArrowRight") {
      const same = blocks.filter(function (b) { return b.ring === current.ring; });
      const n = same.indexOf(current);
      next = same[(n + (key === "ArrowLeft" ? 1 : -1) + same.length) % same.length];
    } else {
      const ring = current.ring + (key === "ArrowUp" ? 1 : -1);
      const mid = (current.bounds.t0 + current.bounds.t1) / 2;
      blocks.forEach(function (b) {
        if (b.ring !== ring) return;
        const d = Math.abs(wrapAngle((b.bounds.t0 + b.bounds.t1) / 2 - mid));
        if (!next || d < next.d) next = { block: b, d: d };
      });
      next = next ? next.block : null;
    }
    if (next) setHover({ block: next, slab: next.slab, sticky: true });
  }

  // --- Breadcrumb, slab strip, notice ------------------------------------

  function notice(text) {
    if (!els.notice) return;
    els.notice.textContent = text || "";
    els.notice.hidden = !text;
  }

  function renderCrumbs() {
    const nav = els.crumbs;
    if (!nav) return;
    nav.textContent = "";
    const list = document.createElement("ol");
    const items = S.crumbs(stage);
    if (selectedSector) {
      items.push({ label: "Sector " + S.blockLabel({ m: 1, ring: selectedSector.ring, wedge: selectedSector.slot, slab: selectedSector.layer }), last: true, sector: true });
      items[items.length - 2].last = false;
    }
    items.forEach(function (crumb) {
      const item = document.createElement("li");
      if (crumb.last) {
        const here = document.createElement("span");
        here.textContent = crumb.label;
        here.setAttribute("aria-current", "location");
        item.appendChild(here);
      } else {
        const button = document.createElement("button");
        button.type = "button";
        button.className = "galaxy-crumb";
        button.textContent = crumb.label;
        button.addEventListener("click", function () { go(crumb.stage); });
        item.appendChild(button);
      }
      if (crumb.stage && (crumb.stage.at || crumb.stage.slab != null)) {
        const siblings = S.crumbSiblings(crumb.stage, getOutline(), edgePc);
        if (siblings.length > 1) item.appendChild(siblingMenu(crumb, siblings));
      }
      list.appendChild(item);
    });
    nav.appendChild(list);
  }

  function siblingMenu(crumb, siblings) {
    const menu = document.createElement("details");
    menu.className = "galaxy-crumb-menu";
    const summary = document.createElement("summary");
    summary.textContent = "▾";
    summary.setAttribute("aria-label", "Others beside " + crumb.label);
    menu.appendChild(summary);
    const list = document.createElement("ul");
    siblings.forEach(function (sibling) {
      const item = document.createElement("li");
      const button = document.createElement("button");
      button.type = "button";
      button.textContent = sibling.label;
      if (S.sameStage(sibling.stage, crumb.stage)) button.setAttribute("aria-current", "true");
      button.addEventListener("click", function () {
        menu.open = false;
        go(sibling.stage);
      });
      item.appendChild(button);
      list.appendChild(item);
    });
    menu.appendChild(list);
    return menu;
  }

  // The slab strip: one row per child slab of the container, top first,
  // with its generated share. In 3D, hovering a row highlights the slab
  // and clicking pulls it out; top-down, clicking moves to that slab.
  function renderStrip() {
    const box = els.slabs;
    if (!box) return;
    box.textContent = "";
    const heading = document.createElement("h3");
    heading.textContent = S.slabNoun(stage.at) + "s";
    box.appendChild(heading);
    const data = dataFor(stage.at);
    const groups = childrenOf(stage.at).slice().reverse();
    const list = document.createElement("ul");
    groups.forEach(function (g) {
      if (getOutline().shapeless && !g.blocks.some(function (b) { return generatedOf(b, data) > 0; })) return;
      let generated = 0;
      let total = 0;
      g.blocks.forEach(function (block) {
        generated += generatedOf(block, data);
        total += block.total || 0;
      });
      const item = document.createElement("li");
      const button = document.createElement("button");
      button.type = "button";
      button.className = "galaxy-slab-row";
      button.dataset.slab = String(g.slab);
      const layers = S.slabLayers(stage.at, g.slab);
      const name = document.createElement("span");
      name.textContent = S.slabNoun(stage.at) + " " + g.slab
        + (layers.first === layers.last ? "" : " · layers " + layers.first + " to " + layers.last);
      const bar = document.createElement("span");
      bar.className = "galaxy-slab-bar";
      const fill = document.createElement("span");
      fill.style.width = total > 0 ? Math.max(generated > 0 ? 2 : 0, (100 * generated) / total).toFixed(1) + "%" : "0%";
      bar.appendChild(fill);
      const count = document.createElement("span");
      count.className = "galaxy-slab-count";
      count.textContent = S.formatInt(generated) + (total > 0 ? " / " + S.formatInt(total) : "");
      button.appendChild(name);
      button.appendChild(bar);
      button.appendChild(count);
      button.setAttribute("aria-label", slabText(stage.at, g.slab) + ", " + S.formatInt(generated)
        + (total > 0 ? " of " + S.formatInt(total) : "") + " sectors generated");
      if (stage.slab === g.slab) button.setAttribute("aria-current", "true");
      button.addEventListener("mouseenter", function () {
        if (stage.slab == null && !animation) setHover({ slab: g.slab });
      });
      button.addEventListener("focus", function () {
        if (stage.slab == null && !animation) setHover({ slab: g.slab, sticky: true });
      });
      button.addEventListener("click", function () {
        go({ at: stage.at, slab: g.slab });
      });
      item.appendChild(button);
      list.appendChild(item);
    });
    box.appendChild(list);
  }

  function markStripRow(slab) {
    if (!els.slabs) return;
    els.slabs.querySelectorAll(".galaxy-slab-row").forEach(function (row) {
      row.classList.toggle("is-hovered", slab != null && row.dataset.slab === String(slab));
    });
  }

  // --- The address bar (section 9.3) -----------------------------------------

  // Flies to sector {ring, layer, slot}'s stage 8 and selects it.
  function locate(sector) {
    const next = S.sectorStage(sector.ring, sector.layer, sector.slot);
    const problem = S.validStage(next, getOutline(), edgePc);
    if (problem) {
      notice("Sector " + S.blockLabel({ m: 1, ring: sector.ring, wedge: sector.slot, slab: sector.layer })
        + " is outside the galaxy.");
      return false;
    }
    selectedSector = sector;
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

  // --- Turning on and off ------------------------------------------------------

  function setActive(on) {
    if (on === active) return;
    active = on;
    [els.crumbs, els.slabs, els.address].forEach(function (el) { if (el) el.hidden = !on; });
    if (!on) clearMatches();
    if (!on) {
      if (animation) finishAnimation();
      disposeDisplay(display);
      display = null;
      view = null;
      clearOutline();
      showTooltip("", 0, 0);
      notice("");
      return;
    }
    renderCrumbs();
    loadData(stage.at).then(function () {
      if (!active) return;
      display = null;
      flyTo(stage);
      renderStrip();
    });
  }

  // Opens the stage the URL asks for (or the galaxy, with a notice).
  function openFromLocation(push) {
    const parsed = S.parseStageQuery(location.search);
    let next = parsed.stage;
    let problem = parsed.problem || S.validStage(next, getOutline(), edgePc);
    if (problem) next = { at: null, slab: null };
    selectedSector = parsed.sector && !problem ? parsed.sector : null;
    stage = next;
    renderCrumbs();
    renderStrip();
    host.setBlockSize(childM(next.at));
    if (push === false && display) {
      go(next, { push: false, keepSector: true });
    }
    notice(problem || "");
  }

  function onPopState() {
    if (!active) return;
    const parsed = S.parseStageQuery(location.search);
    const problem = parsed.problem || S.validStage(parsed.stage, getOutline(), edgePc);
    selectedSector = parsed.sector && !problem ? parsed.sector : null;
    go(problem ? { at: null, slab: null } : parsed.stage, { push: false, keepSector: true });
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
    onPointerDown: onPointerDown,
    onPointerMove: onPointerMove,
    onPointerUp: onPointerUp,
    onPointerLeave: onPointerLeave,
    onWheel: onWheel,
    onKey: onKey,
    zoomBy: zoomBy,
    home: function () { go({ at: null, slab: null }); },
    up: up,
    setGeneratedOnly: setGeneratedOnly,
    stage: function () { return stage; },
    go: go,
    locate: locate,
  };
}
