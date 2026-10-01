// html/static/galaxystageview.js
//
// The Galaxy Map's drill-down (docs/design/galaxy-drilldown-navigation.md,
// sections 4, 5, 8.1 and 10-11), drawn and driven. Always from straight
// above, with no free camera: the galaxy's quarters, then a layer of the
// disk (picked from the strip beside the map), then an arc of the ring
// band in view (clicked on the map), then a layer, an arc, ... down to a
// sector. galaxystages.js has the rules; this file has the scene, the
// camera moves, the breadcrumb, the layer strip, the tooltip, keys,
// touch, the stage URLs and the map's own Back and Forward.
// galaxymap3d.js creates it (createStageView), hands it the pointer and
// key events and calls step() every frame.
//
// host (from galaxymap3d.js):
// - THREE, scene, camera, canvasEl;
// - edgePc, shape (or null), galaxyRadius, reducedMotion, accentColor;
// - blockScene: galaxyblocks' scene (buildCells);
// - makeBlockMesh(part, translucent): a mesh with the map's block shader;
// - setCamera({target: [x, y, z], dist, theta, phi}): moves the camera
//   (and the tiles, scale bar and wedge lines with it);
// - setWedgeClip(clip): the wedge lines kept to {r0, r1, a0, a1} (null:
//   the whole galaxy, with its bearing labels);
// - fetchStage(query): a Promise of GET /galaxy/stage's JSON for a
//   container ("?at=m.ring.wedge.slab", or "" for the galaxy);
// - setBlockSize(m): the scale readout's "1 block = m sectors";
// - showBlockInfo(info), showPlacedInfo(entry), showCellInfo(cell),
//   showHint(text): the info panel (showBlockInfo's info: {block, total,
//   generated, hint, enter, generate});
// - showPointAt(clientX, clientY): shows a bright star or cloud under a
//   click, true when there was one;
// - canGenerate, courseSectors, sectorUrl(id), locate(name);
// - els: {crumbs, slabs, tooltip, notice, address, matches, controls}
//   (any may be missing).

const VERSION_QUERY = new URL(import.meta.url).search;

const S = await import(`./galaxystages.js${VERSION_QUERY}`);

// Not quite 0: the camera keeps galactic north as its up vector, which
// needs the view direction off vertical by a hair.
const TOP_DOWN_PHI = 1e-3;
// The other choices while one is hovered (MAP.18).
const OTHER_FADE = 0.25;
// A new stage's blocks fade in over the last part of a flight.
const FADE_IN_MS = 200;
const DRAG_CLICK_PX = 6;
// A small cube of sectors (a level-3 block, 27 at most) is shown from a
// fixed slant with its layers pulled apart, so each sector can be picked
// on the map (Boss, 2026-10-01). Nothing rotates. The tilt from straight
// down, and the gap between layers in sector heights.
const CUBE_MAX_SECTORS = 27;
const CUBE_TILT = (55 * Math.PI) / 180;
const CUBE_LAYER_GAP = 2;
const TWO_PI = 2 * Math.PI;

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
        (payload.children || []).forEach(function (child) {
          generated.set(child.ring + "/" + child.wedge + "/" + child.slab, child.generated);
        });
        const sectors = new Map();
        (payload.sectors || []).forEach(function (sector) {
          sectors.set(sector.ring + "/" + sector.slot + "/" + sector.layer, sector);
        });
        entry.value = { generated: generated, sectors: sectors, sectorGenerated: new Map() };
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
        entry.value = { generated: new Map(), sectors: new Map(), sectorGenerated: new Map(), failed: true };
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
    let options = r.kind ? r.options : [{ pick: null, blocks: r.view.blocks, a0: r.view.a0, a1: r.view.a1 }];
    if (isCube(r)) {
      // Every sector is a choice of its own; the strip still offers the
      // layers (r.options).
      options = r.view.blocks.map(function (block) {
        return { pick: null, blocks: [block], a0: r.view.a0, a1: r.view.a1 };
      });
    }
    if (getOutline().shapeless) {
      options = options.map(function (o) {
        return Object.assign({}, o, { blocks: o.blocks.filter(function (b) { return generatedOf(b, data) > 0; }) });
      }).filter(function (o) { return o.blocks.length; });
    }
    return options;
  }

  function isSectorView(r) {
    return !!(r && r.view && r.view.blocks.length && r.view.blocks[0].m === 1);
  }

  // A view of sectors across several layers, few enough to show as a
  // cube and pick one by one.
  function isCube(r) {
    return !!(r && r.kind === "layer" && isSectorView(r) && r.view.blocks.length <= CUBE_MAX_SECTORS);
  }

  // How far a block is lifted in the cube (0 elsewhere): its layer's
  // distance from the cube's middle layer, spread by CUBE_LAYER_GAP.
  function liftOf(r, block) {
    if (!isCube(r)) return 0;
    let lo = Infinity;
    let hi = -Infinity;
    r.view.blocks.forEach(function (b) {
      lo = Math.min(lo, b.slab);
      hi = Math.max(hi, b.slab);
    });
    return (block.slab - (lo + hi) / 2) * CUBE_LAYER_GAP * edgePc;
  }

  function isWholeGalaxy(r) {
    return !r.stage.at && r.view.a1 - r.view.a0 >= TWO_PI - 1e-9;
  }

  // The bounds a set of blocks covers in the plane, bearings counted
  // either way from `mid` (the middle of the view's bearings, so a block
  // reaching a little before the view's first bearing doesn't wrap all the
  // way round): {r0, r1, t0, t1, z1}.
  function spanOf(blocks, mid) {
    const span = { r0: Infinity, r1: 0, t0: Infinity, t1: -Infinity, z1: -Infinity };
    blocks.forEach(function (block) {
      const b = block.bounds;
      const t0 = mid + wrapAngle(b.t0 - mid);
      span.r0 = Math.min(span.r0, b.r0);
      span.r1 = Math.max(span.r1, b.r1);
      span.t0 = Math.min(span.t0, t0);
      span.t1 = Math.max(span.t1, t0 + (b.t1 - b.t0));
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
          filled: generatedOf(block, data), total: block.total, block: block, option: index,
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
      group.position.z = liftOf(r, option.blocks[0]);
      root.add(group);
      groups.push({ option: option, meshes: meshes, fade: 1 });
    });
    host.scene.add(root);
    root.updateMatrixWorld(true);
    return { resolved: r, root: root, options: options, groups: groups, fade: 1, data: data };
  }

  function disposeDisplay(d) {
    if (!d) return;
    host.scene.remove(d.root);
    d.groups.forEach(function (group) {
      group.meshes.forEach(function (mesh) {
        mesh.geometry.dispose();
        mesh.material.dispose();
      });
    });
  }

  function setDisplayFade(d, fade) {
    d.fade = fade;
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

  // The camera for a stage, from straight above with the view's middle
  // bearing pointing up the screen (galactic north up for the whole
  // galaxy): {target, dist, theta, phi}.
  function cameraFor(r) {
    const blocks = r.view.blocks;
    if (!blocks.length) {
      return { target: [0, 0, 0], dist: host.galaxyRadius * 2.4, theta: -Math.PI / 2, phi: TOP_DOWN_PHI };
    }
    const fp = S.footprint(blocks);
    let z0 = Infinity;
    let z1 = -Infinity;
    blocks.forEach(function (block) {
      z0 = Math.min(z0, block.bounds.z0);
      z1 = Math.max(z1, block.bounds.z1);
    });
    const theta = isWholeGalaxy(r) ? -Math.PI / 2 : (r.view.a0 + r.view.a1) / 2 + Math.PI;
    if (isCube(r)) {
      blocks.forEach(function (block) {
        z0 = Math.min(z0, block.bounds.z0 + liftOf(r, block));
        z1 = Math.max(z1, block.bounds.z1 + liftOf(r, block));
      });
      return {
        target: [fp.center[0], fp.center[1], (z0 + z1) / 2],
        dist: (S.FIT_MARGIN * Math.hypot(fp.radius, (z1 - z0) / 2)) / Math.sin(fovHalf()),
        theta: theta, phi: CUBE_TILT,
      };
    }
    return {
      target: [fp.center[0], fp.center[1], (z0 + z1) / 2],
      dist: (S.FIT_MARGIN * fp.radius) / Math.tan(fovHalf()) + (z1 - z0) / 2,
      theta: theta, phi: TOP_DOWN_PHI,
    };
  }

  function applyView() {
    host.setCamera({ target: view.target.slice(), dist: view.dist, theta: view.theta, phi: view.phi });
  }

  function wrapAngle(a) {
    return Math.atan2(Math.sin(a), Math.cos(a));
  }

  // Animates the camera from `from` to `to` along van Wijk and Nuij's
  // path in the plane, easing the height and turn alongside; calls
  // progress(t) each frame (t 0..1) and done() at the end.
  function flyCamera(from, to, progress, done) {
    const tanHalf = Math.tan(THREE.MathUtils.degToRad(camera.fov) / 2);
    const path = S.flightPath(from.target, 2 * from.dist * tanHalf, to.target, 2 * to.dist * tanHalf);
    const turn = wrapAngle(to.theta - from.theta);
    animation = {
      startedAt: performance.now(),
      duration: host.reducedMotion ? 0 : S.flightMs(path.S),
      frame: function (t) {
        const e = S.easeInOut(t);
        const p = path.at(path.S * t);
        view = {
          target: [p.center[0], p.center[1], from.target[2] + (to.target[2] - from.target[2]) * e],
          dist: p.w / (2 * tanHalf),
          theta: from.theta + turn * e,
          phi: from.phi + (to.phi - from.phi) * e,
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
    host.setWedgeClip({ r0: span.r0, r1: span.r1, a0: span.t0, a1: span.t1 });
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
    const query = options.query || S.stageQuery(r.stage);
    if (options.push !== false && (!S.sameStage(stage, r.stage) || query !== location.search)) {
      mapIndex += 1;
      maxIndex = mapIndex;
      history.pushState({ galaxyStage: true, mapIndex: mapIndex }, "", location.pathname + query + location.hash);
    }
    stage = r.stage;
    resolved = r;
    renderCrumbs();
    renderStrip();
    renderTravel();
    hover = null;
    showTooltip("", 0, 0);
    host.setBlockSize(r.view.blocks.length ? r.view.blocks[0].m : 243);
    loadData(r.stage.at).then(function () {
      if (!active || token !== goToken) return;
      if (animation) finishAnimation();
      flyTo(r);
    });
  }

  function flyTo(r) {
    const from = display && view ? { target: view.target.slice(), dist: view.dist, theta: view.theta, phi: view.phi } : null;
    const to = cameraFor(r);
    const old = display;
    const incoming = buildDisplay(r);
    setWedgeClip(r);
    setDisplayFade(incoming, from ? 0 : 1);
    display = incoming;
    if (!from) {
      view = to;
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
      view = to;
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
  }

  // --- Picking and hover -----------------------------------------------------

  const raycaster = new THREE.Raycaster();

  // Whether choice `index` can be picked: not while it holds nothing
  // generated with "Generated only" on.
  function pickable(index) {
    if (!display || !display.options[index]) return false;
    if (!generatedOnly) return true;
    const data = display.data;
    return display.options[index].blocks.some(function (b) { return generatedOf(b, data) > 0; });
  }

  // The choice under a screen point (an index into display.options), or
  // -1. Layers are picked from the strip, not the map: from above, one
  // covers the others.
  function optionAt(clientX, clientY) {
    if (!display || !resolved || (resolved.kind === "layer" && !isCube(display.resolved))) return -1;
    const rect = canvasEl.getBoundingClientRect();
    if (!rect.width || !rect.height) return -1;
    const ndc = new THREE.Vector2(((clientX - rect.left) / rect.width) * 2 - 1, -((clientY - rect.top) / rect.height) * 2 + 1);
    raycaster.setFromCamera(ndc, camera);
    const meshes = [];
    display.groups.forEach(function (group) { Array.prototype.push.apply(meshes, group.meshes); });
    const hits = raycaster.intersectObjects(meshes, false);
    for (let n = 0; n < hits.length; n++) {
      const hit = hits[n];
      const cell = hit.object.userData.cells[hit.object.userData.part.owners[hit.face.a]];
      if (!cell) continue;
      if (generatedOnly && !(cell.filled > 0)) continue;
      return cell.option;
    }
    return -1;
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

  // An accent outline around a choice's area, seen from above (MAP.18).
  function outlineSpan(span) {
    clearOutline();
    const steps = Math.max(2, Math.ceil((span.t1 - span.t0) / (Math.PI / 90)));
    const points = [];
    for (let k = 0; k <= steps; k++) {
      const t = span.t0 + ((span.t1 - span.t0) * k) / steps;
      points.push(new THREE.Vector3(span.r1 * Math.cos(t), span.r1 * Math.sin(t), span.z1));
    }
    for (let k = steps; k >= 0; k--) {
      const t = span.t0 + ((span.t1 - span.t0) * k) / steps;
      points.push(new THREE.Vector3(span.r0 * Math.cos(t), span.r0 * Math.sin(t), span.z1));
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

  // The hovered choice stays as it is and the others fade (a layer
  // hovered in the strip too: in the cube, its sectors stay); a choice on
  // the map gets an outline.
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
    if (option && (resolved.kind !== "layer" || isCube(resolved))) {
      const span = spanOf(option.blocks, (option.a0 + option.a1) / 2);
      if (option.blocks.length === 1) {
        const b = option.blocks[0].bounds;
        span.t0 = b.t0;
        span.t1 = b.t1;
      }
      span.z1 += liftOf(resolved, option.blocks[0]);
      outlineSpan(span);
    } else {
      clearOutline();
    }
    markStripRow(layer ? layer : option && isCube(resolved) ? { lo: option.blocks[0].slab, hi: option.blocks[0].slab } : null);
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

  // The tooltip beside a choice picked by keyboard: at its middle on
  // screen.
  function tooltipAtOption(index) {
    const option = display.options[index];
    const fp = S.footprint(option.blocks);
    const point = new THREE.Vector3(fp.center[0], fp.center[1], option.blocks[0].bounds.z1 + liftOf(resolved, option.blocks[0])).project(camera);
    const rect = canvasEl.getBoundingClientRect();
    showTooltip(optionText(index), rect.left + ((point.x + 1) / 2) * rect.width, rect.top + ((1 - point.y) / 2) * rect.height);
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
    if (!r || !r.kind) return "";
    if (isCube(r)) {
      return "Click a sector to open it (a sector that isn't generated yet shows where it is"
        + (host.canGenerate ? " and how to generate it" : "") + "), or pick a layer from the list to see just that one.";
    }
    if (r.kind === "layer") {
      return "Pick a " + S.slabNoun(r.stage.at).toLowerCase() + " (a layer of the disk) from the list beside the map.";
    }
    if (r.kind === "quadrant") return "Click a quarter of the galaxy to look at it more closely.";
    if (isSectorView(r)) {
      return "Click a sector to open it; a sector that isn't generated yet shows where it is"
        + (host.canGenerate ? " and how to generate it." : ".");
    }
    return "Click an arc of the map to zoom into it.";
  }

  // --- Acting on a choice ----------------------------------------------------

  function act(index) {
    if (!pickable(index)) return;
    const option = display.options[index];
    if (option.blocks.length === 1 && option.blocks[0].m === 1) {
      const block = option.blocks[0];
      const sector = sectorRecord(block);
      if (sector && host.sectorUrl(sector.id)) {
        window.location.assign(host.sectorUrl(sector.id));
        return;
      }
      selectedSector = { ring: block.ring, layer: block.slab, slot: block.wedge };
      if (isCube(resolved)) {
        // Picked in the cube: on to that sector's own layer, selected.
        const layerPick = resolved.options.find(function (o) {
          return block.slab >= o.pick.lo && block.slab <= o.pick.hi;
        });
        if (layerPick) {
          go({ at: stage.at, picks: stage.picks.concat([layerPick.pick]) }, { keepSector: true });
          return;
        }
      }
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

  let pointer = null;
  let touchPending = -1;

  function onPointerDown(event) {
    if (event.button !== 0) return;
    pointer = { x0: event.clientX, y0: event.clientY, id: event.pointerId, type: event.pointerType };
  }

  function onPointerMove(event) {
    if (animation || event.pointerType === "touch" || (pointer && event.pointerId === pointer.id)) return;
    const index = optionAt(event.clientX, event.clientY);
    if (index < 0) {
      if (hover && !hover.sticky) setHover(null);
      showTooltip("", 0, 0);
      return;
    }
    if (!hover || hover.option !== index) setHover({ option: index });
    showTooltip(optionText(index), event.clientX, event.clientY);
  }

  function onPointerUp(event) {
    if (!pointer || event.pointerId !== pointer.id) return;
    const moved = Math.abs(event.clientX - pointer.x0) + Math.abs(event.clientY - pointer.y0);
    const type = pointer.type;
    pointer = null;
    if (moved > DRAG_CLICK_PX || animation) return;
    // Inside a container, a bright star or cloud under the click is shown
    // rather than the block picked (over the whole galaxy and its
    // quarters the stars are too thick for that).
    if (stage.at && host.showPointAt && host.showPointAt(event.clientX, event.clientY)) {
      showTooltip("", 0, 0);
      return;
    }
    const index = optionAt(event.clientX, event.clientY);
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

  function onPointerLeave() {
    showTooltip("", 0, 0);
    if (hover && !hover.sticky) setHover(null);
  }

  // Arrow keys move among the choices: Left and Right through them in
  // order, Up and Down to the next one along the same arc (or the next
  // layer up or down); Enter takes it.
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
    } else if (isCube(resolved) && (key === "ArrowUp" || key === "ArrowDown")) {
      // The same column, one layer up or down.
      const here = options[current].blocks[0];
      const slab = here.slab + (key === "ArrowUp" ? 1 : -1);
      next = options.findIndex(function (o) {
        const b = o.blocks[0];
        return b.ring === here.ring && b.wedge === here.wedge && b.slab === slab;
      });
    } else if (resolved.kind === "layer" || key === "ArrowLeft" || key === "ArrowRight") {
      const forward = key === "ArrowRight" || key === "ArrowUp";
      next = (current + (forward ? 1 : -1) + options.length) % options.length;
    } else {
      const a0 = options[current].a0;
      const order = key === "ArrowUp" ? 1 : -1;
      for (let n = current + order; n >= 0 && n < options.length; n += order) {
        if (Math.abs(options[n].a0 - a0) < 1e-9) {
          next = n;
          break;
        }
      }
    }
    if (next < 0) return;
    setHover({ option: next, sticky: true });
    tooltipAtOption(next);
  }

  // --- Breadcrumb, layer strip, notice ---------------------------------------

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
    const items = S.crumbs(stage, getOutline(), edgePc);
    if (selectedSector && !(resolved && resolved.sector)) {
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
      list.appendChild(item);
    });
    nav.appendChild(list);
  }

  // The layer strip (MAP.17): while the next pick is a layer, one row per
  // choice, top first, with its generated share; hovering a row fades the
  // others on the map, clicking takes it. Otherwise it says which layers
  // the view holds.
  function renderStrip() {
    const box = els.slabs;
    if (!box || !resolved || !resolved.view) return;
    box.textContent = "";
    const heading = document.createElement("h3");
    // A thin block's view is sectors, so its "slabs" are sector layers.
    const noun = isSectorView(resolved) ? "Layer" : S.slabNoun(stage.at);
    heading.textContent = noun + "s";
    box.appendChild(heading);
    const slabs = Array.from(new Set(resolved.view.blocks.map(function (b) { return b.slab; }))).sort(function (p, q) { return p - q; });
    if (resolved.kind !== "layer" || !display || display.resolved !== resolved) {
      const note = document.createElement("p");
      note.className = "galaxy-slab-note";
      if (slabs.length) {
        const lo = slabs[0];
        const hi = slabs[slabs.length - 1];
        note.textContent = "Showing " + (lo === hi ? noun.toLowerCase() + " " + lo : noun.toLowerCase() + "s " + lo + " to " + hi)
          + (resolved.kind === "layer" ? "." : "; pick an arc on the map.");
      }
      box.appendChild(note);
      return;
    }
    const data = display.data;
    const list = document.createElement("ul");
    resolved.options.map(function (option, index) { return index; }).reverse().forEach(function (index) {
      const option = resolved.options[index];
      const sum = sumOf(option.blocks, data);
      const pick = option.pick;
      const item = document.createElement("li");
      const button = document.createElement("button");
      button.type = "button";
      button.className = "galaxy-slab-row";
      button.dataset.lo = String(pick.lo);
      button.dataset.hi = String(pick.hi);
      const name = document.createElement("span");
      name.textContent = S.pickLabel(option.pick, stage.at, resolved.view);
      const bar = document.createElement("span");
      bar.className = "galaxy-slab-bar";
      const fill = document.createElement("span");
      fill.style.width = sum.total > 0 ? Math.max(sum.generated > 0 ? 2 : 0, (100 * sum.generated) / sum.total).toFixed(1) + "%" : "0%";
      bar.appendChild(fill);
      const count = document.createElement("span");
      count.className = "galaxy-slab-count";
      count.textContent = S.formatCount(sum.generated) + (sum.total > 0 ? " / " + S.formatCount(sum.total) : "");
      button.appendChild(name);
      button.appendChild(bar);
      button.appendChild(count);
      button.setAttribute("aria-label", layerText(pick) + ", " + S.formatInt(sum.generated)
        + (getOutline().shapeless ? "" : " of " + S.formatInt(sum.total)) + " sectors generated");
      button.disabled = generatedOnly && !(sum.generated > 0);
      button.addEventListener("mouseenter", function () { if (!animation) setHover({ layer: pick }); });
      button.addEventListener("mouseleave", function () { if (hover && !hover.sticky) setHover(null); });
      button.addEventListener("focus", function () { if (!animation) setHover({ layer: pick, sticky: true }); });
      button.addEventListener("click", function () { go({ at: stage.at, picks: stage.picks.concat([pick]) }); });
      item.appendChild(button);
      list.appendChild(item);
    });
    box.appendChild(list);
  }

  // Marks the strip row holding `layer` ({lo, hi}, or null for none).
  function markStripRow(layer) {
    if (!els.slabs) return;
    els.slabs.querySelectorAll(".galaxy-slab-row").forEach(function (row) {
      row.classList.toggle("is-hovered", !!layer && Number(row.dataset.lo) <= layer.lo && Number(row.dataset.hi) >= layer.hi);
    });
  }

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
    host.setBlockSize(r.view.blocks.length ? r.view.blocks[0].m : 243);
    notice(problem || "");
  }

  function onPopState(event) {
    if (!active) return;
    const state = event.state && event.state.galaxyStage ? event.state : null;
    mapIndex = state && state.mapIndex != null ? state.mapIndex : 0;
    maxIndex = Math.max(maxIndex, mapIndex);
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
    onPointerDown: onPointerDown,
    onPointerMove: onPointerMove,
    onPointerUp: onPointerUp,
    onPointerLeave: onPointerLeave,
    onKey: onKey,
    home: home,
    up: up,
    travel: travel,
    setGeneratedOnly: setGeneratedOnly,
    stage: function () { return stage; },
    go: go,
    locate: locate,
  };
}
