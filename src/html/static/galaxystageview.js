// html/static/galaxystageview.js
//
// The Galaxy Map's drill-down (docs/design/interface/galaxy-drilldown-navigation.md,
// sections 4, 5, 8.1 and 10-11, as MAP.150 and MAP.151 changed it), drawn
// and driven: the whole galaxy in 3D with every block of it to click, then,
// for each click, the window of blocks one size down centred on the block
// clicked, ... down to a window of sectors. galaxystages.js has the rules;
// this file has the scene, the camera moves, the breadcrumb, the tooltip,
// keys, touch, the stage URLs and the map's own Back and Forward.
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
//   (and the nodes, scale bar and wedge lines with it);
// - fetchStage(query): a Promise of GET /galaxy/stage's JSON for a
//   container ("?at=m.shell.band.wedge", or "" for the galaxy);
// - fetchWindow(query): a Promise of GET /galaxy/window's JSON for a window
//   ("?at=m.shell.band.wedge&g=9");
// - showBlockInfo(info), showPlacedInfo(entry), showCellInfo(cell),
//   showHint(text): the info panel (showBlockInfo's info: {block, total,
//   generated, hint, enter});
// - picker: the map's picker (mappick.js; a new one when missing), whose
//   other layers (phenomena, clouds) take a click or hover inside a
//   container before the blocks: their tooltip(entry), hover(entry) and
//   select(entry) are called; the drill-down adds its choices to it;
// - clearSelection(): a block was picked, so a phenomenon's ring goes;
// - restingHint: what a pinned sector's panel says before a pick, back when it is cleared;
// - canGenerate, courseSectors, sectorUrl(id), locate(name);
// - fetchSectorScene(id): a Promise of GET /sector/<id>/scene's JSON, the
//   sector opened in place (MAP.66; galaxysector.js);
// - fetchUnchartedScene({shell, band, wedge}) (optional): a Promise of the
//   scene JSON of a cell nothing was generated in, which opens in place
//   with an "Uncharted" mark when the scatters left objects in it (MAP.162);
// - fetchSystemScene(href): a Promise of GET /system/<id>/scene's JSON, a
//   system of the open sector opened in place (MAP.125; galaxysystem.js); lightBackground,
//   pixelRatio(), showInfo(spec), hideCell(bounds or null), preview(facts, entry) (the hover preview of a thing of the open
//   sector: `facts()` is its card's spec, null with no entry when the pointer leaves it; MAP.348);
// - els: {crumbs, tooltip, notice, address, matches, controls}
//   (any may be missing).

const VERSION_QUERY = new URL(import.meta.url).search;

const S = await import(`./galaxystages.js${VERSION_QUERY}`);
const P = await import(`./galaxyprisms.js${VERSION_QUERY}`);
const MC = await import(`./mapcontrol.js${VERSION_QUERY}`);
const { worldUnitsPerPixel, nearestOnScreen } = await import(`./mapcore.js${VERSION_QUERY}`);
const B = await import(`./bookmarks.js${VERSION_QUERY}`);
const { createSelection } = await import(`./picker.js${VERSION_QUERY}`);
const { createBreadcrumb } = await import(`./breadcrumb.js${VERSION_QUERY}`);
const { createPicker, createTooltip } = await import(`./mappick.js${VERSION_QUERY}`);
const { createSectorStage, DEFAULT_HIDDEN_KINDS } = await import(`./galaxysector.js${VERSION_QUERY}`);
const { createLayerRegistry } = await import(`./maplayers.js${VERSION_QUERY}`);
const { KINDS } = await import(`./sectorscene.js${VERSION_QUERY}`);
const { createSystemStage } = await import(`./galaxysystem.js${VERSION_QUERY}`);
const { nearestInDirection } = await import(`./systemstep.js${VERSION_QUERY}`);
const CROSS = await import(`./crossstars.js${VERSION_QUERY}`);
const NF = await import(`./nearfield.js${VERSION_QUERY}`);
const GO = await import(`./gridopacity.js${VERSION_QUERY}`);
const SR = await import(`./sectorreach.js${VERSION_QUERY}`);
const SYS = await import(`./systemreach.js${VERSION_QUERY}`);
const FIT = await import(`./flytofit.js${VERSION_QUERY}`);
const RL = await import(`./regionlanding.js${VERSION_QUERY}`);
const FPLACE = await import(`./findplace.js${VERSION_QUERY}`);
const { createFindBox } = await import(`./findplacebox.js${VERSION_QUERY}`);
const { centreOf, offsetsFrom } = await import(`./camerarelative.js${VERSION_QUERY}`);

// The other choices while one is hovered (MAP.18).
const OTHER_FADE = 0.25;
const KEY_MOVE = 0.15; // MAP.200: a step of the center by key, a share of the map's width
const REST_BLOCK_EDGES = 0.1; // MAP.202: the share of the lines between blocks drawn at rest
const HOVER_SHADE_OPACITY = 0.16; // MAP.198: the shade over the block a click would open
// A new stage's blocks fade in over the last part of a flight.
const FADE_IN_MS = 200;
// The block of the sector opened in place, round its own scene.
const ENTERED_BLOCK_FADE = 0.06;
// The other sectors' blocks while one is opened: only the sector being viewed stands out.
const ENTERED_OTHERS_FADE = 0.12;
// The sector opened fills the window as closely as its corners allow (the usual fit leaves a tenth).
const SECTOR_FIT_MARGIN = 1.02;
const DRAG_CLICK_PX = 6;
// What the + and - buttons of a sector page's map zoom by (MAP.68).
const BUTTON_ZOOM = 1.25;
// Each zoom step flies the camera to a preset for what it shows (MAP.97,
// Boss 2026-10-02: "when it zooms to a block it moves as isometric ...
// The only exception is the galaxy strtas out top-down so we can see the
// spiral arms."): the whole galaxy straight down, a window (a ball of cells
// on the shells, MAP.205) at the isometric slant, so its cells show side by
// side and one above another. Cells touch: no space between them (Boss,
// 2026-10-01).
export const ISO_TILT = Math.atan(Math.SQRT2);
// Every view can be turned, moved and zoomed (Boss, 2026-10-01; the
// whole galaxy too since MAP.85): drag turns it, right-drag or
// Shift-drag moves it, the wheel or a pinch zooms. It turns any way by
// any amount, through edge-on and under the plane, trackball style
// (MAP.96: "rotate the contents of the galaxy map any direction any
// amount"); zoom runs from MIN_ZOOM to MAX_ZOOM times the stage's own fit
// (the whole galaxy: GALAXY_MIN_ZOOM, about twice as close, out to its
// fit and no further, MAP.58); the view's middle can go anywhere (MAP.263).
const ROTATE_PER_PX = (0.4 * Math.PI) / 180;
// Shift and an arrow key turn the view this far.
const KEY_TURN = (5 * Math.PI) / 180;
export const MIN_ZOOM = 1 / 8;
export const MAX_ZOOM = 2.5;
export const GALAXY_MIN_ZOOM = 0.5;
// A sector page's map (pinned) zooms out no further than this times the sector's own fit: the sector, not its surroundings.
export const PINNED_MAX_ZOOM = 1.1;
// An opened system zooms from its fit down to a moon, nine orders closer.
export const SYSTEM_MIN_ZOOM = 1e-9;
// MAP.333: the nearest the camera may be asked to come to its target by any move, and the farthest, the whole galaxy in view (the same
// 2.4 galaxy radii the whole-galaxy fit uses). MAP.338: 30 km, the camera held as the target plus an exact offset (galaxymap3d.js draws an
// opened system from the eye's own frame once the camera is that close), where it was 3,000 km, about a moon framed (Boss, 6:19 PM EDT 2026-10-28).
export const GLOBAL_FLOOR_PC = 1e-12;
// A sector selected and left stays in the address while the camera is within this many sector edges of its target.
const SELECTION_REACH_EDGES = 81;
export const GALAXY_FAR_RADII = 2.4;
// MAP.334: the nearest the camera comes outside an opened system, until the scene has the precision for finer (MAP.338): the single-precision
// positions of the sector and the windows are good to about a thousandth of a parsec a few thousand parsecs out.
export const FIELD_FLOOR_PC = 1e-3;
const WHEEL_ZOOM_PER_PX = 0.0025;
// MAP.225: how close to the thing under the pointer a wheel step in leaves the eye, as a share of the view's distance.
const ZOOM_CLEARANCE = 0.02;
// A star counts as the one zoomed on when it is within this share of the view's distance from its middle.
const NEAREST_STAR_REACH = 0.35;
// MAP.292: how far past the end of its range (a share of the nearest zoom) the wheel carries the view while the next window arrives.
const ARRIVING_OVERSHOOT = 0.25;

// MAP.245: how near a star is to the pointer, in pixels, for the wheel to zoom toward it.
const STAR_ANCHOR_PX = 14;

// MAP.247: the wheel's point on a cursor ray that meets nothing pickable is where the stars' column along the ray
// passes COLUMN_SHARE of the view's distance in local-average density (1 is the average): in space of that
// density half way to the focus, in a dense region (a bulge, an arm) nearer, in a sparse one not at all (the
// point of the focus plane is taken). The ray is read out to COLUMN_REACH times the eye's distance to the focus,
// and a point is never nearer than COLUMN_NEAREST of the way to the focus, so a step in still moves the eye.
const COLUMN_SHARE = 0.5;
const COLUMN_REACH = 2;
const COLUMN_NEAREST = 0.25;
// MAP.301: a wheel step heads for the densest layer of drawn stars in front of the eye along the cursor ray when there are at least
// CROWD_MIN_STARS of them within CROWD_CONE (the tangent of the angle, about 24 pixels at the default lens) of it, in front of the
// thing the wheel would otherwise head for and no nearer than twice the clearance, so the approach slows and stops in a cluster instead of flying through it.
const CROWD_CONE = 0.025;
const CROWD_MIN_STARS = 12;
const CROSS_PX = 12; // MAP.247: a pointer this near the middle of the map is on the cross

// The zoom policies (mapcontrol.js): a short range around the fit, and
// on the whole galaxy only closer than its fit.
const FREE_VIEW_ZOOM = MC.zoomPolicy(MC.ZOOM_RANGE, MIN_ZOOM, MAX_ZOOM);
const GALAXY_ZOOM = MC.zoomPolicy(MC.ZOOM_RANGE, GALAXY_MIN_ZOOM, 1);
const SYSTEM_ZOOM = MC.zoomPolicy(MC.ZOOM_RANGE, SYSTEM_MIN_ZOOM, MAX_ZOOM);
// MAP.152: in a window of sectors the camera may come this many sector edges from the sector at the
// middle (the window's own fit and 1/8 of it stop about 1.4 edges away, short of where a sector's
// own scene is whole, sectorreach.js FULL_EDGES).
const SECTOR_NEAREST_EDGES = 0.5;
// A stage of more choices than this is drawn as one mesh: a hovered choice
// is outlined, and the others don't fade (a mesh and a draw call each would
// be too many).
const FADE_GROUP_LIMIT = 800;
// The outlines of charted blocks with "Charted only" on (MAP.111).
const CHARTED_LINE_OPACITY = 0.9;

// The free view's limits: the camera's distance from MIN_ZOOM to MAX_ZOOM times the stage's
// fit ...
export function clampZoom(dist, fitDist) {
  return MC.clampDistance(FREE_VIEW_ZOOM, dist, fitDist);
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
  let stage = S.galaxyStage();
  let resolved = null;
  let selectedSector = null;
  // MAP.265: the sector a click on a coarser block read under it (closed form): outlined and described once the next
  // window is in, but no part of the selection (the URL, the breadcrumb, the bookmark); a click on the
  // sector itself, once it is a choice, selects it as ever.
  let pickedSector = null;
  // The sector opened in place: {sector: {shell, band, wedge}, id, center,
  // halfEdge, fitPoints}, while galaxysector.js has its scene up.
  // A sector page's map (MAP.68): the one sector, shell/band/wedge, it is
  // locked to; it has no steps, history or URL of its own.
  const pinned = host.pinned || null;
  let entered = null;
  let enterToken = 0;
  // The system opened in place inside the open sector (MAP.125): {id, href,
  // name, center, radiusPc, fitPoints}.
  let enteredSystem = null;
  let enterSystemToken = 0;
  // A system to open once its sector has (a reload or Back to its URL).
  let pendingSystem = null;
  // MAP.329: how to open that system (history entry and scale) when a step to it crossed a sector's edge.
  let pendingSystemOptions = null;
  // The body of that system the URL named (?object=planet:12).
  let pendingObject = null;
  // The body of the open system the camera rides along with (a moving target).
  let followBody = null;
  // A sector to open once the stage it sits in has arrived (a reload or
  // Back to an opened sector's URL).
  let pendingOpen = null;
  let display = null;
  let leaving = [];
  let animation = null;
  // {option (index into display.options), sticky}
  let hover = null;
  let chartedOnly = false;
  // MAP.265: the blocks (translucent prisms) are drawn only when the reader asks: the cell grid tool, a "Color by" fill
  // or "Charted only". They are always there to be picked; the one cell the reader sees otherwise is the one hovered
  // or selected.
  let cellGrid = false;
  // What the blocks' fill shows (MAP.131): a key of galaxyblocks' COLOR_MODES.
  let colorBy = "default";
  // Choosing a NAV end (setNeedGenerated): only a choice holding
  // something generated can be taken, since endpoints live in generated
  // sectors. "Charted only" itself never blocks picking (MAP.112).
  let needGenerated = false;
  let view = null;
  let goToken = 0;
  let arriving = false;
  // The map's own history (Back and Forward buttons): this entry's index,
  // and the furthest one known.
  let mapIndex = 0;
  let maxIndex = 0;

  // --- The outline and the stage API's counts ----------------------------

  function getOutline() {
    if (!outline) {
      outline = S.galaxyOutline(edgePc, host.shape, host.galaxyRadius);
      galaxyFarPc = 0;
      outline.shapeless = !host.shape;
    }
    return outline;
  }

  function resolve(s) {
    return S.resolveStage(s, getOutline(), edgePc);
  }

  function containerKey(win) {
    return win ? "w" + S.formatDrillKey(win.centre) + "g" + win.g : "galaxy";
  }

  // The stage API's answer for window `win` (null: the galaxy), as
  // {generated: Map(key -> count), sectors: Map("shell/band/wedge" ->
  // sector)}; null while it loads (and on an error, which leaves counts
  // at 0).
  function dataFor(win) {
    const cached = dataCache.get(containerKey(win));
    return cached && cached.ready ? cached.value : null;
  }

  function loadData(win) {
    const key = containerKey(win);
    let entry = dataCache.get(key);
    if (!entry) {
      entry = { ready: false, value: null };
      const fetched = win
        ? host.fetchWindow("?at=" + S.formatDrillKey(win.centre) + "&g=" + win.g)
        : host.fetchStage("");
      entry.promise = fetched.then(function (payload) {
        const generated = new Map();
        const stats = new Map();
        (payload.children || []).forEach(function (child) {
          generated.set(cellKey(child), child.generated);
          stats.set(cellKey(child), child.stats || null);
        });
        const sectors = new Map();
        (payload.sectors || []).forEach(function (sector) {
          sectors.set(cellKey(sector), sector);
        });
        entry.value = { generated: generated, stats: stats, sectors: sectors, sectorGenerated: new Map() };
        payload.sectors && payload.sectors.forEach(function (sector) {
          entry.value.sectorGenerated.set(cellKey(sector), 1);
        });
        entry.ready = true;
        return entry.value;
      }).then(null, function () {
        entry.value = { generated: new Map(), stats: new Map(), sectors: new Map(), sectorGenerated: new Map(), failed: true };
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
      loadData(S.holder(stage)).then(function () {
        if (display && !animation) rebuildDisplay();
      });
    }
  }

  // A cell's key among the data's: "shell/band/wedge" (the level is the container's).
  function cellKey(cell) {
    return cell.shell + "/" + cell.band + "/" + cell.wedge;
  }

  function generatedOf(block, data) {
    if (!data) return 0;
    const key = cellKey(block);
    if (block.m === 1) return data.sectorGenerated.get(key) || 0;
    return data.generated.get(key) || 0;
  }

  // DB.14, MAP.128: what a block's generated sectors hold, for its color
  // and opacity ({systems, expected_systems, stars, mean_age_gy,
  // luminosity_sol} from the stage API), or null.
  function statsOf(block, data) {
    if (!data) return null;
    const key = cellKey(block);
    if (block.m === 1) {
      const sector = data.sectors.get(key);
      return sector ? sector.stats || null : null;
    }
    return data.stats.get(key) || null;
  }

  // The choices a stage offers: [{pick, blocks, a0, a1}], one per block.
  // Without an outline (no shape yet), only blocks holding generated sectors.
  function choicesOf(r, data) {
    let options = r.options;
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

  // Whether the view can be turned and moved: every stage, since the
  // whole galaxy turns too (MAP.85).
  function isFree(r) {
    return !!(r && r.stage);
  }

  // The view's zoom policy: a short range round the fit; the whole galaxy
  // only zooms in from its own.
  function zoomPolicyFor(r, at) {
    if (enteredSystem) return SYSTEM_ZOOM;
    if (r && r.view && isWholeGalaxy(r)) return GALAXY_ZOOM;
    if (pinned) return MC.zoomPolicy(MC.ZOOM_RANGE, MIN_ZOOM, PINNED_MAX_ZOOM);
    if (entered && !pinned && view && view.fit && host.fetchSystemScene) {
      // MAP.152: in an opened sector the camera may come right up to a star near the middle of the view, to
      // where its system's disc is whole, but no closer to the middle than keeps that star in reach.
      const fit = fitDistance(entered.fitPoints, view.fit.target, view.quat);
      const near = nearestStarTo(at || view.target);
      if (fit > 0 && near) {
        const radius = systemRadii.get(near.entry.href) || SYS.NOMINAL_RADIUS_PC;
        const discPc = SYS.distanceForPixels(SYS.FULL_PX * 1.5, radius, canvasEl.clientHeight || 600, camera.fov);
        const nearest = Math.min(MIN_ZOOM, Math.max(discPc, near.distance / NEAREST_STAR_REACH) / fit);
        return MC.zoomPolicy(MC.ZOOM_RANGE, nearest, MAX_ZOOM);
      }
    }
    if (!entered && !pinned && view && view.fit && isSectorView(r)) {
      const fit = fitDistance(fitPointsOf(r), view.fit.target, view.quat);
      if (fit > 0) return MC.zoomPolicy(MC.ZOOM_RANGE, Math.min(MIN_ZOOM, (SECTOR_NEAREST_EDGES * edgePc) / fit), MAX_ZOOM);
    }
    return FREE_VIEW_ZOOM;
  }

  function isWholeGalaxy(r) {
    return !r.stage.win;
  }

  // --- Drawing a stage ---------------------------------------------------------

  // How many sectors `block` can hold: 1 for a sector, else the count under the outline (null without a shape).
  function totalOf(block) {
    return block.m === 1 ? 1 : S.blockTotal(block, getOutline());
  }

  // One stage's meshes: a group per choice (so the hovered one stays and
  // the others fade), with each cell's choice behind the meshes for
  // picking.
  function buildDisplay(r) {
    const data = dataFor(S.holder(r.stage));
    const options = choicesOf(r, data);
    const root = new THREE.Group();
    const groups = [];
    const eye = camera.position.toArray();
    // Every choice's cells, ranked against each other (MAP.131) so one
    // color means the same thing across the stage, and the legend can say it.
    const cellsOf = options.map(function (option, index) {
      return option.blocks.map(function (block) {
        const b = block.bounds;
        const filled = generatedOf(block, data);
        return {
          shell: block.shell, band: block.band, wedge: block.wedge,
          r0: b.r0, r1: b.r1, c0: b.c0, c1: b.c1, t0: b.t0, t1: b.t1,
          filled: filled, total: filled > 0 ? totalOf(block) : 0, stats: statsOf(block, data), block: block,
          option: index,
        };
      });
    });
    const stageRanges = host.blockScene.rangesOf([].concat.apply([], cellsOf));
    if (host.showLegend) host.showLegend(host.blockScene.legendOf(colorBy, stageRanges));
    // One group per choice, so the hovered one stays and the others fade; a
    // stage of very many choices is one group, its hovered choice outlined.
    const single = options.length > FADE_GROUP_LIMIT;
    const partCells = single ? [[].concat.apply([], cellsOf)] : cellsOf;
    partCells.forEach(function (cells, index) {
      const dim = chartedOnly ? function (cell) { return !(cell.filled > 0); } : null;
      const built = host.blockScene.buildCells(cells, eye, dim, colorBy, stageRanges);
      const group = new THREE.Group();
      const meshes = [];
      // The whole galaxy shows no sector or block lines (MAP.85); a window
      // shows the lines between its blocks.
      const gridEdges = isWholeGalaxy(r) ? 0 : REST_BLOCK_EDGES;
      [built.solid, built.glass].forEach(function (part, n) {
        if (!part.vertexCount) return;
        const mesh = host.makeBlockMesh(part, n === 1);
        if (mesh.material.uniforms && mesh.material.uniforms.gridEdges) mesh.material.uniforms.gridEdges.value = gridEdges;
        mesh.userData.cells = built.cells[n];
        mesh.visible = cellsShown();
        mesh.renderOrder = n;
        group.add(mesh);
        meshes.push(mesh);
      });
      root.add(group);
      groups.push({ option: single ? null : options[index], meshes: meshes, fade: 1, restEdges: gridEdges });
    });
    const charted = chartedOnly ? chartedLines(options, data) : null;
    if (charted) root.add(charted);
    host.scene.add(root);
    root.updateMatrixWorld(true);
    const gridStats = GO.cellStats([].concat.apply([], cellsOf));
    return { resolved: r, root: root, options: options, groups: groups, fade: 1, data: data, charted: charted, gridStats: gridStats };
  }

  // MAP.111: with "Charted only" on, every block holding generated
  // sectors (one sector, once the view is down to sectors) outlined, so
  // where the charted space lies in the wedge is plain to see.
  function chartedLines(options, data) {
    const points = [];
    options.forEach(function (option) {
      option.blocks.forEach(function (block) {
        if (!(generatedOf(block, data) > 0)) return;
        Array.prototype.push.apply(points, S.cellEdgePoints(block.bounds));
      });
    });
    if (!points.length) return null;
    const line = lineObject(
      points,
      new THREE.LineBasicMaterial({ color: accent, transparent: true, opacity: CHARTED_LINE_OPACITY, depthTest: false, depthWrite: false }),
      "charted lines",
    );
    line.renderOrder = 4;
    line.frustumCulled = false;
    return line;
  }

  // A line-segment geometry from a flat [x, y, z, ...] list of end points, as offsets from `origin` (MAP.240).
  function lineGeometry(points, origin) {
    const geometry = new THREE.BufferGeometry();
    geometry.setAttribute("position", new THREE.BufferAttribute(offsetsFrom(points, origin), 3));
    return geometry;
  }

  // The line segments through `points`, an object at the middle of them with offsets from it (MAP.240).
  function lineObject(points, material, layer) {
    const origin = centreOf(points);
    const line = new THREE.LineSegments(lineGeometry(points, origin), material);
    line.position.set(origin[0], origin[1], origin[2]);
    line.userData.relativeLayer = layer;
    return line;
  }

  function disposeDisplay(d) {
    if (!d) return;
    host.scene.remove(d.root);
    if (d.charted) {
      d.charted.geometry.dispose();
      d.charted.material.dispose();
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
    if (d.charted) d.charted.material.opacity = CHARTED_LINE_OPACITY * fade;
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
    if (view) updateGridScale();
    disposeDisplay(old);
    applyHover();
  }

  // --- Framing -----------------------------------------------------------------

  // The points a stage's view is fitted round, worked out once per stage.
  const fitPointsCache = new WeakMap();
  function fitPointsOf(r) {
    let points = fitPointsCache.get(r);
    if (!points) {
      points = S.boundsPoints(r.view.blocks);
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
    const margin = entered && points === entered.fitPoints ? SECTOR_FIT_MARGIN : S.FIT_MARGIN;
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
  // galaxy, the isometric slant for a window.
  function presetTilt(r) {
    return isWholeGalaxy(r) ? 0 : ISO_TILT;
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
    const extent = S.extentOf(blocks);
    // A window looks from outside its middle, the core's bearing up the screen.
    const theta = isWholeGalaxy(r) ? -Math.PI / 2 : Math.atan2(extent.center[1], extent.center[0]) + Math.PI;
    const quat = presetQuat(theta, presetTilt(r));
    const fit = centeredFit(fitPointsOf(r), extent.center, quat);
    return { target: fit.target, dist: fit.dist, quat: quat };
  }

  // The camera for the sector opened in place: its cell fitted at the
  // isometric slant, the core's bearing up the screen like a block's.
  function cameraForSector(open) {
    const quat = presetQuat(Math.atan2(open.center[1], open.center[0]) + Math.PI, ISO_TILT);
    const fit = centeredFit(open.fitPoints, open.center.slice(), quat);
    return { target: fit.target, dist: fit.dist, quat: quat };
  }

  // The camera for the system opened in place: its extent fitted at the
  // isometric slant, like the sector it sits in.
  function cameraForSystem(open) {
    const quat = presetQuat(Math.atan2(open.center[1], open.center[0]) + Math.PI, ISO_TILT);
    const fit = centeredFit(open.fitPoints, open.center.slice(), quat);
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

  // Where a point of the scene (x, y, z in parsecs, as the map's own
  // frame) falls on the page, in client pixels, or null behind the camera.
  // Read by the browser tests, which click on the charted space.
  canvasEl.galaxyScreenOf = function (x, y, z) {
    camera.updateMatrixWorld();
    const v = new THREE.Vector3(x, y, z).project(camera);
    if (v.z > 1) return null;
    const rect = canvasEl.getBoundingClientRect();
    return { x: rect.left + ((v.x + 1) / 2) * rect.width, y: rect.top + ((1 - v.y) / 2) * rect.height };
  };

  // MAP.245: the point the wheel zooms toward at a screen point (`anchorUnder`), [x, y, z] in the map's frame, or
  // null. Read by the browser tests.
  canvasEl.galaxyAnchorUnder = function (clientX, clientY) { return anchorUnder(clientX, clientY); };

  // MAP.243: whether an object of the map (anything with x, y, z in parsecs) lies beyond the galaxy's outline, in the
  // halo or a satellite galaxy, where no block of the whole-galaxy view holds it. Worked out once per object.
  function outsideOutline(item) {
    if (item.beyondOutline === undefined) {
      const at = P.sectorAddressAt(item.x, item.y, item.z, edgePc);
      item.beyondOutline = !P.sectorInOutline(at.shell, at.band, at.wedge, getOutline());
    }
    return item.beyondOutline;
  }

  // For the browser tests: whether a point of the map (x, y, z in parsecs) lies beyond the galaxy's outline.
  canvasEl.galaxyOutsideOutline = function (x, y, z) { return outsideOutline({ x: x, y: y, z: z }); };

  // Whether the map has caught up with its stage: what is drawn is the
  // stage the breadcrumb names and no flight is running. The stage changes
  // at once, its drawing when the container's data has loaded. Read by
  // the browser tests, which wait on it instead of a clock (TEST.96,
  // TEST.98, TEST.99, TEST.100).
  canvasEl.galaxyReady = function () {
    return !!display && !!resolved && display.resolved === resolved && !animation && leaving.length === 0;
  };

  // The system opened in place, if any: {id, name, mode, refs (its bodies),
  // flying, following}. Read by the browser tests.
  canvasEl.galaxySystem = function () {
    if (!enteredSystem) return null;
    return {
      id: enteredSystem.id, name: enteredSystem.name, mode: systemStage.mode(), refs: systemStage.entries().map(function (e) { return e.ref; }),
      spins: systemStage.spins, flying: !!animation, following: followBody, trailOpacity: systemStage.trailOpacity, trailRadius: systemStage.trailRadius, distance: view ? view.dist : null, selected: systemStage.selected(),
      // Where a body is on screen, in client pixels (null: behind the camera).
      screen: function (ref) {
        const entry = systemStage.entry(ref);
        if (!entry) return null;
        camera.updateMatrixWorld();
        const v = new THREE.Vector3(entry.x, entry.y, entry.z).project(camera);
        if (v.z > 1) return null;
        const rect = canvasEl.getBoundingClientRect();
        return { x: rect.left + ((v.x + 1) / 2) * rect.width, y: rect.top + ((1 - v.y) / 2) * rect.height };
      },
    };
  };

  // For the browser tests (TEST.448): what the stage drew for each charted block, {shell, band, wedge, m, filled, habitable,
  // noFill, opacity}: `noFill` is a block the chosen colour mode leaves unfilled (no habitable worlds under "Habitable worlds"),
  // which a screenshot's few faint pixels cannot tell from a filled one.
  canvasEl.galaxyBlockLooks = function () {
    const looks = [];
    if (!display) return looks;
    display.groups.forEach(function (group) {
      group.meshes.forEach(function (mesh) {
        (mesh.userData.cells || []).forEach(function (cell) {
          if (!(cell.filled > 0) || !cell.block) return;
          looks.push({
            shell: cell.block.shell, band: cell.block.band, wedge: cell.block.wedge, m: cell.block.m, filled: cell.filled,
            habitable: cell.stats ? cell.stats.habitable_worlds : null, noFill: !!cell.noFill, opacity: cell.opacity,
          });
        });
      });
    });
    return looks;
  };

  // Which lines the stage draws: {blockEdges} (the blocks' own edges on or
  // off), {blocksDrawn} (how many block meshes are drawn: none unless the cell grid, a fill or "Charted only" is on,
  // MAP.265) and {chartedLines} (how many line pieces outline the charted
  // blocks, MAP.111). Read by the browser tests.
  canvasEl.galaxyLines = function () {
    if (!display) return null;
    let blockEdges = 0;
    let blocksDrawn = 0;
    display.groups.forEach(function (group) {
      group.meshes.forEach(function (mesh) {
        if (mesh.visible) blocksDrawn += 1;
        if (mesh.material.uniforms && mesh.material.uniforms.gridEdges) blockEdges = Math.max(blockEdges, mesh.material.uniforms.gridEdges.value);
      });
    });
    const chartedLines = display.charted ? display.charted.geometry.getAttribute("position").count / 2 : 0;
    // MAP.198: the hovered block's outline: its line pieces, the heights
    // they reach, and the shade over its faces.
    const hovered = { lines: 0, zMin: null, zMax: null, shade: 0 };
    if (outlineGroup) {
      outlineGroup.children.forEach(function (part) {
        if (part.userData.hoverShade) {
          hovered.shade = part.material.opacity;
          return;
        }
        const position = part.geometry.getAttribute("position");
        hovered.lines += position.count / 2;
        for (let i = 0; i < position.count; i++) {
          const z = position.getZ(i) + part.position.z;
          hovered.zMin = hovered.zMin === null ? z : Math.min(hovered.zMin, z);
          hovered.zMax = hovered.zMax === null ? z : Math.max(hovered.zMax, z);
        }
      });
    }
    return { blockEdges: blockEdges, blocksDrawn: blocksDrawn, chartedLines: chartedLines, hovered: hovered };
  };

  // The camera's tilt from straight down, in degrees (0 top-down, 90
  // edge-on, 180 from under the plane), and whether a flight is running.
  // Read by the browser tests.
  // MAP.334: the range the camera may zoom over on the stage on show, in parsecs ({nearestPc, farthestPc}), and the nearest thing along the
  // ray through a screen point ({point, kind, distance, radius?}); read by the browser tests.
  canvasEl.galaxyZoomLimits = function () {
    const fit = fitNow();
    if (!view || !(fit > 0)) return null;
    const limits = zoomLimits();
    return { nearestPc: limits.nearest * fit, farthestPc: limits.farthest * fit };
  };
  canvasEl.galaxyNearestThing = function (x, y) { return nearestThingAlong(x, y); };

  canvasEl.galaxyCamera = function () {
    if (!view) return null;
    const back = new THREE.Vector3(0, 0, 1).applyQuaternion(view.quat);
    return {
      tilt: THREE.MathUtils.radToDeg(Math.acos(Math.max(-1, Math.min(1, back.z)))), flying: !!animation,
      zoom: view.zoom, target: view.target.slice(), dist: view.dist, eye: eyeNow(),
    };
  };

  // Keeps the view fitted round the stage as it turns or the map changes
  // size (MAP.78): the distance is the user's zoom times the fit at the
  // view's own turn.
  function reframe() {
    const fit = fitNow();
    if (fit > 0) view.dist = view.zoom * fit;
  }

  // The distance at which the stage on show (the system, the sector or the window) fits the view as it is turned: what `view.zoom` is a ratio of.
  // 0 when there is no stage to fit.
  function fitNow() {
    if (!view || !view.fit || !resolved || !resolved.view || !resolved.view.blocks.length) return 0;
    const around = enteredSystem ? enteredSystem.fitPoints : entered ? entered.fitPoints : fitPointsOf(resolved);
    return fitDistance(around, view.fit.target, view.quat);
  }

  // The block of the stage the camera is inside, if any (MAP.149): its
  // center, for the block shader to draw its walls from the inside.
  function containerCenter(eye) {
    if (!display || !display.options) return null;
    for (let n = 0; n < display.options.length; n++) {
      const blocks = display.options[n].blocks;
      for (let k = 0; k < blocks.length; k++) {
        if (NF.blockContains(blocks[k].bounds, eye[0], eye[1], eye[2])) return NF.blockCenter(blocks[k].bounds);
      }
    }
    return null;
  }

  // MAP.188: the share the blocks' faces and lines are drawn at so the cells
  // a ray crosses add up to no more than GRID_ACCUMULATION_CAP; worked out
  // from the stage's cells and the way the camera looks. Read by the tests.
  let gridScaleNow = 1;
  let gridReading = { crossings: 1, scale: 1, maxOpacity: 0, worst: 0 };

  function updateGridScale() {
    const stats = display ? display.gridStats : null;
    const forward = new THREE.Vector3(0, 0, -1).applyQuaternion(view.quat);
    const n = GO.crossings(stats, [forward.x, forward.y, forward.z]);
    const maxOpacity = stats ? stats.maxOpacity : 0;
    gridScaleNow = GO.gridScale(n, maxOpacity);
    gridReading = {
      crossings: n, scale: gridScaleNow, maxOpacity: maxOpacity, worst: GO.accumulated(n, maxOpacity * gridScaleNow),
    };
    if (host.setGridScale) host.setGridScale(gridScaleNow);
  }

  canvasEl.galaxyGridOpacity = function () {
    return Object.assign({ cap: GO.GRID_ACCUMULATION_CAP, distance: view ? view.dist : null }, gridReading);
  };

  function applyView() {
    updateGridScale();
    const back = new THREE.Vector3(0, 0, view.dist).applyQuaternion(view.quat);
    const eye = [view.target[0] + back.x, view.target[1] + back.y, view.target[2] + back.z];
    host.setCamera({
      target: view.target.slice(), dist: view.dist, quaternion: view.quat.toArray(),
      focusRadius: entered ? entered.halfEdge : 0, inside: containerCenter(eye),
    });
    refreshWhere();
    scheduleCameraAddress();
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

  // A stage's query with what a URL keeps beyond the stage: the NAV pick
  // (host.pickQuery(), "?pick=...") and the kinds of object left off the map
  // ("hide=nebula,roguePlanet", MAP.79).
  function withKept(query) {
    const extra = [];
    const pickQuery = host.pickQuery ? host.pickQuery() : "";
    if (pickQuery) extra.push(pickQuery.slice(1));
    return joinQuery(query, extra.concat(displayParts()));
  }

  // What the map shows beyond the stage (kinds, star classes and floor, and
  // the fill's color), as query parts; a bookmark of a view keeps them too.
  function displayParts() {
    const extra = [];
    layers.writeUrl(new URLSearchParams()).forEach(function (value, name) { extra.push(name + "=" + value); });
    if (sectorStage.markedKinds().length) extra.push("mark=" + sectorStage.markedKinds().join(","));
    if (sectorStage.hiddenClasses().length) extra.push("stars=" + sectorStage.hiddenClasses().join(","));
    if (sectorStage.minLuminosity() > 0) extra.push("lum=" + sectorStage.minLuminosity());
    if (colorBy !== "default") extra.push("color=" + colorBy);
    return extra;
  }

  function joinQuery(query, extra) {
    if (!extra.length) return query;
    return (query ? query + "&" : "?") + extra.join("&");
  }

  function displayQuery(query) {
    return joinQuery(query, displayParts());
  }

  // The color mode a URL's `color` names ("default" for none or an unknown one).
  function colorFromLocation() {
    const raw = new URLSearchParams(location.search).get("color");
    return raw && host.colorModes && host.colorModes.indexOf(raw) >= 0 ? raw : "default";
  }

  // Takes the URL's color mode up (the map's select follows, host.colorChanged).
  function setColorFromLocation() {
    colorBy = colorFromLocation();
    if (host.colorChanged) host.colorChanged(colorBy);
  }

  // MAP.217: the kinds of object the map shows are layers of the registry the maps share (maplayers.js): the
  // Menu's buttons and the URL's `hide` (kinds turned off) and `show` (kinds that begin hidden, rogue planets,
  // MAP.137, when turned on) come from it, and the sector stage follows what it says.
  const layers = createLayerRegistry(pinned ? "sector" : "galaxy", { forgetUnknown: true });
  let layersFollowing = false; // the registry is being brought up to the stage or the URL, so its changes are not echoed back
  KINDS.forEach(function (pair) {
    layers.register({
      id: pair[0], label: pair[1], modes: ["galaxy", "sector"], group: "Show on the map", data: { kind: pair[0] },
      defaultOn: DEFAULT_HIDDEN_KINDS.indexOf(pair[0]) < 0, urlKey: "hide", urlKeyOn: "show",
      available: function () { return kindCount(pair[0]) > 0; },
      title: function () { return "Show or hide the " + pair[1].toLowerCase() + " (" + kindCount(pair[0]) + ")"; },
    });
  });
  layers.onChange(function (id, on) {
    if (layersFollowing) return;
    sectorStage.setKindHidden(id, !on);
    writeHiddenToUrl();
  });

  function kindCount(kind) {
    const found = sectorStage.kinds().find(function (entry) { return entry.kind === kind; });
    return found ? found.count : 0;
  }

  // The registry follows the URL, then the stage follows the registry.
  // The registry takes the address again (for layers the map registered since), with nothing echoed back to it.
  function readLayersFromLocation() {
    layersFollowing = true;
    try {
      layers.readUrl(location.search);
    } finally {
      layersFollowing = false;
    }
  }

  function takeLayersFromLocation() {
    readLayersFromLocation();
    sectorStage.setHiddenKinds(KINDS.map(function (pair) { return pair[0]; }).filter(function (kind) { return !layers.isOn(kind); }));
  }

  // The registry follows the stage (a selection showing its kind).
  function takeLayersFromStage() {
    layersFollowing = true;
    try {
      KINDS.forEach(function (pair) { layers.setOn(pair[0], !sectorStage.kindHidden(pair[0])); });
    } finally {
      layersFollowing = false;
    }
  }

  // Shows or hides a kind of object on the open sector, and keeps the choice
  // in the URL (changed in place: it is no step of its own).
  function setKindHidden(kind, hidden) {
    layers.setOn(kind, !hidden);
  }

  // MAP.123: the star classes left off and the luminosity floor, from the URL.
  function applyStarFiltersFromLocation() {
    const params = new URLSearchParams(location.search);
    const classes = (params.get("stars") || "").split(",").filter(Boolean);
    sectorStage.setHiddenClasses(classes);
    sectorStage.setMarkedKinds((params.get("mark") || "").split(",").filter(Boolean));
    const lum = Number(params.get("lum"));
    sectorStage.setMinLuminosity(lum > 0 && isFinite(lum) ? lum : 0);
  }

  function setKindMarked(kind, marked) {
    sectorStage.setKindMarked(kind, marked);
    writeHiddenToUrl();
  }

  function setStarClassHidden(starClass, hidden) {
    sectorStage.setStarClassHidden(starClass, hidden);
    writeHiddenToUrl();
  }

  function setMinLuminosity(value) {
    sectorStage.setMinLuminosity(value);
    writeHiddenToUrl();
  }

  function writeHiddenToUrl() {
    const params = new URLSearchParams(location.search);
    layers.writeUrl(params);
    const marked = sectorStage.markedKinds();
    if (marked.length) params.set("mark", marked.join(","));
    else params.delete("mark");
    const classes = sectorStage.hiddenClasses();
    if (classes.length) params.set("stars", classes.join(","));
    else params.delete("stars");
    if (sectorStage.minLuminosity() > 0) params.set("lum", String(sectorStage.minLuminosity()));
    else params.delete("lum");
    const query = params.toString().replace(/%2C/gi, ",");
    history.replaceState(history.state, "", location.pathname + (query ? "?" + query : "") + location.hash);
  }

  // The body picked in the open system is the URL's ?object=, in place.
  function writeObjectToUrl(ref) {
    if (!enteredSystem) return;
    const params = new URLSearchParams(location.search);
    if (ref) params.set("object", ref);
    else params.delete("object");
    const query = params.toString().replace(/%2C/gi, ",").replace(/%3A/gi, ":");
    history.replaceState(history.state, "", location.pathname + (query ? "?" + query : "") + location.hash);
  }

  // The NAV pick moved on (navpick.js, NAV.29): the URL carries the new one,
  // in place, so the view stays where it is and a reload keeps the pick.
  function pickChanged() {
    const params = new URLSearchParams(location.search);
    ["pick", "from", "to"].forEach(function (name) { params.delete(name); });
    const pickQuery = host.pickQuery ? host.pickQuery() : "";
    new URLSearchParams(pickQuery).forEach(function (value, name) { params.set(name, value); });
    const query = params.toString().replace(/%2C/gi, ",").replace(/%3A/gi, ":");
    history.replaceState(history.state, "", location.pathname + (query ? "?" + query : "") + location.hash);
  }

  // A kind shown or hidden by the sector stage itself (a selection of
  // something hidden shows its kind): the URL and the buttons follow.
  function kindsChanged() {
    takeLayersFromStage();
    writeHiddenToUrl();
    if (host.kindsChanged) host.kindsChanged();
  }

  // A new entry in the browser history (the map's Back and Forward).
  function pushEntry(query) {
    if (pinned) return;
    mapIndex += 1;
    maxIndex = mapIndex;
    history.pushState({ galaxyStage: true, mapIndex: mapIndex, maxIndex: maxIndex }, "", location.pathname + query + location.hash);
  }

  // MAP.337: goes to stage `next` by moving the camera: it flies to the stage's own camera (or `options.camera`) and the
  // stage the crumbs, Up and the steps work from becomes `next`; the display drawn stays as it is, and the scenes come up from
  // where the camera is. push: record it in the browser history (default yes).
  // MAP.336: the "densest" landing a region's window is to go on to once it has arrived (goToRegion).
  let landingAfter = null;
  let landingToken = 0;
  function flyToStage(next, options) {
    options = options || {};
    if (pinned || !view) return;
    landingAfter = null;
    landingToken += 1;
    if (animation) finishAnimation();
    const r = resolve(next);
    if (r.problem) {
      notice(r.problem);
      return;
    }
    notice("");
    if (entered) closeSector();
    else {
      // MAP.352: a sector or system still loading for the view being left is not wanted in the one taken up now.
      enterToken += 1;
      enterSystemToken += 1;
      entering = false;
    }
    if (!options.keepSector) selectedSector = null;
    if (!options.keepPicked) pickedSector = null;
    const token = ++goToken;
    arriving = true;
    const to = options.camera ? asCamera(options.camera) : cameraFor(r);
    const query = withKept(options.query || cameraQueryFor(r, to));
    if (options.push !== false && (!S.sameStage(stage, r.stage) || withoutCam(query) !== withoutCam(location.search))) {
      pushEntry(query);
    }
    stage = r.stage;
    renderCrumbs();
    renderTravel();
    hover = null;
    showTooltip("", 0, 0);
    // The windows' data the info panels and a sector's opening read, asked for with the move.
    loadData(S.holder(r.stage)).then(function () {
      if (!active || token !== goToken) return;
      if (animation) finishAnimation();
      const kept = view;
      flyCamera(viewNow(), to, null, function () {
        if (token !== goToken) return;
        landCamera(kept, to);
        afterFlight();
      });
    }).finally(function () {
      if (token === goToken) arriving = false;
    });
  }

  // The view lands at camera `to`: the fit it flew from (`kept`, the display's) stays, the zoom is what keeps the distance.
  function landCamera(kept, to) {
    view = Object.assign({}, kept, { target: to.target.slice(), dist: to.dist, quat: to.quat.clone() });
    viewAtCamera(resolved, { target: to.target, dist: to.dist, quat: to.quat.toArray() });
    applyView();
  }

  // A camera {target, dist, quat} (the quaternion as an array or as a THREE one) with a THREE quaternion.
  function asCamera(camera) {
    return { target: camera.target.slice(), dist: camera.dist, quat: camera.quat.toArray ? camera.quat.clone() : new THREE.Quaternion().fromArray(camera.quat) };
  }

  // The camera has landed at a stage's frame: the sector selected or picked on the way is described, a sector asked to open opens
  // in place, and a region's densest landing goes on.
  function afterFlight() {
    hover = null;
    if (selectedSector) {
      showSectorInfo(sectorCell(selectedSector));
    } else if (pickedSector && !pinned && !pendingOpen) {
      describePickedSector();
    } else {
      showStageInfo();
    }
    if (pendingOpen) {
      const asked = pendingOpen;
      pendingOpen = null;
      const block = sectorCell(asked);
      const record = sectorRecord(block);
      if (record) enterSector(block, record, { push: false });
      else if (host.fetchUnchartedScene) enterSector(block, null, { push: false });
    }
    if (landingAfter) landDensest();
  }

  // MAP.124: the whole galaxy opens as close as its zoom range allows
  // while every charted sector shows, centered on them; zooming out still
  // reaches the whole galaxy. Only on the first view, with no sector asked.
  function closeOnCharted(r, v) {
    if (!isWholeGalaxy(r) || pendingOpen || selectedSector) return;
    const data = dataFor(S.holder(r.stage));
    const charted = r.view.blocks.filter(function (block) { return generatedOf(block, data) > 0; });
    if (!charted.length) return;
    const fit = centeredFit(S.boundsPoints(charted), v.target.slice(), v.quat);
    const full = fitDistance(fitPointsOf(r), v.fit.target, v.quat);
    const zoom = Math.max(GALAXY_MIN_ZOOM, Math.min(1, fit.dist / full));
    if (zoom >= 1) return;
    v.zoom = zoom;
    v.target = fit.target.slice();
    v.dist = zoom * full;
  }

  function flyTo(r) {
    closeAmbient();
    const from = display && view ? viewNow() : null;
    const to = cameraFor(r);
    // MAP.228: the camera the address asked for, taken up once the stage is in place.
    const want = pendingCamera;
    pendingCamera = null;
    const old = display;
    const incoming = buildDisplay(r);
    setDisplayFade(incoming, from ? 0 : 1);
    display = incoming;
    if (!from) {
      view = settledView(to);
      if (want) viewAtCamera(r, want);
      else closeOnCharted(r, view);
      applyView();
      afterArrival();
      return;
    }
    if (old) leaving.push(old);
    flyCamera(from, want ? { target: want.target.slice(), dist: want.dist, quat: new THREE.Quaternion().fromArray(want.quat) } : to, function (t) {
      const ms = animation ? animation.duration : 0;
      const fadeStart = ms > 0 ? Math.max(0, 1 - FADE_IN_MS / ms) : 0;
      const f = fadeStart >= 1 ? 1 : Math.max(0, (t - fadeStart) / (1 - fadeStart));
      setDisplayFade(incoming, f);
      if (old) setDisplayFade(old, 1 - f);
    }, function () {
      setDisplayFade(display, 1);
      view = settledView(to);
      if (want) viewAtCamera(r, want);
      applyView();
      afterArrival();
    });
  }

  function afterArrival() {
    leaving.forEach(function (d) { if (d && d !== display) disposeDisplay(d); });
    leaving = [];
    hover = null;
    const index = selectedSector ? sectorOption(selectedSector) : -1;
    if (index >= 0) {
      setHover({ option: index, sticky: true });
      showSectorInfo(display.options[index].blocks[0]);
    } else if (pickedSector && !selectedSector && !pinned && !pendingOpen) {
      // Picked from a coarser view (MAP.265): not a choice of this stage, but outlined and described all the same.
      applyHover();
      describePickedSector();
    } else {
      applyHover();
      showStageInfo();
    }
    if (pendingOpen) {
      const asked = pendingOpen;
      pendingOpen = null;
      const at = sectorOption(asked);
      const block = at >= 0 ? display.options[at].blocks[0] : null;
      const record = block ? sectorRecord(block) : null;
      if (record) enterSector(block, record, { push: false });
      else if (block && host.fetchUnchartedScene) enterSector(block, null, { push: false });
    }
    if (landingAfter) landDensest();
  }

  function finishAnimation() {
    const a = animation;
    animation = null;
    a.frame(1);
    a.done();
  }

  // Every frame: the running move, if any.
  function step(now) {
    if (!active) return;
    if (camChangedAt && now - camChangedAt >= CAM_URL_DELAY_MS) writeCameraAddress();
    sectorStage.update();
    systemStage.update(Date.now());
    updateAmbient(now);
    updateSystemAmbient(now);
    if (followBody && !animation && view) {
      const at = systemStage.where(followBody);
      if (at) {
        view.target = at.center;
        applyView();
      }
    }
    if (zoomGoal) chaseZoomGoal();
    if (!animation) return;
    const a = animation;
    const t = a.duration > 0 ? Math.min(1, (now - a.startedAt) / a.duration) : 1;
    a.frame(t);
    if (t >= 1) {
      animation = null;
      a.done();
    }
  }

  // The map's selection (NAV.13, picker.js): the stage it shows, or the
  // sector picked on it. Keys, the Up and Reset buttons and the breadcrumb
  // move it, and the map follows (below); renderCrumbs records where the
  // map is.
  const selection = createSelection({
    parentOf: function (ref) {
      if (ref.kind === "sector") return { kind: "stage", stage: ref.stage };
      const parent = S.parentStage(ref.stage, getOutline(), edgePc);
      return parent ? { kind: "stage", stage: parent } : null;
    },
    same: function (a, b) {
      if (a.kind !== b.kind) return false;
      if (a.kind === "stage") return S.sameStage(a.stage, b.stage);
      return a.shell === b.shell && a.band === b.band && a.wedge === b.wedge;
    },
  });
  selection.onChange(function (change) {
    if (change.ref.kind === "stage") flyToStage(change.ref.stage);
  });

  // Whether the selected sector is the one the window is centred on (the
  // window's own crumb names it).
  function isSelectedCentre(of) {
    const win = (of || stage).win;
    return !!(selectedSector && win && win.centre.m === 1 && win.centre.shell === selectedSector.shell
      && win.centre.band === selectedSector.band && win.centre.wedge === selectedSector.wedge);
  }

  function syncSelection() {
    if (selectedSector && !isSelectedCentre()) {
      selection.select({ kind: "sector", shell: selectedSector.shell, band: selectedSector.band, wedge: selectedSector.wedge, stage: stage }, { quiet: true });
    } else {
      selection.select({ kind: "stage", stage: stage }, { quiet: true });
    }
  }

  function up(options) {
    if (pinned) return;
    const keep = !!(options && options.keepCamera && view);
    if (enteredSystem) {
      // Out of the system, to the sector it sits in.
      closeSystem();
      pushEntry(withKept(sectorOpenQuery()));
      const back = cameraForSector(entered);
      if (keep) keepFit(back, entered.fitPoints);
      else {
        flyCamera(viewNow(), back, null, function () {
          view = settledView(back);
          applyView();
        });
      }
      renderCrumbs();
      return;
    }
    if (entered) {
      // Out of the sector, to the stage it sits in with it still selected.
      closeSector();
      pushEntry(withKept(sectorQuery(selectedSector, false)));
      const back = cameraFor(resolve(stage));
      if (!keep) {
        const kept = view;
        flyCamera(viewNow(), back, null, function () { landCamera(kept, back); });
      }
      applyHover();
      renderCrumbs();
      renderTravel();
      return;
    }
    selection.up();
  }

  // MAP.263: the camera stays as it is and `fit` (a stage's own camera, round `points`) becomes its fit: the zoom is
  // what keeps the distance.
  function keepFit(fit, points) {
    const full = fitDistance(points, fit.target, view.quat);
    view.fit = fit;
    if (full > 0) view.zoom = view.dist / full;
    applyView();
  }

  function home() {
    if (pinned) return;
    const top = S.galaxyStage();
    // Already there: go again, so the camera returns to the whole galaxy.
    if (S.sameStage(stage, top) && !selectedSector) flyToStage(top);
    else selection.select({ kind: "stage", stage: top });
  }

  // The map's Back and Forward: the browser's own history, kept to this
  // map's stages.
  function travel(direction) {
    if (direction < 0 && mapIndex > 0) history.back();
    else if (direction > 0 && mapIndex < maxIndex) history.forward();
  }

  // Forward to current (MAP.95): all the way to the newest view in this
  // map's history, in one jump.
  function travelToCurrent() {
    if (mapIndex < maxIndex) history.go(maxIndex - mapIndex);
  }

  // MAP.121: the window one cell over, at the level the window is at, or null where there is nothing to step to (the
  // whole galaxy, past a pole or the core, or where nothing of the galaxy lies). With a sector opened in place the step
  // is to the sector beside it, `{ sector }` (the helper that finds it is the walls' icons', MAP.327); with a system
  // opened, to the nearest system of the sector that lies that way, `{ system }` (the System Map has no grid).
  function stepTarget(direction) {
    if (pinned) return null;
    if (enteredSystem) {
      // MAP.121: from the open system, the nearest system of the open sector that lies that way (systemstep.js).
      // MAP.329: the systems across the sector's edge are candidates too (crossStars), once their sectors' scenes are read.
      const here = sectorStage.entries().find(function (e) { return e.endpoint === "system:" + enteredSystem.id; });
      const to = here ? nearestInDirection(here, sectorStars(true).concat(crossStars()), direction, function (e) { return e === here; }) : null;
      return to ? { system: to } : null;
    }
    if (entered) {
      const here = entered.sector;
      const to = P.neighbourBlock({ m: 1, shell: here.shell, band: here.band, wedge: here.wedge }, direction, edgePc);
      return to && S.sectorStage(to.shell, to.band, to.wedge) ? { sector: { shell: to.shell, band: to.band, wedge: to.wedge } } : null;
    }
    if (!stage.win) return null;
    const to = P.neighbourBlock(stage.win.centre, direction, edgePc);
    if (!to) return null;
    const next = { win: { centre: S.stripBlock(to), g: stage.win.g } };
    return resolve(next).problem ? null : next;
  }

  // Steps the window to the neighbouring cell in `direction` (north, south, east, west, out, in): one camera move to
  // its frame, the address carrying the camera and the breadcrumb following the new centre. With a sector opened in
  // place it opens the sector beside it the same way (MAP.327's icon on the wall).
  function stepToNeighbour(direction) {
    if (animation) return false;
    const next = stepTarget(direction);
    if (!next) return false;
    if (next.sector) return openSector(next.sector);
    if (next.system) {
      const mode = systemStage.mode();
      closeSystem();
      if (next.system.sector) {
        // MAP.329: a system across the sector's edge: its sector opens beside this one (as the wall icon's step does) and its system with it,
        // each its own history entry.
        pendingSystem = next.system.endpoint.slice("system:".length);
        pendingSystemOptions = { push: true, mode: mode };
        return openSector(next.system.sector);
      }
      enterSystem(next.system, { push: true, mode: mode });
      return true;
    }
    flyToStage(next);
    return true;
  }

  // MAP.315: flies the camera to a block (the region an infobox names), as a click on its crumb does: its window,
  // at the width the window is now (or the default). False while the camera is moving or when the block holds
  // nothing to show.
  function openBlock(block) {
    if (!view || animation || !P.isBlock(block)) return false;
    const next = S.windowStage(block, stage.win ? stage.win.g : P.DEFAULT_WINDOW);
    if (resolve(next).problem) return false;
    flyToStage(next);
    return true;
  }

  // MAP.319: flies to a region: its 9-block's window, with the camera where the region's bounding sphere
  // (`target` = {center: [x, y, z] pc, radius pc}) fills the view, keeping the viewing direction. Leaves a sector
  // opened in place. False while the camera is moving or when the block holds nothing to show.
  function openRegion(block, target) {
    if (!view || animation || !P.isBlock(block) || !target || !target.center) return false;
    const next = S.windowStage(block, P.DEFAULT_WINDOW);
    if (resolve(next).problem) return false;
    const want = {
      target: target.center.slice(), quat: view.quat.toArray(),
      dist: FIT.fitDistance(target.radius, camera.fov, aspectNow(), { minRadius: FIT.MIN_VIEW_RADIUS_PC }),
    };
    flyToStage(next, { camera: want });
    return true;
  }

  function renderTravel() {
    if (!els.controls) return;
    els.controls.querySelectorAll("[data-step]").forEach(function (button) {
      button.disabled = !stepTarget(button.dataset.step);
    });
    const back = els.controls.querySelector('[data-action="back"]');
    const forward = els.controls.querySelector('[data-action="forward"]');
    const upButton = els.controls.querySelector('[data-action="up"]');
    if (back) back.disabled = mapIndex <= 0;
    if (forward) forward.disabled = mapIndex >= maxIndex;
    newestButton.disabled = mapIndex >= maxIndex;
    if (upButton) upButton.disabled = !stage.win && !selectedSector;
    const resetButton = els.controls.querySelector('[data-action="reset-view"]');
    if (resetButton) resetButton.disabled = !isFree(resolved);
  }

  // --- Picking and hover -----------------------------------------------------

  // The stage's choices are the picker's last layer (mappick.js): a
  // phenomenon the map puts in front takes the pointer first.
  const picker = host.picker || createPicker(camera, canvasEl);
  const sectorStage = createSectorStage({
    THREE: THREE, scene: host.scene, camera: camera, canvasEl: canvasEl, picker: picker,
    accentColor: host.accentColor || "#4f5fe8", lightBackground: !!host.lightBackground,
    pixelRatio: host.pixelRatio || function () { return 1; },
    fetchScene: host.fetchSectorScene, fetchUnchartedScene: host.fetchUnchartedScene, showInfo: host.showInfo || function () {}, navPick: host.navPick || null,
    hideCell: host.hideCell || function () {}, preview: host.preview || null,
    viewport: host.viewport || null, opened: host.opened, closed: host.closed, kindsChanged: kindsChanged,
    deselected: function () { host.showHint(pinned && host.restingHint ? host.restingHint : hintFor(resolved)); },
    closeness: function () { return view && view.zoom > 0 ? 1 / view.zoom : 1; },
    openSystem: !pinned && host.fetchSystemScene ? openSystemFromAmbient : null,
    flyToEntry: !pinned ? flyToEntry : null,
  });
  const systemStage = createSystemStage({
    THREE: THREE, scene: host.scene, camera: camera, canvasEl: canvasEl, picker: picker,
    accentColor: host.accentColor || "#4f5fe8",
    fetchScene: host.fetchSystemScene, showInfo: host.showInfo || function () {}, viewport: host.viewport || null,
    navPick: host.navPick || null,
    flyTo: function (entry) { flyToBody(entry); },
    setMode: function (mode) { setSystemMode(mode); },
    deselected: function () { host.showHint(hintFor(resolved)); },
    picked: function (ref) { writeObjectToUrl(ref); },
  });
  const tooltip = createTooltip(els.tooltip);
  const choiceLayer = picker.addLayer({
    name: "choices", priority: 10,
    enabled: function () { return !!display && !!resolved; },
    meshes: function () {
      const meshes = [];
      display.groups.forEach(function (group) { Array.prototype.push.apply(meshes, group.meshes); });
      return meshes;
    },
    // Only what is drawn can be picked (MAP.149): not the back faces (the
    // walls of the block the camera is inside) and not what the near field
    // has faded out, so the pick falls through to what is behind it.
    visible: function (hit, ray) {
      if (hit.face && ray && hit.face.normal.dot(ray.direction) > 0) return false;
      return NF.visibleEnough(host.nearFieldAt(hit.point.x, hit.point.y, hit.point.z));
    },
    entryOf: function (hit) {
      const cell = hit.object.userData.cells[hit.object.userData.part.owners[hit.face.a]];
      return cell ? columnChoice(cell.option, hit.point) : null;
    },
  });

  // MAP.297: the blocks are not picked by the pointer on the Galaxy Map (the map is about locations, not sector boxes):
  // a hover, click or double click goes to the stars, objects and markers, or to nothing. They are still the layer the
  // wheel's crossing from one window to the next reads (blockUnderWheel). A sector page's map (pinned) keeps its
  // neighbouring sectors pickable: a click there opens that sector's page. Choosing a NAV end (needGenerated) is
  // the one thing on the Galaxy Map that is still a choice among blocks.
  function blocksPickable() {
    return !!pinned || needGenerated;
  }
  function pointerLayers(layer) {
    return blocksPickable() || layer !== choiceLayer;
  }
  // What the wheel's anchoring and crossings look at: not the map's own stars (MAP.304), which take the pointer's
  // click and hover only, so a zoom aims where it did.
  function notMapStars(layer) {
    return layer.name !== "stars";
  }

  // Whether choice `index` can be picked: not while choosing a NAV end
  // and it holds nothing generated.
  function pickable(index) {
    if (!display || !display.options[index]) return false;
    if (!needGenerated) return true;
    const data = display.data;
    return display.options[index].blocks.some(function (b) { return generatedOf(b, data) > 0; });
  }

  // The whole galaxy is seen from above, where a ray meets the topmost of
  // a column's stacked blocks: the choice is the block of the column
  // (the blocks above and below the point `point` of the plane) holding the
  // most generated sectors, else the one nearest the midplane.
  function columnChoice(index, point) {
    if (!display || !isWholeGalaxy(display.resolved) || !display.options[index] || !point) return index;
    const off = function (b) { return Math.abs((b.c0 + b.c1) / 2 - Math.PI / 2); };
    let best = index;
    let bestGenerated = generatedOf(display.options[index].blocks[0], display.data);
    let bestOff = off(display.options[index].blocks[0].bounds);
    display.options.forEach(function (option, n) {
      const candidate = option.blocks[0];
      if (n === index || !S.columnContains(candidate.bounds, point.x, point.y)) return;
      const generated = generatedOf(candidate, display.data);
      const away = off(candidate.bounds);
      if (generated > bestGenerated || (generated === bestGenerated && away < bestOff)) {
        best = n;
        bestGenerated = generated;
        bestOff = away;
      }
    });
    return best;
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

  function addOutline(points, opacity) {
    if (!outlineGroup) {
      outlineGroup = new THREE.Group();
      outlineGroup.renderOrder = 6;
      host.scene.add(outlineGroup);
    }
    const line = lineObject(
      points,
      new THREE.LineBasicMaterial({ color: accent, depthTest: false, transparent: true, opacity: opacity }),
      "block outline",
    );
    line.renderOrder = 6;
    line.frustumCulled = false;
    outlineGroup.add(line);
  }

  // An accent outline round the whole block a click would open (MAP.18,
  // MAP.52, MAP.198): every edge of its cell, and a faint shade over its
  // faces, on top of the map.
  function outlineOption(index, keep) {
    if (!keep) clearOutline();
    const option = display.options[index];
    const points = [];
    option.blocks.forEach(function (block) { Array.prototype.push.apply(points, S.cellEdgePoints(block.bounds)); });
    addOutline(points, 1);
    // A sector opened in place keeps its outline and no shade: its own scene
    // is drawn inside, and picking its stars must not be dimmed.
    if (entered) return;
    const shadeBounds = option.blocks.map(function (block) { return block.bounds; });
    const shadeOrigin = centreOf([].concat.apply([], shadeBounds.map(function (b) { return P.cellCoordinates(b).cartesian; })));
    const built = P.buildPrismGeometry(shadeBounds, shadeOrigin);
    const geometry = new THREE.BufferGeometry();
    geometry.setAttribute("position", new THREE.BufferAttribute(built.positions, 3));
    geometry.setIndex(new THREE.BufferAttribute(built.indices, 1));
    const fill = new THREE.Mesh(
      geometry,
      new THREE.MeshBasicMaterial({ color: accent, transparent: true, opacity: HOVER_SHADE_OPACITY, depthTest: false, depthWrite: false, side: THREE.DoubleSide }),
    );
    fill.position.set(shadeOrigin[0], shadeOrigin[1], shadeOrigin[2]);
    fill.userData.relativeLayer = "block shade";
    fill.renderOrder = 5;
    fill.frustumCulled = false;
    fill.userData.hoverShade = true;
    outlineGroup.add(fill);
  }

  function setHover(next) {
    hover = blocksPickable() ? next : null;
    applyHover();
  }

  // The hovered choice stays as it is and the others fade, and gets an
  // outline (a stage of very many choices is one mesh: only the outline).
  function applyHover() {
    if (!display) return;
    const index = hover && hover.option != null ? hover.option : -1;
    const ambientIndex = ambientShown() ? sectorOption(ambient.block) : -1;
    display.groups.forEach(function (group, n) {
      if (!group.option) {
        group.fade = group.baseFade = 1;
        return;
      }
      group.fade = index < 0 || n === index ? 1 : OTHER_FADE;
      // The sector opened in place is drawn by its own scene: its block
      // is only a faint cell round it.
      if (entered) group.fade = sectorOption(entered.sector) === n ? ENTERED_BLOCK_FADE : Math.min(group.fade, ENTERED_OTHERS_FADE);
      group.baseFade = group.fade;
      // ... and so is the one the camera is coming up to, as much as its scene shows (MAP.152).
      if (n === ambientIndex) group.fade = Math.min(group.fade, ambientBlockFade());
      // MAP.202: the lines round a block are nearly invisible at rest and
      // show in full on the block being chosen.
      const edges = index === n ? 1 : group.restEdges;
      group.meshes.forEach(function (mesh) {
        if (mesh.material.uniforms && mesh.material.uniforms.gridEdges && group.restEdges > 0) mesh.material.uniforms.gridEdges.value = edges;
      });
    });
    if (!animation) setDisplayFade(display, 1);
    // MAP.297: no sector box is outlined or shaded on the Galaxy Map, hovered, selected or entered.
    if (!blocksPickable()) {
      clearOutline();
      return;
    }
    // MAP.201: the sector selected (one that is not opened, as an
    // unfilled one cannot be) stays outlined until something else is
    // selected, whether or not the pointer is on it.
    const selected = selectedSector && !entered ? sectorOption(selectedSector) : -1;
    if (index >= 0 && display.options[index]) outlineOption(index);
    else clearOutline();
    if (selected >= 0 && selected !== index) outlineOption(selected, true);
    else if (selected < 0 && !selectedSector && pickedSector && !entered && !pinned) addOutline(S.cellEdgePoints(sectorCell(pickedSector).bounds), 1);
  }

  function showTooltip(text, clientX, clientY) {
    tooltip.show(text, clientX, clientY);
  }

  // The tooltip beside a choice picked by keyboard: at its middle on
  // screen.
  function tooltipAtOption(index) {
    const option = display.options[index];
    const middle = S.extentOf(option.blocks).center;
    tooltip.showAt(optionText(index), camera, canvasEl, middle[0], middle[1], middle[2]);
  }

  function blockText(block, data) {
    const b = block.bounds;
    const generated = generatedOf(block, data);
    const longitude = S.bearingRange(b);
    const latitude = S.latitudeRange(b);
    const parts = [(block.m === 1 ? "Sector " : "Block ") + S.blockLabel(block),
      P.isFullTurn(b) ? "all longitudes" : "longitude " + longitude[0].toFixed(1) + "°–" + longitude[1].toFixed(1) + "°",
      "latitude " + latitude[1].toFixed(1) + "°–" + latitude[0].toFixed(1) + "°",
      Math.round(b.r0) + "–" + Math.round(b.r1) + " pc from the core"];
    if (block.m === 1) parts.push(generated > 0 ? "charted" : "uncharted");
    else parts.push(S.formatCount(generated) + (block.total != null ? " of " + S.formatCount(block.total) : "") + " sectors charted");
    return parts.join(", ");
  }

  function optionText(index) {
    return blockText(display.options[index].blocks[0], display.data) + (pickable(index) ? "" : " (nothing charted here to pick)");
  }

  // The block the panel on arrival describes: the cell the window is
  // centred on, or for a window of sectors the level-3 block holding its
  // centre (the sector's own panel opens when it is selected).
  function infoBlock() {
    const centre = stage.win.centre;
    return centre.m === 1 ? S.drillParent(centre) : centre;
  }

  // How many sectors `block` (the infoBlock) holds that are generated, from
  // the window's data.
  function generatedIn(block, data) {
    if (!data) return null;
    if (stage.win.centre.m !== 1) return data.generated.get(cellKey(block)) || 0;
    let generated = 0;
    data.sectors.forEach(function (sector) {
      const parent = S.drillParent({ m: 1, shell: sector.shell, band: sector.band, wedge: sector.wedge });
      if (S.sameBlock(parent, block)) generated += 1;
    });
    return generated;
  }

  // The panel on arrival: the hint alone for the whole galaxy, else the
  // window's middle block (its hint under it).
  function showStageInfo() {
    const hint = hintFor(resolved);
    if (!stage.win) {
      host.showHint(hint);
      return;
    }
    const block = infoBlock();
    const data = dataFor(S.holder(stage));
    const generated = generatedIn(block, data);
    const total = getOutline().shapeless ? null : S.drillBlockTotal(block, getOutline());
    host.showBlockInfo({ block: block, total: total, generated: generated, hint: hint, enter: null });
  }

  function sectorEntry(block, sector) {
    const middle = P.cellCenter(block.bounds);
    return {
      id: sector.id, name: sector.name, system_count: sector.system_count,
      shell_index: block.shell, band_index: block.band, wedge_index: block.wedge,
      designation: S.sectorDesignation(block.shell, block.band, block.wedge),
      galactic_radius_pc: Math.hypot(middle[0], middle[1], middle[2]),
      position: middle,
    };
  }

  function sectorRecord(block) {
    const data = dataFor(S.holder(stage));
    return data ? data.sectors.get(cellKey(block)) : null;
  }

  function showSectorInfo(block) {
    // A sector page names its sector itself.
    if (pinned) return;
    const sector = sectorRecord(block) || (isEnteredCell(block) ? nearRecords.get(cellKey(block)) : null);
    if (sector) host.showPlacedInfo(sectorEntry(block, sector));
    else host.showCellInfo({
      m: 1, bounds: block.bounds, address: { shell: block.shell, band: block.band, wedge: block.wedge }, filled: 0,
      uncharted: isEnteredCell(block) && entered.uncharted ? sectorStage.unchartedSummary() : null,
    });
  }

  // Whether `block` is the sector open in place.
  function isEnteredCell(block) {
    return !!entered && entered.sector.shell === block.shell && entered.sector.band === block.band
      && entered.sector.wedge === block.wedge;
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
    if (!blocksPickable()) {
      return "Click a star or an object for its details, double-click to fly to it, or scroll to zoom toward it.";
    }
    if (isSectorView(r)) {
      return "Click a sector to open it; a sector that isn't charted yet shows where it is"
        + (host.canGenerate ? " and how to chart it." : ".");
    }
    return "Click a block to zoom in on it; Up zooms out again.";
  }

  // --- Acting on a choice ----------------------------------------------------

  function act(index, at, how) {
    if (!pickable(index)) return;
    const option = display.options[index];
    if (pinned) {
      // A sector page keeps to its sector: another sector's page opens
      // instead, and nothing else is a step to take.
      const other = option.blocks.length === 1 && option.blocks[0].m === 1 ? sectorRecord(option.blocks[0]) : null;
      if (other && host.openSectorPage && (!entered || entered.id !== other.id)) host.openSectorPage(other.id);
      return;
    }
    // A sector is opened.
    if (option.blocks[0].m === 1) {
      const block = option.blocks[0];
      const record = sectorRecord(block);
      if (record && host.fetchSectorScene) {
        // Opened in place, unless it is the one open already (MAP.66).
        if (!entered || entered.id !== record.id) enterSector(block, record, { push: true, fly: !(how && how.keepCamera) });
        return;
      }
      if (!record && host.fetchUnchartedScene && !unchartedNone.has(cellKey(block))) {
        // Nothing generated here: what the scatters left in it still opens
        // (MAP.162); an empty cell is only selected, as below.
        if (!isEnteredCell(block)) enterSector(block, null, { push: true, fly: !(how && how.keepCamera) });
        return;
      }
      if (entered) closeSector();
      selectedSector = { shell: block.shell, band: block.band, wedge: block.wedge };
      setHover({ option: index, sticky: true });
      showSectorInfo(block);
      renderCrumbs();
      renderTravel();
      return;
    }
    takePick(option.pick, at, how);
  }

  // Takes `pick` (a cell) of the choices on screen: the window one size down,
  // centred on it. A choice is made on the stage that is DRAWN, not the one
  // `go` has already moved to while the next window's data loads: the cell
  // clicked is a cell of the drawn stage, whichever stage that is (MAP.107).
  function takePick(pick, at, how) {
    if (pinned) return;
    const base = display && display.resolved ? display.resolved.stage : stage;
    // The window opens round the middle of what the click went through, not the face it met first.
    if (at) at = S.chordMiddle(pick.block.bounds || P.drillBlockBounds(pick.block, edgePc), at, host.camera.position);
    const into = S.windowInto(base, pick.block, at, edgePc);
    if (!into) return;
    // MAP.265: any sector is selectable from any distance: the sector holding the middle of what the click went
    // through is found in closed form (no block of sectors need be drawn), and is the one selected on arrival.
    const sector = at ? P.cellAt(1, at.x, at.y, at.z, edgePc) : null;
    const chosen = sector && P.cellExists(sector.shell, sector.band, sector.wedge)
      ? { shell: sector.shell, band: sector.band, wedge: sector.wedge } : null;
    // (A crossing the wheel makes only carries what is picked on.)
    const clicked = !(how && how.keepCamera);
    if (clicked) pickedSector = chosen;
    flyToStage(into, Object.assign({ keepPicked: clicked && !!chosen }, how));
  }

  // A sector cell {shell, band, wedge} as the info panels and outlines take it.
  function sectorCell(sector) {
    const cell = { m: 1, shell: sector.shell, band: sector.band, wedge: sector.wedge };
    cell.bounds = P.drillBlockBounds(cell, edgePc);
    return cell;
  }

  // The panel of the selected sector when it is no choice on screen: its own record, asked for as the smallest window
  // round it (the stage's counts are of coarser cells), and the cell's place when nothing is charted there.
  function describePickedSector() {
    const asked = pickedSector;
    const cell = sectorCell(asked);
    loadData({ centre: cell, g: P.WINDOW_WIDTHS[0] }).then(function (data) {
      if (!active || pickedSector !== asked || selectedSector || entered || pinned) return;
      const record = data.sectors.get(cellKey(cell));
      if (record) host.showPlacedInfo(sectorEntry(cell, record));
      else host.showCellInfo({ m: 1, bounds: cell.bounds, address: asked, filled: data.failed ? null : 0, uncharted: null });
    });
  }

  canvasEl.galaxySelectedSector = function () {
    const sector = selectedSector || pickedSector;
    return sector ? Object.assign({}, sector) : null;
  };

  function sectorOption(sector) {
    if (!display) return -1;
    return display.options.findIndex(function (option) {
      return option.blocks.length === 1 && option.blocks[0].m === 1 && option.blocks[0].shell === sector.shell
        && option.blocks[0].band === sector.band && option.blocks[0].wedge === sector.wedge;
    });
  }

  // --- The sector the camera comes up to (MAP.152) ----------------------------------

  // In a window of sectors a sector's own scene (every star and phenomenon of it) fades in as
  // the camera comes within sectorreach.js's FADE_IN_EDGES of its centre, with no stage change,
  // flight or history entry, and is dropped past UNLOAD_EDGES. It is not "the sector opened":
  // `entered` stays null, the panels and the URL are the window's, and a click on the sector
  // opens it as ever (taking the scene over as it is). `ambient` is {key, id, block, record,
  // loading} while its scene is up or being fetched.
  let ambient = null;
  // The sector in range that has not stayed there for SR.DWELL_MS yet: {key, since}.
  let ambientWant = null;
  let ambientToken = 0;
  let ambientShare = 0;
  let ambientEye = null;
  let ambientSeenIn = null;
  // A sector being opened in place (enterSector) until it has arrived.
  let entering = false;

  function ambientShown() {
    return !!ambient && !ambient.loading && ambientShare > 0;
  }

  // The share of its block left drawn once its scene shows: ENTERED_BLOCK_FADE at the most.
  function ambientBlockFade() {
    return 1 - ambientShare * (1 - ENTERED_BLOCK_FADE);
  }

  // The camera's position (parsecs).
  function eyeNow() {
    const back = new THREE.Vector3(0, 0, view.dist).applyQuaternion(view.quat);
    return [view.target[0] + back.x, view.target[1] + back.y, view.target[2] + back.z];
  }

  // The sectors of the window on screen that have a scene to show: [{key, center, block, record}].
  function ambientCandidates() {
    if (!display.ambientCandidates) {
      const failed = display.ambientFailed || (display.ambientFailed = new Set());
      display.ambientCandidates = [];
      display.options.forEach(function (option) {
        const block = option.blocks[0];
        const record = block && block.m === 1 ? sectorRecord(block) : null;
        if (!record || failed.has(cellKey(block))) return;
        display.ambientCandidates.push({ key: cellKey(block), center: P.cellCenter(block.bounds), block: block, record: record });
      });
    }
    return display.ambientCandidates;
  }

  // MAP.335: in a window of bigger cells there are no sector records to hold scenes by, so the sectors the eye is near are asked of the
  // server (`host.fetchNear`, GET /galaxy/near): the generated ones within NEAR_ASK_EDGES sector edges of the eye, read again when the eye has
  // gone NEAR_REUSE_EDGES from where it was asked. The scene is then chosen from the camera alone, with no window of sectors in front of it.
  const NEAR_ASK_EDGES = 4;
  const NEAR_REUSE_EDGES = 1;
  // Past this view distance, in sector edges, no sector's scene can be in range of the eye and the server is not asked.
  const NEAR_MAX_VIEW_EDGES = 60;
  const NEAR_RETRY_MS = 5000;
  let nearRead = null;
  let nearAsk = null;
  let nearToken = 0;
  let nearRetryAt = 0;
  // Every sector the server has named, by cell: what the sector's panel and scene need (its id and name).
  const nearRecords = new Map();
  const nearFailed = new Set();

  function distanceBetween(a, b) {
    return Math.hypot(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
  }

  function askNear(eye) {
    const token = ++nearToken;
    const at = eye.slice();
    nearAsk = { at: at };
    const query = "?x=" + at[0].toFixed(3) + "&y=" + at[1].toFixed(3) + "&z=" + at[2].toFixed(3) + "&r=" + (NEAR_ASK_EDGES * edgePc).toFixed(3);
    host.fetchNear(query).then(function (payload) {
      if (!active || token !== nearToken) return;
      nearAsk = null;
      const candidates = [];
      ((payload && payload.sectors) || []).forEach(function (sector) {
        const block = { m: 1, shell: sector.shell, band: sector.band, wedge: sector.wedge };
        block.bounds = P.drillBlockBounds(block, edgePc);
        const key = cellKey(block);
        nearRecords.set(key, sector);
        if (!nearFailed.has(key)) candidates.push({ key: key, center: P.cellCenter(block.bounds), block: block, record: sector });
      });
      nearRead = { at: at, candidates: candidates };
      ambientEye = null;
    }, function () {
      if (!active || token !== nearToken) return;
      nearAsk = null;
      nearRetryAt = performance.now() + NEAR_RETRY_MS;
      ambientEye = null;
    });
  }

  // The sectors of the near read for an eye at `eye`: [{key, center, block, record}], asking again when the last read was made from further
  // than NEAR_REUSE_EDGES away (the last read stands meanwhile). Null when the map has no such read.
  function nearCandidates(eye, now) {
    if (!host.fetchNear) return null;
    const reuse = NEAR_REUSE_EDGES * edgePc;
    if ((!nearRead || distanceBetween(eye, nearRead.at) > reuse) && (!nearAsk || distanceBetween(eye, nearAsk.at) > reuse) && now >= nearRetryAt) askNear(eye);
    return nearRead ? nearRead.candidates : [];
  }

  // Drops the sector scene the camera came up to, and what was waiting for its dwell.
  function closeAmbient() {
    const was = ambient;
    ambientToken += 1;
    ambientWant = null;
    ambientEye = null;
    ambient = null;
    ambientShare = 0;
    // MAP.335: a star's system comes up out of the sector's scene, so it goes with it.
    if (was && !entered && systemAmbient && !enteredSystem) closeSystemAmbient();
    if (was && !entering && sectorStage.isAmbient()) sectorStage.close();
    if (was && display) applyHover();
  }

  // The state forgotten without touching the scene (it is the sector opened in place now, or
  // already closed).
  function forgetAmbient() {
    ambientToken += 1;
    ambientWant = null;
    ambientEye = null;
    ambient = null;
    ambientShare = 0;
  }

  function setAmbientShare(share) {
    if (share === ambientShare || (share > 0 && share < 1 && Math.abs(share - ambientShare) < 0.004)) return;
    ambientShare = share;
    sectorStage.setFade(share);
    // Only the block of the sector coming in changes: not a whole redraw of the hover.
    if (!display || animation) return;
    const n = sectorOption(ambient.block);
    const group = n >= 0 ? display.groups[n] : null;
    if (!group) return;
    group.fade = Math.min(group.baseFade == null ? 1 : group.baseFade, ambientBlockFade());
    setDisplayFade(display, 1);
  }

  function loadAmbient(candidate) {
    const token = ++ambientToken;
    ambientWant = null;
    ambientShare = 0;
    ambient = { key: candidate.key, id: candidate.record.id, block: candidate.block, record: candidate.record, loading: true };
    const mine = ambient;
    // The block of the scene that was up (if any) is drawn whole again while this one loads.
    if (display) applyHover();
    const failed = function () {
      if (!active || token !== ambientToken || ambient !== mine) return;
      // Left out of this window's candidates, so it is not asked for again and again.
      if (display && display.ambientFailed) display.ambientFailed.add(candidate.key);
      if (display) display.ambientCandidates = null;
      nearFailed.add(candidate.key);
      if (nearRead) nearRead.candidates = nearRead.candidates.filter(function (c) { return c.key !== candidate.key; });
      ambient = null;
    };
    sectorStage.open(candidate.record.id, candidate.block.bounds, null, true).then(function (geo) {
      if (!active || token !== ambientToken || ambient !== mine) return;
      if (geo === undefined) {
        // TEST.450: another open took the scene stage's place meanwhile (a cell the wheel tried): not a failure of this sector, which is
        // asked for again.
        ambient = null;
        return;
      }
      if (!geo) {
        failed();
        return;
      }
      mine.loading = false;
      ambientEye = null;
    }, failed);
  }

  // Every frame: which sector scene the camera wants up and how much of it shows.
  function updateAmbient(now) {
    // MAP.335: in a window of sector cells the scenes are held by its records; in any other stage (a window of bigger cells, the whole galaxy,
    // a camera that has flown in) they come from the sectors near the eye, once the view is near enough for one to be in range.
    const sectorWindow = !!(view && display && resolved && isSectorView(resolved));
    const fromEye = !sectorWindow && !!(host.fetchNear && view && display && resolved && view.dist <= NEAR_MAX_VIEW_EDGES * edgePc);
    const eligible = !pinned && (sectorWindow || fromEye) && !entered && !enteredSystem;
    if (!eligible) {
      if (ambient && !entered) closeAmbient();
      return;
    }
    if (animation || entering) return;
    const eye = eyeNow();
    if (ambientEye && !ambientWant && ambientSeenIn === display
      && eye[0] === ambientEye[0] && eye[1] === ambientEye[1] && eye[2] === ambientEye[2]) return;
    ambientEye = eye;
    ambientSeenIn = display;
    const candidates = sectorWindow ? ambientCandidates() : nearCandidates(eye, now);
    const chosen = SR.chooseSector(eye, candidates, edgePc, ambient ? ambient.key : null);
    if (!chosen) {
      if (ambient || ambientWant) closeAmbient();
      return;
    }
    if (ambient && ambient.key === chosen.key) {
      ambientWant = null;
      if (!ambient.loading) setAmbientShare(chosen.share);
      return;
    }
    if (!ambientWant || ambientWant.key !== chosen.key) {
      ambientWant = { key: chosen.key, since: now };
      return;
    }
    if (now - ambientWant.since < SR.DWELL_MS) return;
    const candidate = candidates.find(function (c) { return c.key === chosen.key; });
    if (candidate) loadAmbient(candidate);
  }

  // A star of the scene the camera came up to: its sector opens where the camera is (no
  // flight), then its system.
  function openSystemFromAmbient(entry) {
    if (entered) {
      enterSystem(entry, { push: true });
      return;
    }
    if (!ambient || ambient.loading) return;
    enterSector(ambient.block, ambient.record, {
      push: false, fly: false, then: function () { enterSystem(entry, { push: true }); },
    });
  }

  // For the browser tests: the sector scene up for the camera coming near it, and how much of it shows.
  canvasEl.galaxyAmbient = function () {
    return ambient && !ambient.loading ? { id: ambient.id, key: ambient.key, share: ambientShare } : null;
  };

  // For the browser tests (MAP.329): the name of the system a step in `direction` goes to from the open system, or null.
  canvasEl.galaxyStepTo = function (direction) {
    const next = stepTarget(direction);
    return next && next.system ? next.system.name : null;
  };

  // For the browser tests (TEST.442): the sector opened in place once its scene is up and takes clicks, {id, share}, or null
  // (the panel names the sector from its address before the scene has been fetched, and a click then falls on a cell).
  canvasEl.galaxyOpenSector = function () {
    return sectorStage.isOpen() ? { id: sectorStage.sectorId(), share: sectorStage.fade() } : null;
  };

  // For the browser tests (MAP.360): the share of the map's width and of its height the opened sector's fit points cover, {x, y}, or null.
  canvasEl.galaxySectorFitShare = function () {
    if (!entered) return null;
    camera.updateMatrixWorld();
    const lo = [Infinity, Infinity];
    const hi = [-Infinity, -Infinity];
    const point = new THREE.Vector3();
    for (let k = 0; k + 2 < entered.fitPoints.length; k += 3) {
      point.set(entered.fitPoints[k], entered.fitPoints[k + 1], entered.fitPoints[k + 2]).project(camera);
      lo[0] = Math.min(lo[0], point.x); hi[0] = Math.max(hi[0], point.x);
      lo[1] = Math.min(lo[1], point.y); hi[1] = Math.max(hi[1], point.y);
    }
    return { x: (hi[0] - lo[0]) / 2, y: (hi[1] - lo[1]) / 2 };
  };

  // For the browser tests (MAP.302): every entry of the sector opened in place, {key, kind, name, x, y, z, hidden}, in the map's frame
  // (parsecs) as galaxyScreenOf takes it, so a test can look for one of every kind the sector holds and click where it is drawn.
  canvasEl.galaxySectorEntries = function () {
    if (!sectorStage.isOpen()) return [];
    return sectorStage.entries().map(function (entry) {
      return {
        key: entry.key || null, kind: sectorStage.kindOf(entry), name: entry.name || "", x: entry.x, y: entry.y, z: entry.z,
        hidden: sectorStage.hides(entry),
      };
    });
  };

  // --- The sector opened in place (MAP.66) ------------------------------------

  // Closes the open sector's scene; the stage it sits in stays.
  function closeSector() {
    enterToken += 1;
    closeSystem();
    sectorStage.close();
    entered = null;
    entering = false;
    forgetAmbient();
  }

  // Opens a generated sector in place: its scene joins the map, the
  // camera flies to fit it, and (push) it takes its own history entry,
  // ?sector=<designation>&open=1. A sector that can't be opened (not
  // placed, or its scene didn't load) is only selected, as an ungenerated
  // one is. `block` is its cell, `record` the stage's sector record; with no
  // record the cell is one nothing was generated in, which opens (as
  // "Uncharted") only when the scatters left objects in it (MAP.162).
  // TEST.450: the cells nothing was generated in that the server has said hold nothing to open (a 404, or a scene with no objects), so the
  // wheel is not left asking for the same cell on every tick.
  const unchartedNone = new Set();

  function enterSector(block, record, options) {
    const token = ++enterToken;
    const at = { shell: block.shell, band: block.band, wedge: block.wedge };
    entering = true;
    notice("Opening " + (record ? record.name || "the sector" : "the sector") + "…");
    sectorStage.open(record ? record.id : null, block.bounds, record ? null : at).then(function (geo) {
      if (!active || token !== enterToken) return;
      entering = false;
      notice("");
      selectedSector = at;
      if (!geo) {
        if (!record && geo === null) unchartedNone.add(cellKey(block));
        if (!record && entered) closeSector();
        const index = sectorOption(at);
        if (index >= 0) setHover({ option: index, sticky: true });
        showSectorInfo(block);
        renderCrumbs();
        renderTravel();
        return;
      }
      forgetAmbient();
      entered = {
        sector: at, id: record ? record.id : null, uncharted: !record, center: geo.center, halfEdge: geo.halfEdge,
        fitPoints: geo.fitPoints,
      };
      if (options.push) pushEntry(withKept(sectorQuery(at, true)));
      showSectorInfo(block);
      applyHover();
      renderCrumbs();
      renderTravel();
      if (animation) finishAnimation();
      if (options.fly === false) {
        // Opened where the camera is (the sector the camera came up to): no flight. Its fit is the
        // view's now, the zoom what keeps the camera where it is.
        const fit = cameraForSector(entered);
        view.fit = fit;
        view.zoom = view.dist / fitDistance(entered.fitPoints, fit.target, view.quat);
        applyView();
        if (options.then) options.then();
        return;
      }
      const to = cameraForSector(entered);
      flyCamera(viewNow(), to, null, function () {
        view = settledView(to);
        applyView();
        openPendingSystem();
      });
    }, function () {
      if (!active || token !== enterToken) return;
      entering = false;
      if (!record) {
        // Only a bonus on a cell nothing was generated in: it is just selected.
        unchartedNone.add(cellKey(block));
        notice("");
        if (entered) closeSector();
        selectedSector = at;
        const index = sectorOption(at);
        if (index >= 0) setHover({ option: index, sticky: true });
        showSectorInfo(block);
        renderCrumbs();
        renderTravel();
        return;
      }
      notice("The sector could not be loaded. Please try again shortly.");
    });
  }

  // --- The system the camera comes up to (MAP.152) --------------------------------------

  // In an opened sector a star's system (its planets, moons and orbits) fades in as its disc comes to about
  // LOAD_PX pixels across (systemreach.js), with no flight or history entry, and is dropped below UNLOAD_PX. It is
  // not "the system opened": `enteredSystem` stays null and the scene cannot be picked; one more wheel tick at
  // the nearest zoom opens it from the scene already there (the camera then flies to fit it as ever).
  // `systemAmbient` is {key (the system's href), entry, loading} while its scene is up or being fetched.
  let systemAmbient = null;
  let systemWant = null;
  let systemToken = 0;
  let systemShare = 0;
  let systemEye = null;
  const systemFailed = new Set();
  // The real radius (pc) of each system seen so far, by href: a system is taken for 50 AU across until its scene
  // has been loaded once. Without this a small system would be loaded at 20 px of 50 AU, found too small to stay
  // (under 10 px), dropped, and loaded again, over and over.
  const systemRadii = new Map();

  function systemCandidates() {
    // MAP.335: of the sector opened, or of the sector scene the camera came up to.
    const holder = entered || ambient;
    if (!holder.systemCandidates) {
      holder.systemCandidates = sectorStage.sceneEntries().filter(function (e) {
        return e.href && e.endpoint && e.endpoint.indexOf("system:") === 0 && !systemFailed.has(e.href);
      }).map(function (e) {
        return { key: e.href, center: [e.x, e.y, e.z], radiusPc: null, entry: e };
      });
    }
    holder.systemCandidates.forEach(function (c) {
      c.radiusPc = systemRadii.get(c.key) || null;
    });
    return holder.systemCandidates;
  }

  // Drops the system scene the camera came up to.
  function closeSystemAmbient() {
    const was = systemAmbient;
    forgetSystemAmbient();
    if (was && systemStage.isAmbient()) systemStage.close();
  }

  function forgetSystemAmbient() {
    systemToken += 1;
    systemWant = null;
    systemEye = null;
    systemAmbient = null;
    systemShare = 0;
  }

  function setSystemShare(share) {
    if (share === systemShare || (share > 0 && share < 1 && Math.abs(share - systemShare) < 0.004)) return;
    systemShare = share;
    systemStage.setFade(share);
  }

  function loadSystemAmbient(candidate) {
    const token = ++systemToken;
    systemWant = null;
    systemShare = 0;
    const mine = systemAmbient = { key: candidate.key, entry: candidate.entry, loading: true };
    const failed = function () {
      if (!active || token !== systemToken || systemAmbient !== mine) return;
      systemFailed.add(candidate.key);
      if (entered) entered.systemCandidates = null;
      if (ambient) ambient.systemCandidates = null;
      systemAmbient = null;
    };
    systemStage.open(candidate.entry.href, candidate.entry.name, [candidate.entry.x, candidate.entry.y, candidate.entry.z], undefined, true).then(function (geo) {
      if (!active || token !== systemToken || systemAmbient !== mine) return;
      if (geo === undefined) {
        // TEST.450: another open took the scene stage's place meanwhile (a cell the wheel tried): not a failure of this sector, which is
        // asked for again.
        ambient = null;
        return;
      }
      if (!geo) {
        failed();
        return;
      }
      mine.loading = false;
      systemRadii.set(candidate.key, geo.radiusPc);
      systemEye = null;
    }, failed);
  }

  // Every frame: which system scene the camera wants up and how much of it shows.
  function updateSystemAmbient(now) {
    const eligible = !pinned && view && (entered || ambientShown()) && !enteredSystem && host.fetchSystemScene;
    if (!eligible) {
      if (systemAmbient && !enteredSystem) closeSystemAmbient();
      return;
    }
    if (animation) return;
    const eye = eyeNow();
    if (systemEye && !systemWant && eye[0] === systemEye[0] && eye[1] === systemEye[1] && eye[2] === systemEye[2]) return;
    systemEye = eye;
    const chosen = SYS.chooseSystem(eye, systemCandidates(), { heightPx: canvasEl.clientHeight || 600, fovDeg: camera.fov }, systemAmbient ? systemAmbient.key : null);
    if (!chosen) {
      if (systemAmbient || systemWant) closeSystemAmbient();
      return;
    }
    if (systemAmbient && systemAmbient.key === chosen.key) {
      systemWant = null;
      if (!systemAmbient.loading) setSystemShare(chosen.share);
      return;
    }
    if (!systemWant || systemWant.key !== chosen.key) {
      systemWant = { key: chosen.key, since: now };
      return;
    }
    if (now - systemWant.since < SYS.DWELL_MS) return;
    const candidate = systemCandidates().find(function (c) { return c.key === chosen.key; });
    if (candidate) loadSystemAmbient(candidate);
  }

  // For the browser tests: the system scene up for the camera coming near it, and how much of it shows.
  canvasEl.galaxySystemAmbient = function () {
    return systemAmbient && !systemAmbient.loading ? { href: systemAmbient.key, share: systemShare } : null;
  };

  // --- The system opened in place (MAP.125) ----------------------------------

  // The URL of the open sector, ?sector=<designation>&open=1.
  function sectorOpenQuery() {
    return sectorQuery(selectedSector, true);
  }

  // The URL of a sector in the stage on screen: ?sector=<designation>, with &open=1 when it is open. The camera
  // (?cam=, MAP.228) is added by the address writer once the view is at rest.
  function sectorQuery(sector, open) {
    return "?sector=" + S.sectorDesignation(sector.shell, sector.band, sector.wedge) + (open ? "&open=1" : "");
  }

  function closeSystem() {
    enterSystemToken += 1;
    followBody = null;
    systemStage.close();
    enteredSystem = null;
    forgetSystemAmbient();
  }

  // Opens `entry`'s system (a star of the open sector) in place: its scene
  // joins the map at the star, the camera flies to fit it, and (push) it takes
  // its own history entry, ?sector=...&open=1&system=<id>.
  function enterSystem(entry, options) {
    if (!entered || pinned || !host.fetchSystemScene || !entry || !entry.href) return;
    const token = ++enterSystemToken;
    notice("Opening " + (entry.name || "the system") + "…");
    systemStage.open(entry.href, entry.name, [entry.x, entry.y, entry.z], options.mode).then(function (geo) {
      if (!active || token !== enterSystemToken) return;
      notice("");
      if (!geo) return;
      forgetSystemAmbient();
      enteredSystem = {
        id: systemStage.systemId(), href: entry.href, name: entry.name, center: geo.center,
        radiusPc: geo.radiusPc, fitPoints: geo.fitPoints,
      };
      if (options.push) pushEntry(withKept(sectorOpenQuery() + "&system=" + enteredSystem.id));
      renderCrumbs();
      renderTravel();
      if (animation) finishAnimation();
      const to = cameraForSystem(enteredSystem);
      const object = pendingObject;
      pendingObject = null;
      if (options.fly === false) {
        // MAP.263: opened where the camera is, by the wheel: no flight, the camera's turn and distance kept.
        keepFit(to, enteredSystem.fitPoints);
        return;
      }
      flyCamera(viewNow(), to, null, function () {
        view = settledView(to);
        applyView();
        const entry = object ? systemStage.entry(object) : null;
        if (entry) {
          systemStage.select(object);
          flyToBody(entry);
        }
      });
    }, function () {
      if (!active || token !== enterSystemToken) return;
      notice("The system could not be loaded. Please try again shortly.");
    });
  }

  // The system a URL named, once its sector is open.
  function openPendingSystem() {
    if (!pendingSystem || !entered) return;
    const wanted = pendingSystem;
    pendingSystem = null;
    const entry = sectorStage.entries().find(function (e) { return e.endpoint === "system:" + wanted; });
    const options = pendingSystemOptions || { push: false };
    pendingSystemOptions = null;
    if (entry) enterSystem(entry, options);
  }

  // The map's scale of the open system, true or compressed.
  function setSystemMode(mode) {
    const geo = systemStage.setMode(mode);
    if (!geo || !enteredSystem) return;
    enteredSystem.radiusPc = geo.radiusPc;
    enteredSystem.fitPoints = geo.fitPoints;
    const to = cameraForSystem(enteredSystem);
    flyCamera(viewNow(), to, null, function () {
      view = settledView(to);
      applyView();
    });
  }

  // Flies the camera to a body of the open system, close enough to see it
  // with its moons about it.
  function flyToBody(entry) {
    const at = systemStage.where(entry.ref);
    if (!at || !view) return;
    const from = viewNow();
    const fit = view.fit;
    const to = { target: at.center.slice(), dist: Math.max(at.radius * 6, 1e-12), quat: view.quat.clone() };
    flyCamera(from, to, null, function () {
      const fitDist = fitDistance(enteredSystem.fitPoints, fit.target, to.quat);
      view = Object.assign({}, to, { target: to.target.slice(), quat: to.quat.clone(), fit: fit, zoom: to.dist / fitDist });
      followBody = entry.ref;
      applyView();
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
      followBody = null;
      const perPx = worldUnitsPerPixel(camera, view.dist, canvasEl.clientHeight);
      const before = view.target.slice();
      MC.panInScreenPlane(THREE, camera, view.target, dx, dy, perPx);
      // A sector page's map keeps its target inside the sector: a move that would leave it is not taken.
      if (pinned && !NF.blockContains(sectorCell(pinned).bounds, view.target[0], view.target[1], view.target[2])) {
        view.target[0] = before[0]; view.target[1] = before[1]; view.target[2] = before[2];
      }
    }
    applyView();
  }

  // MAP.245: the stars of the open sector the wheel and the clicks can take for a thing: those with a system, not
  // left off the map.
  function sectorStars(inSystem) {
    if (!(entered || ambientShown()) || (enteredSystem && !inSystem)) return [];
    return sectorStage.sceneEntries().filter(function (e) {
      return e.href && e.endpoint && e.endpoint.indexOf("system:") === 0 && !sectorStage.sceneHides(e);
    });
  }

  // MAP.329: the systems of the sectors round the open one, for stepping from an open system across the sector's edge: the generated sectors
  // near the system are asked of the server (`host.fetchNear`, as the ambient scenes are) and each one's scene read for its systems
  // (crossstars.js), nothing drawn. Read once for each system opened; the step buttons come on when it is in (renderTravel).
  let cross = null;
  function crossStars() {
    if (!enteredSystem || !entered || !host.fetchNear || !host.fetchSectorScene) return [];
    if (!cross || cross.id !== enteredSystem.id) loadCross();
    return cross.stars;
  }

  function loadCross() {
    const mine = cross = { id: enteredSystem.id, stars: [] };
    const here = sectorStage.entries().find(function (e) { return e.endpoint === "system:" + enteredSystem.id; });
    if (!here) return;
    const centreOf = function (record) {
      const block = { m: 1, shell: record.shell, band: record.band, wedge: record.wedge };
      return P.cellCenter(P.drillBlockBounds(block, edgePc));
    };
    const query = "?x=" + here.x.toFixed(3) + "&y=" + here.y.toFixed(3) + "&z=" + here.z.toFixed(3) + "&r=" + (2 * edgePc).toFixed(3);
    const openId = entered.id;
    host.fetchNear(query).then(function (payload) {
      if (!active || cross !== mine) return null;
      const records = CROSS.neighbourRecords((payload && payload.sectors) || [], centreOf, [here.x, here.y, here.z], edgePc, openId);
      return Promise.all(records.map(function (record) {
        return host.fetchSectorScene(record.id).then(function (data) { return CROSS.starsOfScene(data, record); }, function () { return []; });
      }));
    }).then(function (lists) {
      if (!lists || !active || cross !== mine) return;
      mine.stars = [].concat.apply([], lists);
      renderTravel();
    }, function () { /* the neighbours are only a bonus: the step stays to the open sector */ });
  }

  // MAP.245: the star of the open sector under a screen point (within STAR_ANCHOR_PX of it on the screen), as
  // [x, y, z], or null. The stars are not picked (MAP.101), so this is what lets the wheel know one is there.
  // MAP.300: with no sector open it is a star the map draws from its nodes' lists (the brightest of those level with the point), so the
  // wheel comes in on any star of a dense cluster that is drawn, not on the focus plane or the column behind it.
  function starUnder(clientX, clientY) {
    const sector = sectorStars();
    const rect = canvasEl.getBoundingClientRect();
    camera.updateMatrixWorld();
    const reach = { reach: function () { return STAR_ANCHOR_PX; } };
    const found = sector.length ? nearestOnScreen(sector, camera, rect, clientX, clientY, reach) : null;
    if (found) return [found.entry.x, found.entry.y, found.entry.z];
    const drawn = host.drawnStars ? host.drawnStars() : [];
    if (!drawn.length) return null;
    const drawnFound = nearestOnScreen(drawn, camera, rect, clientX, clientY, { reach: reach.reach, lastWins: true });
    return drawnFound ? [drawnFound.entry.x, drawnFound.entry.y, drawnFound.entry.z] : null;
  }

  // MAP.150: the world point under a screen point: a star of the open sector there (MAP.245), else the nearest
  // thing the picker finds there, else the point on the plane through the view's middle (the focus) at right
  // angles to the way the camera looks. In the open sector the picker's cells and walls lie behind the stars and
  // are never taken for the point: the way to the point must not end on a mesh behind what is pointed at.
  function anchorUnder(clientX, clientY) {
    const thing = nearestThingAlong(clientX, clientY);
    return thing ? thing.point : null;
  }

  // MAP.336: the nearest named thing of the scene the camera is in, for the ruler's readout: a body of the open system, else a star of the sector
  // scene (up, ambient or open), by its distance from the eye; null with none.
  function nearestNamed() {
    const eye = camera.position;
    let best = null;
    const consider = function (entry, kind, radius) {
      if (!entry || !entry.name || entry.x == null) return;
      const distance = Math.max(0, Math.hypot(entry.x - eye.x, entry.y - eye.y, entry.z - eye.z) - (radius || 0));
      if (!best || distance < best.distance) best = { name: entry.name, kind: kind, distance: distance };
    };
    if (enteredSystem) systemStage.entries().forEach(function (e) { consider(e, "body", e.r); });
    if (!best) sectorStage.sceneEntries().forEach(function (e) { if (!sectorStage.sceneHides(e)) consider(e, "star", 0); });
    return best;
  }

  // MAP.334: the nearest thing along the ray through a screen point, which the zoom steps are shares of: {point, kind, distance} with the point
  // in the world (parsecs), the distance `h` from the eye to it and the kind "star", "body" (a picked body of the open system, with its
  // `radius`), "pick" (any other thing the picker found), "crowd" (the densest layer of drawn stars in front of the thing the wheel would otherwise head for, MAP.301), "column" (the galaxy's stars thickening along the ray) or "plane" (nothing drawn
  // there: the focus plane's point). The same law in every stage: a wheel step covers a share of h and never takes the eye within the clearance
  // of the thing (`limitedStep`), so there is no wall and no landing on a surface.
  function nearestThingAlong(clientX, clientY) {
    const found = thingAlong(clientX, clientY);
    if (!found || found.kind === "body") return found;
    // MAP.301: a crowd of drawn stars in front of that thing is what the wheel heads for first.
    const rect = canvasEl.getBoundingClientRect();
    const ray = new THREE.Raycaster();
    ray.setFromCamera(new THREE.Vector2(((clientX - rect.left) / rect.width) * 2 - 1, -((clientY - rect.top) / rect.height) * 2 + 1), camera);
    const crowd = crowdAlong(ray.ray.origin, ray.ray.direction, found.distance);
    if (!crowd) return found;
    const eye = camera.position;
    return { point: crowd, kind: "crowd", distance: Math.hypot(crowd[0] - eye.x, crowd[1] - eye.y, crowd[2] - eye.z) };
  }

  function thingAlong(clientX, clientY) {
    const rect = canvasEl.getBoundingClientRect();
    if (!view || !rect.width || !rect.height) return null;
    const eye = camera.position;
    const thing = function (point, kind, extra) {
      return Object.assign({ point: point, kind: kind, distance: Math.hypot(point[0] - eye.x, point[1] - eye.y, point[2] - eye.z) }, extra);
    };
    const star = starUnder(clientX, clientY);
    if (star) return thing(star, "star");
    const ray = new THREE.Raycaster();
    ray.setFromCamera(new THREE.Vector2(((clientX - rect.left) / rect.width) * 2 - 1, -((clientY - rect.top) / rect.height) * 2 + 1), camera);
    const origin = ray.ray.origin;
    const direction = ray.ray.direction;
    const found = picker.pick(clientX, clientY, notMapStars);
    const forward = new THREE.Vector3(view.target[0], view.target[1], view.target[2]).sub(origin).normalize();
    const onPlane = MC.pointOnFocusPlane(origin.toArray(), direction.toArray(), view.target, forward.toArray());
    // A pointer on the cross zooms about the cross (MAP.199), whatever the column behind it holds.
    const onCross = Math.hypot(clientX - (rect.left + rect.width / 2), clientY - (rect.top + rect.height / 2)) <= CROSS_PX;
    // MAP.292: the blocks are picking geometry unless they are drawn (MAP.265), and the first face of one nobody sees
    // (the top of a column of them lies kiloparsecs above the disc the eye is looking at) is no place to zoom toward:
    // the point is the focus plane's under the pointer, which a step in and a step out leave where it was.
    if (found && found.layer === choiceLayer && !entered && !cellsShown() && onPlane) return thing(onPlane, "plane");
    if (found && found.distance > 0 && !(entered && found.layer === choiceLayer)) {
      const body = enteredSystem && found.entry && found.entry.r > 0 ? found.entry : null;
      return thing([origin.x + direction.x * found.distance, origin.y + direction.y * found.distance, origin.z + direction.z * found.distance],
        body ? "body" : "pick", body ? { radius: body.r } : null);
    }
    const column = !onCross && columnAnchor(origin, direction, forward, onPlane);
    if (column) return thing(column, "column");
    return onPlane ? thing(onPlane, "plane") : null;
  }

  // MAP.301: the densest layer of the drawn stars along the cursor ray, nearer the eye than `farthest` (where the wheel would otherwise head),
  // or null: no stars supplied, a system open, or no crowd of CROWD_MIN_STARS in front of that.
  function crowdAlong(origin, direction, farthest) {
    if (!host.drawnStars || enteredSystem || !(farthest > 0)) return null;
    const found = MC.densestLayerAlong(origin.toArray(), direction.toArray(), host.drawnStars(),
      { cone: CROWD_CONE, nearest: 2 * ZOOM_CLEARANCE * view.dist, farthest: farthest, minPoints: CROWD_MIN_STARS });
    return found ? found.point : null;
  }

  // MAP.247: the first point of the cursor ray where the column of the galaxy's stars passes the threshold
  // (COLUMN_SHARE), or null: no shape yet, a system open, nothing dense enough within reach.
  function columnAnchor(origin, direction, forward, onPlane) {
    if (!host.shape || enteredSystem || !onPlane) return null;
    const focus = Math.hypot(onPlane[0] - origin.x, onPlane[1] - origin.y, onPlane[2] - origin.z);
    const found = MC.columnDensityPoint(origin.toArray(), direction.toArray(),
      function (x, y, z) { return P.relativeDensity(x, y, z, host.shape); }, COLUMN_SHARE * view.dist, COLUMN_REACH * focus);
    if (!found) return null;
    const nearest = COLUMN_NEAREST * focus;
    const away = Math.hypot(found[0] - origin.x, found[1] - origin.y, found[2] - origin.z);
    if (!(away > 0)) return null;
    if (away >= nearest) return found;
    const k = nearest / away;
    return [origin.x + (found[0] - origin.x) * k, origin.y + (found[1] - origin.y) * k, origin.z + (found[2] - origin.z) * k];
  }

  // Zooms by `factor` (under 1 comes closer) about the view's middle: the
  // one centre the view turns about too (MAP.199).
  // MAP.292: while the next window's data is on its way (`arriving`, the camera kept where it is) the wheel is not held
  // to the end of the range it carried the view to: the view goes on coming closer meanwhile, so a scroll is not stopped
  // for as long as the data takes, and the next window takes the camera over where it has got to.
  // MAP.334: the range at which the stage on show hands over to the next (or the one before): the ends of its own zoom, which the wheel
  // carries the view to and across (`zoomAcross`), not the ends of what the camera may do (`zoomLimits`).
  function handoverPolicy(at) {
    const policy = zoomPolicyFor(resolved, at);
    return arriving && MC.canZoom(policy) ? MC.zoomPolicy(policy.kind, policy.nearest * ARRIVING_OVERSHOOT, policy.farthest) : policy;
  }

  // MAP.334: where the wheel, plus and minus, the slider and a pinch may take the camera: one absolute range in every stage, from the whole
  // galaxy in view down to FIELD_FLOOR_PC outside an opened system and GLOBAL_FLOOR_PC (a moon framed) in one, as a ratio to the stage's fit like
  // a policy. No stage ends it (it ends where the stage hands over, `handoverPolicy`, and the next one takes the camera where it is); only
  // while the next stage's data is on its way is the view held to the end of the range it was carried to (MAP.292).
  function zoomLimits(at) {
    const handover = handoverPolicy(at);
    if (!MC.canZoom(handover) || arriving) return relaxed(handover);
    const fit = fitNow();
    if (!(fit > 0)) return relaxed(handover);
    // The nearest is the floor, except that outside a system it is never short of the stage's own hand-over (a floor above where a star's
    // system opens would be a wall); a system has nothing nearer to hand over to.
    const nearest = enteredSystem ? GLOBAL_FLOOR_PC / fit : Math.min(handover.nearest, FIELD_FLOOR_PC / fit);
    if (pinned) return MC.zoomPolicy(MC.ZOOM_RANGE, nearest, PINNED_MAX_ZOOM);
    return MC.zoomPolicy(MC.ZOOM_RANGE, nearest, Math.max(handover.farthest, globalFarPc() / fit));
  }

  function zoomBy(factor) {
    if (!view || animation || !MC.canZoom(zoomPolicyFor(resolved))) return;
    view.zoom = MC.clampDistance(zoomLimits(), view.zoom * factor, 1);
    reframe();
    applyView();
  }

  // MAP.333: plus and minus ask for a distance, and cross the stages as the wheel does (`zoomToDistance`).
  function zoomButton(factor) {
    if (!view || !MC.canZoom(zoomPolicyFor(resolved))) return;
    zoomToDistance(view.dist * factor);
  }

  // MAP.333: the policy with its range opened to where the view is: a view an explicit move left beyond its stage's range (a fly-to
  // that came closer than the stage's nearest) is not pushed back into it, only held from going further out of it.
  function relaxed(policy) {
    if (!view || !policy || policy.kind !== MC.ZOOM_RANGE) return policy;
    return MC.zoomPolicy(policy.kind, Math.min(policy.nearest, view.zoom), Math.max(policy.farthest, view.zoom));
  }

  // MAP.333: the farthest the camera goes, the whole galaxy in view: the whole-galaxy stage's own camera distance.
  let galaxyFarPc = 0;
  function globalFarPc() {
    if (pinned) return (fitNow() || cameraFor(resolved).dist) * PINNED_MAX_ZOOM;   // a sector page's map ends at its sector
    if (!galaxyFarPc) {
      const whole = S.resolveStage(S.galaxyStage(), getOutline(), edgePc);
      galaxyFarPc = whole && !whole.problem && whole.view && whole.view.blocks.length ? cameraFor(whole).dist : host.galaxyRadius * GALAXY_FAR_RADII;
    }
    return galaxyFarPc;
  }

  // MAP.333: the controls that ask for a camera distance (the slider, plus and minus, a pinch) take it as a goal, whole decades at a time if
  // need be: within the stage's range the zoom is set; beyond it the zoom goes to the end of the range and then crosses to the next stage as
  // a wheel tick at the end does, again and again (each frame, while the next stage's data is on its way) until the distance is reached.
  let zoomGoal = 0;
  let zoomGoalInward = true;
  let zoomStalledAt = 0;
  const ZOOM_GOAL_PATIENCE_MS = 3000;
  function zoomToDistance(dist) {
    if (!view || !(dist > 0)) return;
    finishSoftFlight();
    zoomStalledAt = 0;
    zoomGoal = Math.max(GLOBAL_FLOOR_PC, Math.min(globalFarPc(), dist));
    zoomGoalInward = zoomGoal < view.dist;
    chaseZoomGoal();
  }

  // TEST.458: a move the person asks for (Reset view, a fly-to, a centring) ends a distance goal still being chased: left alive, it
  // took the view back out to where the plus or minus presses had been taking it as soon as the move was over.
  function dropZoomGoal() {
    zoomGoal = 0;
    zoomStalledAt = 0;
  }

  // A camera-only flight (a fly to a point, a centring) is finished at once where it was going, so the controls take the view from there.
  function finishSoftFlight() {
    if (!animation || !animation.soft) return;
    const flight = animation;
    animation = null;
    flight.frame(1);
    flight.done();
  }

  function chaseZoomGoal() {
    if (!zoomGoal || !view || !view.fit || animation) return;
    // Reached; or passed, which a crossing can do (the next stage's own view sits at another distance: a system's, once a sector's star is
    // gone in to): the goal is one way, so the chase never turns back on it.
    if (Math.abs(Math.log(zoomGoal / view.dist)) < 1e-3 || (zoomGoalInward ? view.dist < zoomGoal : view.dist > zoomGoal)) {
      zoomGoal = 0;
      zoomStalledAt = 0;
      return;
    }
    if (arriving) return;
    const policy = handoverPolicy();
    const fit = fitNow();
    if (!(fit > 0)) {
      zoomGoal = 0;
      return;
    }
    const want = zoomGoal / fit;
    const inward = zoomGoalInward;
    if (!MC.canZoom(policy) || (want >= policy.nearest && want <= policy.farthest)) {
      view.zoom = want;
      reframe();
      applyView();
      zoomGoal = 0;
      return;
    }
    if (!entered && !enteredSystem) {
      // Out of an opened sector or system the stage has no range's end to stop at: the zoom is the camera's alone.
      view.zoom = MC.clampDistance(zoomLimits(), want, 1);
      reframe();
      applyView();
      if (Math.abs(Math.log(zoomGoal / view.dist)) < 1e-3) zoomGoal = 0;
      return;
    }
    const end = inward ? policy.nearest : policy.farthest;
    if (inward ? view.zoom > end * 1.001 : view.zoom < end * 0.999) {
      view.zoom = end;
      reframe();
      applyView();
      return;
    }
    // At the end of the range and no crossing yet: what the next stage shows may still be on its way (a sector's scene waits for a
    // moment's dwell), so the goal is kept for a while, as a wheel kept turning would keep it.
    if (zoomAcross(inward ? -1 : 1, null)) {
      zoomStalledAt = 0;
    } else {
      const now = performance.now();
      if (!zoomStalledAt) zoomStalledAt = now;
      else if (now - zoomStalledAt > ZOOM_GOAL_PATIENCE_MS) zoomGoal = 0;
    }
  }

  // MAP.225, MAP.245: the zoom factor a step toward `anchor` may take: no more than half the way to it (`MC.limitApproach`),
  // and never through a star of the open sector (`MC.limitPathClearance`; the clearance is a share of the view's
  // distance, ZOOM_CLEARANCE).
  function limitedStep(anchor, factor) {
    const eye = camera.position;
    const away = Math.hypot(anchor[0] - eye.x, anchor[1] - eye.y, anchor[2] - eye.z);
    const clearance = ZOOM_CLEARANCE * view.dist;
    let asked = MC.limitApproach(factor, away, clearance);
    if (asked < 1) {
      const stars = sectorStars();
      if (stars.length) {
        const to = [anchor[0] + asked * (eye.x - anchor[0]), anchor[1] + asked * (eye.y - anchor[1]), anchor[2] + asked * (eye.z - anchor[2])];
        const allowed = MC.limitPathClearance([eye.x, eye.y, eye.z], to,
          stars.map(function (e) { return [e.x, e.y, e.z]; }), clearance);
        asked = 1 - allowed * (1 - asked);
      }
    }
    return asked;
  }

  // MAP.225: the wheel zooms toward the point under the pointer, which stays where it is on the screen: the
  // picked star, remnant or cell, else the point of the focus plane under it (`anchorUnder`). A step in moves the
  // eye by no more than half the way to that point (`MC.limitApproach`), so it is a logarithmic approach and the
  // camera never lands on a star by accident. A view that cannot be moved (a locked one, a stage fitted whole), or
  // a wheel with no pointer over the map, zooms about the middle as before.
  function zoomTowardPointer(factor, event) {
    if (!view || animation || !MC.canZoom(zoomPolicyFor(resolved))) return;
    const anchor = event && isFree(resolved) && view.fit ? anchorUnder(event.clientX, event.clientY) : null;
    if (!anchor) {
      zoomBy(factor);
      return;
    }
    const asked = limitedStep(anchor, factor);
    const before = view.dist;
    // MAP.326: how near the camera may come is judged from the point being zoomed toward (a star under the pointer),
    // not from the view's middle, which can be a sector's width from that star: judged from the middle the wheel
    // stopped short of the star and never reached its system.
    view.zoom = MC.clampDistance(zoomLimits(anchor), view.zoom * asked, 1);
    reframe();
    const taken = before > 0 ? view.dist / before : 1;
    view.target = MC.targetAboutAnchor(view.target, anchor, taken);
    applyView();
  }

  // MAP.200: the zoom as a share of the policy's log range, 0 at the
  // farthest the view goes and 1 at the nearest.
  // MAP.333: the slider is the camera's distance on one log scale, 0 at the whole galaxy and 1 at the floor, the same wherever the view is.
  function zoomShare() {
    if (!view || !(view.dist > 0)) return 0;
    const far = Math.log(globalFarPc());
    const share = (far - Math.log(view.dist)) / (far - Math.log(GLOBAL_FLOOR_PC));
    return Math.max(0, Math.min(1, share));
  }

  function setZoomShare(share) {
    if (!view || !MC.canZoom(zoomPolicyFor(resolved))) return;
    const s = Math.max(0, Math.min(1, share));
    const far = Math.log(globalFarPc());
    zoomToDistance(Math.exp(far - s * (far - Math.log(GLOBAL_FLOOR_PC))));
  }

  // Moves the center `dx` across and `dy` up by shares of the map's width.
  function moveCentre(dx, dy) {
    if (!view || animation || !resolved || !isFree(resolved)) return;
    const px = canvasEl.clientWidth || 1;
    drag(-dx * px, dy * px, true);
  }

  // The center back to the middle of the region viewed, keeping the angle
  // and the zoom.
  function resetCentre() {
    if (!view || animation || !resolved || !view.fit || !isFree(resolved)) return;
    centreOn(view.fit.target.slice());
  }

  // The wheel zooms within the view's zoom policy (a locked view leaves
  // it to scroll the page). True when it was used.
  function onWheel(event) {
    if (!MC.canZoom(zoomPolicyFor(resolved)) || !view) return false;
    const deltaPx = MC.wheelPixels(event, canvasEl.clientHeight, 200);
    if (deltaPx && !zoomAcross(deltaPx, event)) zoomTowardPointer(Math.exp(deltaPx * WHEEL_ZOOM_PER_PX), event);
    return true;
  }

  // MAP.125, MAP.152: the wheel at the end of an opened sector's or system's own range carries the zoom on, with no stop and no flight: right
  // up to a sector whose scene is whole it opens where the camera is, zoomed right in on an opened sector the star picked (else the one
  // nearest the middle of the view, its system opens), and zoomed right out of a system or a sector the sector or the window comes back.
  // The windows are not crossed here (MAP.337): the camera's zoom is global and the display window stays as it is. True when this tick moved across.
  function zoomAcross(deltaPx, event) {
    if (animation || pinned || arriving) return false;
    // MAP.326: in an opened sector the range and the star to go in to are those of the point under the pointer.
    const under = event && entered && !enteredSystem && isFree(resolved) && view.fit ? anchorUnder(event.clientX, event.clientY) : null;
    const policy = zoomPolicyFor(resolved, under);
    if (deltaPx > 0 && view.zoom >= policy.farthest * 0.999 && (enteredSystem || entered)) {
      up({ keepCamera: true });
      return true;
    }
    if (deltaPx >= 0 || view.zoom > policy.nearest * 1.001) return false;
    if (entered && !enteredSystem) {
      if (!host.fetchSystemScene) return false;
      const picked = sectorStage.selectedEntry();
      const star = picked && picked.endpoint && picked.endpoint.indexOf("system:") === 0 && picked.href
        ? picked : (under && nearestStar(under, view.dist * NEAREST_STAR_REACH)) || nearestStar(view.target, view.dist * NEAREST_STAR_REACH);
      if (!star) return false;
      enterSystem(star, { push: true, fly: false });
      return true;
    }
    if (enteredSystem || !display || !resolved) return false;
    if (!entered && ambientShown() && ambientShare >= 0.999) {
      // Right up to a sector whose scene is whole: it opens where the camera is, and the next
      // tick can go on in to a star (MAP.152).
      enterSector(ambient.block, ambient.record, { push: true, fly: false });
      return true;
    }
    return false;
  }

  // The star of the open sector nearest `point` (parsecs): {entry, distance}, or null.
  function nearestStarTo(point) {
    let best = null;
    sectorStage.entries().forEach(function (e) {
      if (!e.href || !e.endpoint || e.endpoint.indexOf("system:") !== 0) return;
      const d = Math.hypot(e.x - point[0], e.y - point[1], e.z - point[2]);
      if (!best || d < best.distance) best = { entry: e, distance: d };
    });
    return best;
  }

  // The star of the open sector nearest `point` (parsecs) within `reach`.
  function nearestStar(point, reach) {
    let best = null;
    let bestDistance = reach;
    sectorStage.entries().forEach(function (e) {
      if (!e.href || !e.endpoint || e.endpoint.indexOf("system:") !== 0) return;
      const d = Math.hypot(e.x - point[0], e.y - point[1], e.z - point[2]);
      if (d < bestDistance) {
        best = e;
        bestDistance = d;
      }
    });
    return best;
  }

  // Back to the stage's own view, after turning or moving it.
  function resetView() {
    if (!resolved || animation || !display) return;
    dropZoomGoal();
    const from = viewNow();
    const own = enteredSystem || entered;
    const to = enteredSystem ? cameraForSystem(enteredSystem) : entered ? cameraForSector(entered) : cameraFor(resolve(stage));
    const kept = view;
    flyCamera(from, to, null, function () {
      if (own) {
        view = settledView(to);
        applyView();
      } else landCamera(kept, to);
    });
  }

  // Moves the orbit center to the selection (a star, cloud or body), or to
  // the middle of what the stage shows when nothing is selected, keeping the
  // angle and the zoom, so the mouse and the keys then turn and zoom about it
  // (MAP.138).
  function recenter() {
    if (!view || animation || !resolved || !view.fit) return;
    let goal = null;
    const selected = systemStage.selected();
    const where = selected ? systemStage.where(selected) : null;
    const entry = sectorStage.selectedEntry();
    if (where) goal = where.center.slice();
    else if (entry) goal = [entry.x, entry.y, entry.z];
    else goal = view.fit.target.slice();
    centreOn(goal);
  }

  // MAP.146, MAP.319: how far from `entry` the camera comes to look at it: a star with a system or a cloud or body of
  // a radius, to where its bounding sphere fills the view (flytofit.js; the zoom the view allows clamps it, so a
  // system comes no nearer than its disc whole); anything else keeps the distance the view has.
  function aspectNow() {
    return (canvasEl.clientWidth || 800) / (canvasEl.clientHeight || 600);
  }

  function distanceToLookAt(entry) {
    if (entry.href && entry.endpoint && entry.endpoint.indexOf("system:") === 0) {
      return FIT.fitDistance(systemRadii.get(entry.href) || SYS.NOMINAL_RADIUS_PC, camera.fov, aspectNow());
    }
    if (entry.r > 0) return FIT.fitDistance(entry.r, camera.fov, aspectNow(), { minRadius: FIT.MIN_VIEW_RADIUS_PC });
    // MAP.304: a star or point of the map comes to a sector and a half away, no farther than the view is now.
    if (entry.luminosity_sol != null) return Math.min(view.dist, 1.5 * edgePc);
    return view.dist;
  }

  // Flies the camera to look at `goal` (a world point) from `wanted` away, as near as the view's zoom allows, keeping
  // the angle: the view then turns and zooms about it (MAP.199).
  function flyToPoint(goal, wanted) {
    if (!view || animation || !resolved || !view.fit || !isFree(resolved)) return;
    dropZoomGoal();
    const point = goal.slice();
    const around = enteredSystem ? enteredSystem.fitPoints : entered ? entered.fitPoints : fitPointsOf(resolved);
    const fit = fitDistance(around, view.fit.target, view.quat);
    if (!(fit > 0)) return;
    // MAP.333: the distance asked for is the one the flight lands at; only the global limits hold it, not the range of the stage the
    // view is in (the wheel carries on across the stage's end from there).
    const zoom = Math.max(GLOBAL_FLOOR_PC, Math.min(globalFarPc(), wanted)) / fit;
    const to = { target: point, dist: zoom * fit, quat: view.quat.clone() };
    const kept = view;
    followBody = null;
    flyCamera(viewNow(), to, null, function () {
      view = Object.assign({}, kept, { target: point.slice(), dist: to.dist, quat: to.quat.clone(), zoom: zoom });
      applyView();
    });
    animation.soft = true;
  }

  // MAP.319: flies the camera to look at a bounding sphere, `target` = {center: [x, y, z] pc, radius pc}, from where it
  // fills the view (the region, system, moon system or nebula it bounds), keeping the viewing direction.
  function flyToFit(target, arrive) {
    if (!target || !target.center) return;
    // MAP.333: inside an opened system a sphere as small as a moon system is framed as closely as it asks; elsewhere the map's minimum view
    // radius holds a framing back from a point-sized sphere.
    flyToPoint(target.center, (arrive > 0 ? arrive : 1) * FIT.fitDistance(target.radius, camera.fov, aspectNow(), { minRadius: enteredSystem ? GLOBAL_FLOOR_PC : FIT.MIN_VIEW_RADIUS_PC }));
  }

  // Flies the camera to an object of the sector scene (a star, a phenomenon): a double click on it, its Go to
  // button, Enter with it selected.
  function flyToEntry(entry) {
    if (!view || !entry || entry.x == null) return;
    flyToPoint([entry.x, entry.y, entry.z], distanceToLookAt(entry));
  }

  // MAP.336: the one way to go to anything, the verb "Go to" on every button, name, link and result. `target` says where, by one of:
  // - {center: [x, y, z] pc, radius pc, arrive}: a sphere, framed to fill the view (`arrive` times that distance, 1 by default);
  // - {center, dist}: a point, from `dist` parsecs away;
  // - {entry}: an object of the sector scene or the map (a star, a point, a crowd {x, y, z, r}), from where it is looked at;
  // - {body: entry or ref}: a body of the open system, close enough to see it with its moons;
  // - {sector: {shell, band, wedge}, open, system, object}: a sector, opened in place with `open`, and with it its system (its printed
  //   ID) and a body of that (its `<kind>:<id>` reference);
  // - {region: {block, center, radius}, landing}: a charted region, landing "frame" (its charted part fills the view) or "densest" (the
  //   cluster of its heaviest sectors), the reader's setting `map.regionLanding` when `landing` is not given;
  // - {block}, {stage}, {up: true}, {home: true}: a block's window, a stage, one level out, the whole galaxy.
  // `label` (a name) puts the place in the visit log. True when the move began. Nothing here is held to a stage's range; only the
  // global limits hold (MAP.333), and the instant move under reduced motion is the flight's own.
  function goTo(target) {
    if (!target || !view) return false;
    let began = true;
    if (target.up) {
      up();
    } else if (target.home) {
      home();
    } else if (target.stage) {
      selection.select({ kind: "stage", stage: target.stage });
    } else if (target.region) {
      began = goToRegion(target.region, target.landing);
    } else if (target.block) {
      began = openBlock(target.block);
    } else if (target.sector) {
      pendingSystem = target.open && target.system ? target.system : null;
      pendingObject = pendingSystem ? target.object || null : null;
      began = target.open ? openSector(target.sector) : locate(target.sector);
    } else if (target.body) {
      const entry = typeof target.body === "string" ? systemStage.entry(target.body) : target.body;
      began = !!entry && !!systemStage.where(entry.ref);
      if (began) flyToBody(entry);
    } else if (target.entry) {
      flyToEntry(target.entry);
    } else if (target.center && target.dist > 0) {
      flyToPoint(target.center, target.dist);
    } else if (target.center && target.radius > 0) {
      flyToFit({ center: target.center, radius: target.radius }, target.arrive);
    } else {
      began = false;
    }
    if (began && target.label) recents.add(target.label, plainTarget(target));
    return began;
  }

  // The part of a target that is plain data, which the visit log keeps and goTo takes again.
  function plainTarget(target) {
    const plain = {};
    ["center", "radius", "dist", "arrive", "sector", "open", "system", "object", "home"].forEach(function (key) {
      if (target[key] !== undefined && target[key] !== null) plain[key] = target[key];
    });
    if (target.region) plain.region = { block: target.region.block, center: target.region.center, radius: target.region.radius };
    return plain;
  }

  // A region of the map: its window with the camera framing its extent, and, landing "densest", on from there to the cluster of its
  // heaviest sectors once their ranking is in (a region the window holds no ranking of stays framed).
  function goToRegion(region, landing) {
    if (!region.block) return false;
    const mode = landing || (host.regionLanding ? host.regionLanding() : "frame");
    const began = region.center && region.radius > 0 ? openRegion(region.block, { center: region.center, radius: region.radius }) : openBlock(region.block);
    landingAfter = began && mode === "densest" && region.center ? { center: region.center.slice(), radius: region.radius } : null;
    return began;
  }

  // The "densest" landing of a region the window has arrived at: the sectors round its middle are read (the near read, MAP.335, as far as
  // it looks: 4.5 edges, which holds a region as the map frames it) and the camera goes on to the cluster of the heaviest of them. A
  // region the read has no ranked sector for, or a map without the read, stays framed.
  const LANDING_READ_EDGES = 4.5;
  function landDensest() {
    const asked = landingAfter;
    landingAfter = null;
    if (!asked || pinned || !host.fetchNear) return;
    const token = ++landingToken;
    const at = asked.center;
    const query = "?x=" + at[0].toFixed(3) + "&y=" + at[1].toFixed(3) + "&z=" + at[2].toFixed(3) + "&r=" + (LANDING_READ_EDGES * edgePc).toFixed(3);
    host.fetchNear(query).then(function (payload) {
      if (!active || token !== landingToken || animation) return;
      const ranked = ((payload && payload.sectors) || []).map(function (sector) {
        const middle = P.cellCenter(P.drillBlockBounds({ m: 1, shell: sector.shell, band: sector.band, wedge: sector.wedge }, edgePc));
        const mass = sector.stats && sector.stats.mass_kg > 0 ? sector.stats.mass_kg : 0;
        return { x: middle[0], y: middle[1], z: middle[2], sizePc: edgePc, massKg: mass };
      });
      const cluster = RL.densestCluster(ranked, asked, edgePc);
      if (cluster) flyToFit({ center: cluster.center, radius: cluster.radius });
    }, function () {
      // The ranking is an aid; the region stays framed.
    });
  }

  // Flies the orbit center to `goal` (a world point), keeping the angle and
  // the zoom: turning and zooming then both go about it (MAP.199).
  function centreOn(goal) {
    dropZoomGoal();
    const from = viewNow();
    const to = { target: goal, dist: view.dist, quat: view.quat.clone() };
    followBody = null;
    const kept = view;
    flyCamera(from, to, null, function () {
      view = Object.assign({}, kept, { target: goal.slice(), dist: kept.dist, quat: kept.quat.clone() });
      applyView();
    });
    animation.soft = true;
  }

  // For the browser tests: the stage on screen.
  canvasEl.galaxyStage = function () { return stage; };
  // For the browser tests (MAP.337): the zoom level the camera is at, the window stage whose fit the camera's distance reaches (the one it would load
  // on from its address), `win.centre.m` being its block size; the whole-galaxy stage when it is out past them all.
  canvasEl.galaxyZoomLevel = function () {
    if (!view) return stage;
    const now = currentCamera();
    return stageForCamera({ target: now.target, dist: now.dist, quat: now.quat.toArray ? now.quat.toArray() : now.quat });
  };

  // For the browser tests: goTo (MAP.336) with a target as the stage view takes it; true when the move began.
  canvasEl.galaxyGoTo = function (target) { return goTo(target); };
  canvasEl.galaxyFov = function () { return camera.fov; };

  // For the browser tests: whether the camera is on its way to a stage, waiting for the stage's data or flying
  // (the trail follows the camera, MAP.227, so it is read once the camera has landed).
  canvasEl.galaxyMoving = function () { return arriving || !!animation; };

  // For the browser tests: opens the window of `g` cells round the cell {m, shell, band, wedge}, as the
  // breadcrumb's crumbs and a click on a block do.
  canvasEl.galaxyOpenWindow = function (cell, g) {
    if (view && !animation) flyToStage(S.windowStage(cell, g || P.DEFAULT_WINDOW));
  };

  // For the browser tests: moves the view's middle to a point of the map (parsecs).
  canvasEl.galaxyCentreOn = function (x, y, z) {
    if (view && !animation) centreOn([x, y, z]);
  };

  // MAP.199: Shift-click (or a right-click) takes the point under the
  // pointer for the new center.
  function centreAtPointer(event) {
    if (!view || animation || !resolved || !view.fit || !isFree(resolved)) return;
    const point = anchorUnder(event.clientX, event.clientY);
    if (!point) return;
    const goal = point.slice();
    centreOn(goal);
  }

  // The other layer (a phenomenon) under the pointer, while there is one.
  let hoveredOther = null;

  function hoverOther(found) {
    const layer = found && found.layer !== choiceLayer ? found.layer : null;
    if (hoveredOther && hoveredOther !== layer && hoveredOther.hover) hoveredOther.hover(null);
    hoveredOther = layer;
    if (layer && layer.hover) layer.hover(found.entry);
  }

  // MAP.298: the name under the pointer (the label engine's), lit and the pointer a hand; null when none.
  function labelUnder(event) {
    return host.labelAt ? host.labelAt(event.clientX, event.clientY) : null;
  }

  function hoverAt(event) {
    if (animation || event.pointerType === "touch") return;
    const label = labelUnder(event);
    if (host.labelHot) host.labelHot(label ? label.id : null);
    canvasEl.style.cursor = label ? "pointer" : "";
    if (label) {
      hoverOther(null);
      if (hover && !hover.sticky) setHover(null);
      showTooltip("", 0, 0);
      return;
    }
    const found = picker.pick(event.clientX, event.clientY, pointerLayers);
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
    // MAP.298: a click on a name selects its item and flies the camera to it.
    const label = labelUnder(event);
    if (label && host.labelAct) {
      showTooltip("", 0, 0);
      host.labelAct(label);
      return;
    }
    // Inside a container, a phenomenon under the click is shown rather
    // than the block picked. Stars never take the click (MAP.101): in a
    // dense sector they would hide it.
    const found = picker.pick(event.clientX, event.clientY, pointerLayers);
    if (found && found.layer !== choiceLayer) {
      showTooltip("", 0, 0);
      // Picking an object takes the light off a choice the arrow keys left lit, so Enter then goes to the object.
      if (hover && hover.sticky) setHover(null);
      if (found.layer.select) found.layer.select(found.entry);
      return;
    }
    if (host.clearSelection) host.clearSelection();
    const index = found ? found.entry : -1;
    if (index < 0) {
      // MAP.304: a click on nothing selects nothing and puts the panel away.
      if (!blocksPickable() && !sectorStage.deselect()) host.showHint(hintFor(resolved));
      return;
    }
    if (type === "touch" && touchPending !== index) {
      // First tap highlights, a second on the same thing acts.
      touchPending = index;
      setHover({ option: index, sticky: true });
      showTooltip(optionText(index), event.clientX, event.clientY);
      return;
    }
    touchPending = -1;
    showTooltip("", 0, 0);
    act(index, found && found.hit ? { x: found.hit.point.x, y: found.hit.point.y, z: found.hit.point.z } : null);
  }

  // MAP.245: a double click on empty space flies half way toward the point under the pointer (the zoom of a wheel
  // step, as far as it may go, in one flight); false when the view cannot be moved. The point stays where it is on
  // the screen, as it does for the wheel.
  function flyHalfWayTo(event) {
    if (!view || animation || !resolved || !view.fit || !isFree(resolved) || !MC.canZoom(zoomPolicyFor(resolved))) return false;
    const anchor = anchorUnder(event.clientX, event.clientY);
    if (!anchor) return false;
    const asked = limitedStep(anchor, 0.5);
    const zoom = MC.clampDistance(zoomLimits(), view.zoom * asked, 1);
    const taken = view.zoom > 0 ? zoom / view.zoom : 1;
    if (!(taken < 1)) return false;
    const target = MC.targetAboutAnchor(view.target, anchor, taken);
    const to = { target: target, dist: view.dist * taken, quat: view.quat.clone() };
    const kept = view;
    followBody = null;
    showTooltip("", 0, 0);
    flyCamera(viewNow(), to, null, function () {
      view = Object.assign({}, kept, { target: target.slice(), dist: to.dist, quat: to.quat.clone(), zoom: zoom });
      applyView();
    });
    animation.soft = true;
    return true;
  }

  // MAP.146: a double click on something of the scene (a star, a phenomenon, a body of the open system) flies the
  // camera to it. On anything else it is two clicks, as ever (a touch screen's second tap on a block acts on it).
  function doubleClickAt(event, type) {
    if (animation) return;
    if (labelUnder(event)) return;
    const found = picker.pick(event.clientX, event.clientY, pointerLayers);
    // MAP.245: on empty space (nothing picked) it flies half way toward the point under the pointer. A block stays
    // two clicks.
    if (!found && flyHalfWayTo(event)) return;
    if (!found || found.layer === choiceLayer) {
      clickAt(event, type);
      return;
    }
    showTooltip("", 0, 0);
    if (enteredSystem) {
      if (found.entry && found.entry.ref && systemStage.where(found.entry.ref)) {
        if (found.layer.select && systemStage.selected() !== found.entry.ref) found.layer.select(found.entry);
        flyToBody(found.entry);
      }
      return;
    }
    flyToEntry(found.entry);
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
      start: function () { pinchZoom = view.dist; },
      move: function (ratio) {
        if (!view) return;
        zoomToDistance(pinchZoom * ratio);
      },
    },
    onClick: clickAt,
    onDoubleClick: doubleClickAt,
    onPanClick: centreAtPointer,
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

  // Arrow keys move among the choices, through them in order; Enter takes
  // the one lit.
  function onKey(event) {
    const key = event.key;
    if (key === "Escape" || key === "Backspace") {
      event.preventDefault();
      // A selected object in the sector goes first, then a step out.
      if (sectorStage.deselect()) return;
      up();
      return;
    }
    if (key === "Home") {
      event.preventDefault();
      home();
      return;
    }
    if ((key === "[" || key === "]") && !event.ctrlKey && !event.metaKey && !event.altKey) {
      event.preventDefault();
      stepToNeighbour(key === "[" ? "in" : "out");
      return;
    }
    if ((key === "c" || key === "C") && !event.ctrlKey && !event.metaKey && !event.altKey) {
      event.preventDefault();
      recenter();
      return;
    }
    if (["ArrowLeft", "ArrowRight", "ArrowUp", "ArrowDown", "Enter"].indexOf(key) < 0) return;
    event.preventDefault();
    // Alt and an arrow move the center (MAP.200).
    // Alt, Shift and an arrow step to the neighbouring region at this level (MAP.121): up is north, down south, left
    // west, right east; [ and ] go in toward the core and out from it (below).
    if (event.altKey && event.shiftKey && key !== "Enter") {
      stepToNeighbour({ ArrowUp: "north", ArrowDown: "south", ArrowLeft: "west", ArrowRight: "east" }[key]);
      return;
    }
    if (event.altKey && key !== "Enter") {
      const steps = { ArrowLeft: [-KEY_MOVE, 0], ArrowRight: [KEY_MOVE, 0], ArrowUp: [0, KEY_MOVE], ArrowDown: [0, -KEY_MOVE] };
      moveCentre(steps[key][0], steps[key][1]);
      return;
    }
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
    if (animation) return;
    // MAP.297: there are no choices to move among on the Galaxy Map: Enter flies to the selected object.
    if (!blocksPickable()) {
      if (key === "Enter") flyToEntry(sectorStage.selectedEntry());
      return;
    }
    // MAP.146: with nothing lit among the choices, Enter flies to the selected object of the scene.
    if (key === "Enter" && !(hover && hover.option != null && display && display.options.length)) {
      flyToEntry(sectorStage.selectedEntry());
      return;
    }
    if (!display || !display.options.length) return;
    const options = display.options;
    const current = hover && hover.option != null ? hover.option : -1;
    if (key === "Enter") {
      if (current >= 0) act(current);
      return;
    }
    let next = 0;
    if (current >= 0) {
      const forward = key === "ArrowRight" || key === "ArrowUp";
      next = (current + (forward ? 1 : -1) + options.length) % options.length;
    }
    setHover({ option: next, sticky: true });
    tooltipAtOption(next);
  }

  // --- Breadcrumb, notice ----------------------------------------------------

  function notice(text) {
    if (!els.notice) return;
    els.notice.textContent = text || "";
    els.notice.hidden = !text;
  }

  // --- The camera in the address (MAP.228) -------------------------------------

  // The camera the address asked for, until the stage it names is on screen (flyTo takes it up).
  let pendingCamera = null;
  let camChangedAt = 0;
  const CAM_URL_DELAY_MS = 350;
  // A camera written to ten digits and read back is as far as this over or under the stage's own fit.
  const STAGE_FIT_TOLERANCE = 1e-6;

  // The view as {target, dist, quat: [x, y, z, w]}, or null.
  function currentCamera() {
    return view ? { target: view.target.slice(), dist: view.dist, quat: view.quat.toArray() } : null;
  }

  // The stage's own camera as the address holds it.
  function asAddress(camera) {
    return { target: camera.target.slice(), dist: camera.dist, quat: camera.quat.toArray ? camera.quat.toArray() : camera.quat.slice() };
  }

  // Whether `camera` is the stage's own: the galaxy's address then is none.
  function isStageCamera(r, camera) {
    const own = asAddress(cameraFor(r));
    const near = function (a, b, scale) { return Math.abs(a - b) <= 1e-6 * scale; };
    const dot = camera.quat.reduce(function (sum, value, k) { return sum + value * own.quat[k]; }, 0);
    return near(camera.dist, own.dist, own.dist) && Math.abs(Math.abs(dot) - 1) < 1e-9
      && camera.target.every(function (value, k) { return near(value, own.target[k], own.dist); });
  }

  // The query of an address to camera `camera` of stage `r`: "" for the whole galaxy at its own camera.
  function cameraQueryFor(r, camera) {
    const wanted = asAddress(camera);
    if (isWholeGalaxy(r) && isStageCamera(r, wanted)) return "";
    return "?" + S.cameraQuery(wanted);
  }

  function withoutCam(query) {
    const params = new URLSearchParams(query || "");
    params.delete("cam");
    const text = params.toString();
    return text ? "?" + text : "";
  }

  // The window stage that holds what a camera looks at: the finest level whose window, centred on the cell
  // holding the camera's middle, is no smaller than the view (the distance is no more than its own fit), else the
  // coarsest, and the whole galaxy when that is no smaller either. Null when the middle is outside the galaxy.
  function stageForCamera(camera) {
    const quat = new THREE.Quaternion().fromArray(camera.quat);
    const levels = P.STAGE_LEVELS.slice().reverse();
    const outline = getOutline();
    let coarsest = null;
    let coarsestFit = 0;
    for (let n = 0; n < levels.length; n++) {
      const cell = P.cellAt(levels[n], camera.target[0], camera.target[1], camera.target[2], edgePc);
      const candidate = S.windowStage(cell, P.DEFAULT_WINDOW);
      const r = S.resolveStage(candidate, outline, edgePc);
      if (r.problem) continue;
      const extent = S.extentOf(r.view.blocks);
      const fit = centeredFit(fitPointsOf(r), extent.center, quat).dist;
      if (fit * (1 + STAGE_FIT_TOLERANCE) >= camera.dist) return candidate;
      coarsest = candidate;
      coarsestFit = fit;
    }
    return coarsest && camera.dist <= MAX_ZOOM * coarsestFit ? coarsest : S.galaxyStage();
  }

  // Puts the view at `camera` within stage `r` (the view already at the stage's own camera): the distance within
  // the zoom's range.
  function viewAtCamera(r, camera) {
    if (!view || !view.fit) return;
    const quat = new THREE.Quaternion().fromArray(camera.quat);
    const full = fitDistance(fitPointsOf(r), view.fit.target, quat);
    if (!(full > 0)) return;
    view.quat = quat;
    view.target = camera.target.slice();
    // MAP.333: a camera put here (a flight's landing, an address) is held only by the global limits, not by the range of the stage
    // it lands in: the view is not pushed back out after a flight that came closer than the stage's nearest.
    view.zoom = Math.max(GLOBAL_FLOOR_PC, Math.min(globalFarPc(), camera.dist)) / full;
    view.dist = view.zoom * full;
  }

  // Writes the view's camera into the address in place (no history entry of its own) once it has been still for
  // a moment (`step` checks); a stage hop makes its own entry (`go`). An opened sector or system keeps its own address.
  function scheduleCameraAddress() {
    if (pinned || !active) return;
    camChangedAt = performance.now();
  }

  function writeCameraAddress() {
    camChangedAt = 0;
    if (!active || pinned || !view || animation || arriving || entered || enteredSystem || !resolved || !resolved.view) return;
    const params = new URLSearchParams(location.search);
    const query = cameraQueryFor(resolved, currentCamera());
    if (query) params.set("cam", query.slice(5));
    else params.delete("cam");
    // TEST.496: a sector left selected stays in the address while the camera is near it; out past SELECTION_REACH_EDGES edges the address is the camera's alone.
    if (params.has("sector") && view.dist > SELECTION_REACH_EDGES * edgePc) {
      params.delete("sector");
      params.delete("open");
    }
    const text = params.toString().replace(/%2C/gi, ",").replace(/%3A/gi, ":");
    const next = text ? "?" + text : "";
    if (next !== location.search) history.replaceState(history.state, "", location.pathname + next + location.hash);
  }

  // MAP.227: the stage the camera is in, worked out from where the view's middle is: the breadcrumb is the trail
  // of containers of the cell holding it (galaxy, blocks, sector), so panning to a neighbouring cell changes it and
  // "Sector X" in the header names the sector it is in, with no stored state of its own. An opened sector or
  // system is shown in the window round that sector.
  function cameraStage() {
    if (!view || !stage.win) return stage;
    if (entered) {
      // The address of an open sector names only the sector, so the trail is the one of the window round it, as a
      // fresh load of that address shows.
      return S.sectorStage(entered.sector.shell, entered.sector.band, entered.sector.wedge) || stage;
    }
    if (enteredSystem) return stage;
    return S.stageAtPoint(stage, view.target, getOutline(), edgePc);
  }

  let whereKey = "";

  // Redraws the breadcrumb and the header line when the view's middle has moved into another cell or sector.
  function refreshWhere() {
    if (!active || !view || !resolved) return;
    const sector = entered || enteredSystem ? null : S.sectorAtPoint(view.target, edgePc);
    const win = cameraStage().win;
    const key = (win ? win.centre.m + ":" + S.blockLabel(win.centre) : "galaxy") + "|" + (sector ? S.blockLabel(sector) : "");
    if (key === whereKey) return;
    whereKey = key;
    if (els.where) {
      els.where.textContent = sector ? "In Sector " + S.blockLabel(sector) : "";
      els.where.hidden = !sector;
    }
    renderCrumbs();
  }

  // The breadcrumb's ☆ (MAP.23, design doc section 8.2): saves the
  // selected sector, else the stage's own URL, in bookmarks.js.
  let bookmarkButton = null;
  let refreshBookmark = null;

  function bookmarkEntry() {
    if (!resolved) return null;
    if (selectedSector) {
      const designation = S.sectorDesignation(selectedSector.shell, selectedSector.band, selectedSector.wedge);
      const data = dataFor(S.holder(stage));
      const sector = data ? data.sectors.get(cellKey(selectedSector)) : null;
      // The sector page itself, without a pick mode's query.
      const page = sector && host.sectorUrl(sector.id) ? host.sectorUrl(sector.id).split("?")[0] : null;
      return {
        kind: "sector", value: designation,
        name: sector ? sector.name : "Sector " + S.blockLabel(selectedSector),
        url: page || location.pathname + "?sector=" + encodeURIComponent(designation),
        sectorId: sector ? sector.id : null,
      };
    }
    const labels = S.crumbs(cameraStage(), getOutline(), edgePc).map(function (crumb) { return crumb.label; });
    return { kind: "stage", value: location.pathname + displayQuery(resolved ? cameraQueryFor(resolved, currentCamera() || cameraFor(resolved)) : ""), name: labels.join(" › ") };
  }

  // The breadcrumb (MAP.93): one line at any width. When the steps don't
  // all fit, it shows the first, a "…" button (a menu of the steps it
  // hides) and as many of the last ones as fit, then the current one,
  // measured again whenever the line's width changes. At phone width the
  // line gives way to the round Steps button between Back and Forward
  // (MAP.94, els.steps), a menu of every step.
  // UX.59: on a page whose own crumb is just Home > Galaxy, the map's crumb
  // takes its place and starts at Home, so there is one trail.
  const homeStep = (function () {
    const nav = els.pageCrumb;
    const rows = nav && els.crumbs ? nav.querySelectorAll("ol > li") : [];
    const link = rows.length === 2 ? rows[0].querySelector("a") : null;
    return link ? { label: link.textContent.trim(), href: link.getAttribute("href") } : null;
  })();

  function crumbSteps() {
    const shown = cameraStage();
    const items = S.crumbs(shown, getOutline(), edgePc);
    if (homeStep) items.unshift({ label: homeStep.label, href: homeStep.href, last: false });
    if (selectedSector && !isSelectedCentre(entered ? shown : stage)) {
      const uncharted = entered && entered.uncharted && !pinned;
      items.push({
        label: "Sector " + S.blockLabel(selectedSector)
          + (uncharted ? " (uncharted)" : ""),
        last: true, sector: true,
      });
      items[items.length - 2].last = false;
    } else if (selectedSector && entered && entered.uncharted && !pinned) {
      // The window is centred on the sector: its own step carries the mark.
      items[items.length - 1].label += " (uncharted)";
    }
    if (enteredSystem) {
      // Up from the system is to the sector it sits in.
      items[items.length - 1].last = false;
      items[items.length - 1].systemUp = true;
      items.push({ label: enteredSystem.name, last: true, system: true });
    }
    return items;
  }

  // The line itself is the shared breadcrumb (breadcrumb.js, NAV.14); its
  // steps are the stage's, taken through the selection.
  // UX.58: "Jump to newest" (MAP.95) is the last item of the Steps menu.
  const newestButton = document.createElement("button");
  newestButton.type = "button";
  newestButton.className = "galaxy-steps-newest";
  newestButton.dataset.action = "current";
  newestButton.textContent = "Jump to newest";
  newestButton.title = "Jump to the newest view in this history";
  newestButton.disabled = true;
  newestButton.addEventListener("click", function () {
    if (els.steps) els.steps.open = false;
    travelToCurrent();
  });
  const crumbLine = els.crumbs ? createBreadcrumb(els.crumbs, {
    steps: els.steps,
    stepsExtra: newestButton,
    trailing: function () { return bookmarkButton; },
    // MAP.336: each crumb is a Go to: the container it names (the stage), or one level out from an open system.
    onSelect: function (step) {
      if (step.systemUp) goTo({ up: true });
      else if (step.stage) goTo({ stage: step.stage });
    },
  }) : null;

  function renderCrumbs() {
    syncSelection();
    if (!crumbLine) return;
    if (!bookmarkButton) {
      bookmarkButton = document.createElement("button");
      bookmarkButton.type = "button";
      bookmarkButton.className = "btn btn-small btn-secondary btn-bookmark galaxy-bookmark";
      refreshBookmark = B.toggleButton(bookmarkButton, bookmarkEntry);
    }
    crumbLine.set(crumbSteps());
    refreshBookmark();
  }

  // A resized map (the window resized or turned) refits the view to its
  // new size (MAP.53).
  function relayout() {
    if (active && view && !animation) {
      reframe();
      applyView();
    }
  }
  if (typeof ResizeObserver === "function") {
    const observer = new ResizeObserver(relayout);
    observer.observe(canvasEl);
  }
  window.addEventListener("resize", relayout);

  // --- The address bar (section 9.3) -----------------------------------------

  function locate(sector) {
    const next = S.sectorStage(sector.shell, sector.band, sector.wedge);
    if (!next) {
      notice("Sector " + S.blockLabel(sector) + " is outside the galaxy.");
      return false;
    }
    selectedSector = { shell: sector.shell, band: sector.band, wedge: sector.wedge };
    flyToStage(next, { keepSector: true, query: "?sector=" + S.sectorDesignation(sector.shell, sector.band, sector.wedge) });
    return true;
  }

  // MAP.327: the icon on a sector's wall: the camera goes to the neighbouring sector's window and the sector opens in place
  // there (as the address ?sector=...&open=1 does), the address naming it.
  function openSector(sector) {
    const next = S.sectorStage(sector.shell, sector.band, sector.wedge);
    if (!next) {
      notice("Sector " + S.blockLabel(sector) + " is outside the galaxy.");
      return false;
    }
    selectedSector = { shell: sector.shell, band: sector.band, wedge: sector.wedge };
    pendingOpen = selectedSector;
    flyToStage(next, { keepSector: true, query: sectorQuery(selectedSector, true) });
    return true;
  }

  // --- The find-a-place box (MAP.336) ------------------------------------------

  // The places gone to by name, kept in the browser (findplace.js).
  let recentsStorage = null;
  try { recentsStorage = globalThis.localStorage || null; } catch (error) { recentsStorage = null; }
  const recents = FPLACE.createRecents(recentsStorage);

  const findBox = els.address && els.matches ? createFindBox(els.address, els.address.querySelector("input"), els.matches, {
    parse: function (text) { return S.parseAddress(text, edgePc); },
    parseSector: S.parseSectorDesignation,
    sectorExists: function (sector) { return !!S.sectorStage(sector.shell, sector.band, sector.wedge); },
    regions: host.regions || function () { return []; },
    recents: recents,
    locate: function (name) { return host.locate(name); },
    home: { label: "Whole galaxy", target: { home: true } },
    go: goTo,
    notice: notice,
  }) : null;
  if (els.address) {
    els.address.addEventListener("submit", function (event) {
      event.preventDefault();
      if (findBox) findBox.submit();
    });
  }
  function clearMatches() {
    if (findBox) findBox.close();
  }

  // --- Turning on, and the URL -----------------------------------------------

  function setActive(on) {
    if (on === active) return;
    active = on;
    [els.crumbs, els.address].forEach(function (el) { if (el) el.hidden = !on; });
    if (!on) {
      whereKey = "";
      if (els.where) els.where.hidden = true;
    }
    if (homeStep) {
      els.pageCrumb.hidden = on;
      if (on) els.pageCrumb.parentNode.insertBefore(els.crumbs, els.pageCrumb);
    }
    if (!on) {
      clearMatches();
      if (animation) finishAnimation();
      closeSector();
      pendingOpen = null;
      pendingSystem = null;
      disposeDisplay(display);
      display = null;
      view = null;
      clearOutline();
      showTooltip("", 0, 0);
      notice("");
      return;
    }
    const r = resolved || resolve(stage);
    loadData(S.holder(r.stage)).then(function () {
      if (!active) return;
      display = null;
      flyTo(r);
    });
  }

  // The stage a URL's query asks for: {stage, sector, problem}. ?sector=
  // opens the stage showing that sector; with a NAV course and no stage
  // asked for, the smallest stage showing the whole course (section 9.4).
  function stageFromLocation(useCourse) {
    if (pinned) {
      const next = S.sectorStage(pinned.shell, pinned.band, pinned.wedge);
      if (!next) return { stage: S.galaxyStage(), sector: null, problem: "This sector is outside the galaxy." };
      return { stage: next, sector: pinned, open: true, problem: null };
    }
    const parsed = S.parseStageQuery(location.search);
    if (parsed.problem && !parsed.camera && !parsed.sector) return { stage: S.galaxyStage(), sector: null, problem: parsed.problem };
    // The stage the camera asked for is the one whose window holds what it looks at (MAP.228, MAP.337): the camera in the address
    // is all a history entry keeps, so Back, Forward and a reload work the stage out from it.
    const asked = parsed.camera ? stageForCamera(parsed.camera) : null;
    if (parsed.sector) {
      const s = parsed.sector;
      const next = asked || S.sectorStage(s.shell, s.band, s.wedge);
      if (!next) {
        return { stage: S.galaxyStage(), sector: null, problem: "Sector " + S.blockLabel(s) + " is outside the galaxy." };
      }
      return { stage: next, camera: asked ? parsed.camera : null, sector: s, open: parsed.open, system: parsed.system, object: parsed.object, problem: null };
    }
    if (asked) return { stage: asked, camera: parsed.camera, sector: null, problem: null };
    let next = S.galaxyStage();
    if (useCourse && host.courseSectors && host.courseSectors.length) {
      next = S.courseStage(host.courseSectors, getOutline(), edgePc);
    }
    return { stage: next, sector: null, problem: parsed.problem };
  }

  // Opens the stage the URL asks for (or the galaxy, with a notice).
  function openFromLocation() {
    const state = history.state && history.state.galaxyStage ? history.state : null;
    mapIndex = state && state.mapIndex != null ? state.mapIndex : 0;
    maxIndex = Math.max(mapIndex, state && state.maxIndex != null ? state.maxIndex : mapIndex);
    if (!pinned) history.replaceState({ galaxyStage: true, mapIndex: mapIndex, maxIndex: maxIndex }, "");
    // The marks, star classes and floor first: taking the hidden kinds up rewrites the URL
    // from what is set so far, and would drop a `mark` that wasn't set yet.
    applyStarFiltersFromLocation();
    takeLayersFromLocation();
    setColorFromLocation();
    const asked = stageFromLocation(true);
    let r = resolve(asked.stage);
    let problem = asked.problem || r.problem;
    if (r.problem) r = resolve(S.galaxyStage());
    pendingCamera = r.problem ? null : asked.camera || null;
    selectedSector = asked.sector;
    pendingOpen = asked.open ? asked.sector : null;
    pendingSystem = asked.open && asked.system ? asked.system : null;
    pendingObject = pendingSystem ? asked.object : null;
    stage = r.stage;
    resolved = r;
    renderCrumbs();
    renderTravel();
    notice(problem || "");
  }

  function onPopState(event) {
    if (!active || pinned) return;
    // MAP.352: a move a quick press before this one began is finished first, with what it was asked to open, so what this entry asks
    // to open is not handed to the display that move was arriving at.
    if (animation) finishAnimation();
    const state = event.state && event.state.galaxyStage ? event.state : null;
    mapIndex = state && state.mapIndex != null ? state.mapIndex : 0;
    maxIndex = Math.max(maxIndex, mapIndex);
    // Remembered on this entry too, so a reload here still knows how far
    // Forward goes.
    if (state) history.replaceState({ galaxyStage: true, mapIndex: mapIndex, maxIndex: maxIndex }, "");
    // The marks, star classes and floor first: taking the hidden kinds up rewrites the URL
    // from what is set so far, and would drop a `mark` that wasn't set yet.
    applyStarFiltersFromLocation();
    takeLayersFromLocation();
    setColorFromLocation();
    const asked = stageFromLocation(false);
    if (asked.problem) {
      selectedSector = null;
      flyToStage(S.galaxyStage(), { push: false });
      notice(asked.problem);
      return;
    }
    selectedSector = asked.sector;
    pendingOpen = asked.open ? asked.sector : null;
    pendingSystem = asked.open && asked.system ? asked.system : null;
    pendingObject = pendingSystem ? asked.object : null;
    flyToStage(asked.stage, { push: false, keepSector: true, camera: asked.camera || null });
  }
  window.addEventListener("popstate", onPopState);

  function setNeedGenerated(on) {
    needGenerated = on;
    if (display && !animation) rebuildDisplay();
  }

  function setChartedOnly(on) {
    chartedOnly = on;
    if (display && !animation) rebuildDisplay();
  }

  // Whether the blocks are drawn: a sector page's map (pinned) as it always was; else only with the cell grid tool, a
  // fill chosen by "Color by", or "Charted only" (MAP.265).
  function cellsShown() {
    return pinned || cellGrid || chartedOnly || colorBy !== "default";
  }

  function applyCellsShown() {
    [display].concat(leaving).forEach(function (d) {
      if (!d) return;
      d.groups.forEach(function (group) {
        group.meshes.forEach(function (mesh) { mesh.visible = cellsShown(); });
      });
    });
  }

  // MAP.265: the optional grid tool, the cells of the stage drawn (a toggle of the map's Menu).
  function setCellGrid(on) {
    cellGrid = !!on;
    applyCellsShown();
  }
  canvasEl.galaxyCellGrid = function (on) {
    if (on !== undefined) setCellGrid(on);
    return cellGrid;
  };

  // MAP.131: what the blocks' fill shows; the choice rides in the URL in
  // place (it is no step of its own).
  function setColorBy(mode) {
    if (mode === colorBy) return;
    colorBy = mode;
    applyCellsShown();
    const params = new URLSearchParams(location.search);
    if (mode === "default") params.delete("color");
    else params.set("color", mode);
    const query = params.toString().replace(/%2C/gi, ",").replace(/%3A/gi, ":");
    if (!pinned) history.replaceState(history.state, "", location.pathname + (query ? "?" + query : "") + location.hash);
    if (display && !animation) rebuildDisplay();
  }

  // MAP.268: what the label engine (maplabels.js) may name here. The sectors of the stage's window that hold a
  // generated sector, by name, each with its edge (parsecs): {id, text, x, y, z, sizePc}.
  function sectorLabels() {
    if (!display || !display.options || entered || pinned) return [];
    const labels = [];
    display.options.forEach(function (option) {
      const block = option.blocks[0];
      if (option.blocks.length !== 1 || block.m !== 1) return;
      const record = sectorRecord(block);
      if (!record) return;
      const middle = P.cellCenter(block.bounds);
      labels.push({
        id: "sector:" + (record.id !== undefined && record.id !== null ? record.id : cellKey(block)),
        text: record.name || "Sector " + S.sectorDesignation(block.shell, block.band, block.wedge),
        x: middle[0], y: middle[1], z: middle[2], sizePc: edgePc,
        slice: { massKg: record.stats && record.stats.mass_kg > 0 ? record.stats.mass_kg : 0 },    // MAP.322: the heaviest sectors are named first
        ref: { kind: "sector", sector: { shell: block.shell, band: block.band, wedge: block.wedge } },
      });
    });
    return labels;
  }

  // ... and the named systems of the sector the map has open or is coming up to, brightest first by the caller:
  // {id, text, x, y, z, luminositySol}.
  function systemLabels() {
    if (enteredSystem) return [];
    return sectorStage.entries().filter(function (e) {
      return e.name && e.href && e.endpoint && e.endpoint.indexOf("system:") === 0 && !sectorStage.hides(e);
    }).map(function (e) {
      return { id: e.endpoint, text: e.name, x: e.x, y: e.y, z: e.z, luminositySol: e.luminositySol || 0, ref: { kind: "system", entry: e } };
    });
  }

  // ... and the bodies of the system open in place or being come up to: {id, text, x, y, z, sizePc} (the width of
  // the body in parsecs, so its label shows once it is a few pixels across).
  function bodyLabels() {
    return systemStage.entries().filter(function (e) { return e.name && e.r > 0; }).map(function (e) {
      return { id: e.ref, text: e.name, x: e.x, y: e.y, z: e.z, sizePc: 2 * e.r, slice: bodySlice(e), ref: { kind: "body", entry: e } };
    });
  }

  // The label budget's view of a body (MAP.317): a star among the stars (by size, in solar radii: 1 R_sun is 2.2546e-8 pc),
  // a gas giant, a moon, a comet or another world in a category of the uninhabited slice (by size in parsecs across).
  const SOLAR_RADIUS_PC = 2.2546e-8;
  function bodySlice(entry) {
    if (entry.kind === "star") return { domain: "star", radiusSol: entry.r / SOLAR_RADIUS_PC, luminositySol: (entry.body && entry.body.luminosity_sol) || 0 };
    const look = entry.body && entry.body.look;
    const category = look && look.body === "gas" ? "gas giant" : entry.kind === "moon" ? "moon" : entry.kind === "comet" ? "comet" : "rocky or icy world";
    return { domain: "uninhabited", category: category, size: 2 * entry.r };
  }

  // MAP.298: a click on a body's name selects the body and flies to it.
  function selectBody(entry) {
    if (!entry) return;
    systemStage.select(entry.ref);
    flyToBody(entry);
  }

  // MAP.321: the open system for the body symbols: {id, centre, radiusPc, entries (the bodies, kept up to date as they move)},
  // or null when none is open in place.
  function systemMarkerSource() {
    if (!enteredSystem || !systemStage.isOpen()) return null;
    return { id: systemStage.systemId(), centre: systemStage.centre(), radiusPc: systemStage.radiusPc(), entries: systemStage.markerBodies() };
  }

  // MAP.338: the open system's group and place (parsecs), or null when none is open in place.
  function systemFrame() { return enteredSystem ? systemStage.frame() : null; }

  function selectBodyRef(ref) { systemStage.select(ref); }

  return {
    selectBody: selectBody,
    selectBodyRef: selectBodyRef,
    systemMarkerSource: systemMarkerSource,
    systemFrame: systemFrame,
    bodyLabels: bodyLabels,
    sectorLabels: sectorLabels,
    systemLabels: systemLabels,
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
    recenter: recenter,
    home: home,
    up: up,
    stepToNeighbour: stepToNeighbour,
    openBlock: openBlock,
    openRegion: openRegion,
    canStep: function (direction) { return !!stepTarget(direction); },
    travel: travel,
    travelToCurrent: travelToCurrent,
    setChartedOnly: setChartedOnly,
    setCellGrid: setCellGrid,
    setColorBy: setColorBy,
    colorBy: function () { return colorBy; },
    setNeedGenerated: setNeedGenerated,
    stage: function () { return stage; },
    // MAP.336: the one way to go to anything (see goTo).
    goTo: goTo,
    // The opened sector's entries (MAP.68): "Show on map" buttons select
    // one by key.
    // The world point under a screen point (MAP.150); read by the tests.
    anchorUnder: anchorUnder,
    zoomIn: function () { zoomButton(1 / BUTTON_ZOOM); },
    zoomOut: function () { zoomButton(BUTTON_ZOOM); },
    // MAP.200: the on-map controls. The zoom as a share (0 farthest, 1
    // nearest) along the policy's range, on a log scale; a step of the
    // center in the screen's plane (`dx`, `dy` in shares of the map's
    // width, +x right, +y up); the center back to the middle of the region
    // viewed; and where the center is (world parsecs).
    zoomShare: zoomShare,
    setZoomShare: setZoomShare,
    // MAP.336: the log ruler's facts: the camera's distance from its focus, the farthest and nearest the map goes (pc), and the nearest named
    // thing of the scene up ({name, kind, distance}), or null with none.
    scaleInfo: function () {
      if (!view || !(view.dist > 0)) return null;
      return { dist: view.dist, farPc: globalFarPc(), floorPc: GLOBAL_FLOOR_PC, nearest: nearestNamed() };
    },
    moveCentre: moveCentre,
    resetCentre: resetCentre,
    centre: function () { return view ? view.target.slice() : null; },
    controlsEnabled: function () { return !!(view && MC.canZoom(zoomPolicyFor(resolved)) && isFree(resolved)); },
    kinds: sectorStage.kinds,
    kindHidden: sectorStage.kindHidden,
    setKindHidden: setKindHidden,
    // MAP.217: the layer registry (the Menu draws its buttons from it).
    layers: layers,
    // MAP.328: the Spin switch's standing preference, for the Menu's Spin layer (the galaxy map itself turns nothing).
    spinSwitch: systemStage.spinSwitch,
    readLayersFromLocation: readLayersFromLocation,
    kindMarked: sectorStage.kindMarked,
    setKindMarked: setKindMarked,
    // MAP.123: star classes and the luminosity floor.
    starClasses: sectorStage.starClasses,
    sectorLuminosityRange: sectorStage.luminosityRange,
    hiddenClasses: sectorStage.hiddenClasses,
    starClassHidden: sectorStage.starClassHidden,
    setStarClassHidden: setStarClassHidden,
    minLuminosity: sectorStage.minLuminosity,
    setMinLuminosity: setMinLuminosity,
    pickChanged: pickChanged,
    entryByKey: sectorStage.entryByKey,
    selectEntry: sectorStage.select,
    // Inside a window (not the whole galaxy), where phenomena take the
    // pointer before the blocks.
    inContainer: function () { return !!stage.win; },
    outsideOutline: outsideOutline,
    flyToStage: flyToStage,
  };
}
