// html/static/galaxylabels.js
//
// The Galaxy Map's names (MAP.268, MAP.229 to MAP.233): the three kinds of name the map writes through the one label
// engine (maplabels.js) from what it already holds: sectors (a generated sector in view, by name), systems (the
// named stars of the sector open or being come up to), the bodies of the system open (a planet or moon once it is
// a few pixels across, from the far view of its sector down to the planet itself) and nebulae (their name, else their type and designation,
// where the nebula is large enough on screen to read). Each kind is a layer of the registry (maplayers.js), so the
// Menu's "Labels" group, the URL's `hide` list and the legend come from the same declaration as every other layer.
//
// What a kind names comes from a `sources` object the map gives: sectors(), systems() and nebulae(), each a list of
// plain objects, so this file knows nothing of tiles or stages. A source anchors a label at x, y, z and may say how
// large the thing is (`sizePc`, its width in parsecs, turned into `pixels` from the camera) or what its pixels
// are already (`pixels`).
//
// The kinds are ranked so the big things win the room: a region over a nebula over a body over a sector over a system, and a system by
// its apparent brightness (luminosity over distance squared). No DOM at import; three.js comes in as a parameter.

const VERSION_QUERY = new URL(import.meta.url).search;

const { createLabelLayer, inRange } = await import(`./maplabels.js${VERSION_QUERY}`);
const { chooseLabels, labelBudget } = await import(`./labelbudget.js${VERSION_QUERY}`);

export const LABEL_KINDS = [
  { id: "label-region", kind: "region", label: "Region names", title: "Show or hide the names of charted regions and the families of regions" },
  { id: "label-sector", kind: "sector", label: "Sector names", title: "Show or hide the names of generated sectors" },
  { id: "label-system", kind: "system", label: "System names", title: "Show or hide the names of the stars of the open sector" },
  { id: "label-body", kind: "body", label: "Body names", title: "Show or hide the names of the planets and moons of an open system" },
  { id: "label-star", kind: "star", label: "Star names", title: "Show or hide the names of the named stars in view" },
  { id: "label-object", kind: "object", label: "Object names", title: "Show or hide the names of the other named objects in view: comets, rogue planets, black holes, neutron stars and the rest" },
  { id: "label-world", kind: "world", label: "World names", title: "Show or hide the names of the habitable-world markers the map draws" },
  { id: "label-nebula", kind: "nebula", label: "Nebula names", title: "Show or hide the names of nebulae and remnants" },
];

// The size (pixels across) a thing needs before its name is shown, and the width past which its name goes (the
// camera is inside it, or it fills the view).
export const MIN_PIXELS = { sector: 36, nebula: 48, body: 6 };
export const MAX_PIXELS = { sector: 900, nebula: 1400 };
// The share of the label budget (labelbudget.js, MAP.317) a kind of name draws on when its source does not say: a source
// may give an item its own `slice` ({domain, category, and the domain's numbers}). A sector stands for its region and
// is named first (the coarse structure), the rest take their place by share.
const DOMAIN = { region: "region", nebula: "nebula", body: "uninhabited", sector: "region", system: "star", star: "star", object: "uninhabited", world: "inhabited" };
// Priority of a name chosen by the budget: the chosen first-ranked takes FIRST_PRIORITY and each next one a step less.
const FIRST_PRIORITY = 10000;
const REFRESH_MS = 120;

export function createGalaxyLabels(options) {
  const THREE = options.THREE;
  const camera = options.camera;
  const canvasEl = options.canvasEl;
  const sources = options.sources;
  const layers = options.layers;
  const container = document.createElement("div");
  container.className = "map-labels";
  container.setAttribute("aria-hidden", "true");
  if (options.light) container.dataset.light = "";
  options.parent.appendChild(container);
  const point = new THREE.Vector3();
  const reduced = typeof matchMedia === "function" ? matchMedia("(prefers-reduced-motion: reduce)") : null;

  function size() {
    return { width: canvasEl.clientWidth || 1, height: canvasEl.clientHeight || 1 };
  }

  // The anchor's place on the map in CSS pixels and its depth, or null behind the camera.
  function project(x, y, z) {
    point.set(x, y, z);
    point.applyMatrix4(camera.matrixWorldInverse);
    if (point.z >= 0) return null;
    const depth = -point.z;
    point.set(x, y, z).project(camera);
    const view = size();
    return { x: (point.x + 1) / 2 * view.width, y: (1 - point.y) / 2 * view.height, depth: depth };
  }

  const layer = createLabelLayer({
    container: container, project: project, viewport: size, budget: options.budget,
    reducedMotion: function () { return !!(reduced && reduced.matches); },
  });

  // Pixels across for a thing `sizePc` wide at `depth` from the camera.
  function pixelsFor(sizePc, depth) {
    const view = size();
    return sizePc * view.height / (2 * Math.tan(THREE.MathUtils.degToRad(camera.fov) / 2) * Math.max(depth, 1e-9));
  }

  // What the label engine is offered. Every source's things that are in view and in their size range are gathered with
  // the domain they draw their share from (DOMAIN, or the item's own `slice`), the budget (screen area over a name's
  // area, labelbudget.js) chooses which of them to offer, and a chosen one's priority is its rank, so the engine
  // places the first-ranked first and none overlap (MAP.317).
  function candidates() {
    const found = [];
    camera.updateMatrixWorld();
    camera.matrixWorldInverse.copy(camera.matrixWorld).invert();
    const view = size();
    function depthOf(item) {
      return camera.position.distanceTo(point.set(item.x, item.y, item.z));
    }
    function offer(kind, items) {
      items.forEach(function (item) {
        const depth = depthOf(item);
        const pixels = item.pixels !== undefined ? item.pixels : item.sizePc ? pixelsFor(item.sizePc, depth) : undefined;
        const candidate = {
          id: item.id, kind: kind, text: item.text, x: item.x, y: item.y, z: item.z, ref: item.ref,
          pixels: pixels, minPixels: MIN_PIXELS[kind], maxPixels: MAX_PIXELS[kind],
        };
        if (!inRange(candidate, true)) return;
        const spot = project(item.x, item.y, item.z);
        if (!spot || !(spot.x >= 0 && spot.x <= view.width && spot.y >= 0 && spot.y <= view.height)) return;
        const slice = Object.assign({ domain: DOMAIN[kind], first: kind === "sector" }, defaultsOf(kind, item, depth, pixels), item.slice);
        found.push({ candidate: candidate, slice: slice });
      });
    }
    // The numbers a kind's things are ranked by when the source gives no `slice` of its own.
    function defaultsOf(kind, item, depth, pixels) {
      if (kind === "system" || kind === "star") return { luminositySol: item.luminositySol || 0 };
      if (kind === "nebula") return { volume: Math.pow(Math.max(pixels || 0, 1) * depth, 3) };
      if (kind === "object") return { category: "star-like", mass: item.luminositySol || 0 };
      if (kind === "body") return { category: "world", size: item.sizePc || 0 };
      if (kind === "sector") return { shareCharted: 0 };
      return {};
    }
    if (sources.regions && layers.isOn("label-region")) offer("region", sources.regions());
    if (layers.isOn("label-nebula")) offer("nebula", sources.nebulae());
    if (layers.isOn("label-body")) offer("body", sources.bodies());
    if (layers.isOn("label-sector")) offer("sector", sources.sectors());
    if (layers.isOn("label-system")) offer("system", sources.systems());
    if (sources.stars && layers.isOn("label-star")) offer("star", sources.stars());
    if (sources.objects && layers.isOn("label-object")) offer("object", sources.objects());
    if (sources.worlds && layers.isOn("label-world")) offer("world", sources.worlds());
    const out = [];
    const items = found.map(function (entry, n) {
      return Object.assign({ index: n, id: entry.candidate.id }, entry.slice);
    });
    chooseLabels(items, { budget: labelBudget(view) }).forEach(function (chosen) {
      const candidate = found[chosen.index].candidate;
      candidate.priority = FIRST_PRIORITY - chosen.rank;
      out.push(candidate);
    });
    return out;
  }

  let last = 0;
  let seen = null;
  // True when the camera is not where it was at the last placement.
  function cameraMoved() {
    const elements = camera.matrixWorld.elements.concat([camera.projectionMatrix.elements[0], camera.projectionMatrix.elements[5]]);
    const moved = !seen || elements.some(function (value, n) { return value !== seen[n]; });
    seen = elements;
    return moved;
  }

  // Called every frame: the names are placed again whenever the camera has moved, else (the data may have changed)
  // at most every REFRESH_MS.
  function update(now, force) {
    camera.updateMatrixWorld();
    const moved = cameraMoved();
    if (!force && !moved && now - last < REFRESH_MS) return;
    last = now;
    layer.update(candidates());
  }

  return { update: update, shown: layer.shown, element: container, clear: layer.clear, at: layer.at, setHot: layer.setHot };
}
