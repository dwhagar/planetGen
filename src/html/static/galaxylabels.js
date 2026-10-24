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
// The kinds are ranked so the big things win the room: a nebula over a body over a sector over a system, and a system by
// its apparent brightness (luminosity over distance squared). No DOM at import; three.js comes in as a parameter.

const VERSION_QUERY = new URL(import.meta.url).search;

const { createLabelLayer } = await import(`./maplabels.js${VERSION_QUERY}`);

export const LABEL_KINDS = [
  { id: "label-sector", kind: "sector", label: "Sector names", title: "Show or hide the names of generated sectors" },
  { id: "label-system", kind: "system", label: "System names", title: "Show or hide the names of the stars of the open sector" },
  { id: "label-body", kind: "body", label: "Body names", title: "Show or hide the names of the planets and moons of an open system" },
  { id: "label-nebula", kind: "nebula", label: "Nebula names", title: "Show or hide the names of nebulae and remnants" },
];

// The size (pixels across) a thing needs before its name is shown, and the width past which its name goes (the
// camera is inside it, or it fills the view).
export const MIN_PIXELS = { sector: 36, nebula: 48, body: 6 };
export const MAX_PIXELS = { sector: 900, nebula: 1400 };
const TIER = { nebula: 200, body: 150, sector: 100, system: 0 };
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

  function candidates() {
    const found = [];
    camera.updateMatrixWorld();
    camera.matrixWorldInverse.copy(camera.matrixWorld).invert();
    function ranged(kind, item, priority, depth) {
      const pixels = item.pixels !== undefined ? item.pixels : item.sizePc ? pixelsFor(item.sizePc, depth) : undefined;
      found.push({
        id: item.id, kind: kind, text: item.text, x: item.x, y: item.y, z: item.z, priority: TIER[kind] + priority,
        pixels: pixels, minPixels: MIN_PIXELS[kind], maxPixels: MAX_PIXELS[kind],
      });
    }
    function depthOf(item) {
      return camera.position.distanceTo(point.set(item.x, item.y, item.z));
    }
    if (layers.isOn("label-nebula")) {
      sources.nebulae().forEach(function (item) { ranged("nebula", item, Math.log10(Math.max(item.pixels || 1, 1)), depthOf(item)); });
    }
    if (layers.isOn("label-body")) {
      sources.bodies().forEach(function (item) {
        const depth = depthOf(item);
        ranged("body", item, Math.min(49, Math.log10(Math.max(pixelsFor(item.sizePc, depth), 1)) * 10), depth);
      });
    }
    if (layers.isOn("label-sector")) {
      sources.sectors().forEach(function (item) { ranged("sector", item, 0, depthOf(item)); });
    }
    if (layers.isOn("label-system")) {
      sources.systems().forEach(function (item) {
        const depth = Math.max(depthOf(item), 0.01);
        // Apparent brightness on a log scale, 0 to 99 within the tier.
        const priority = Math.max(0, Math.min(99, 50 + 10 * (Math.log10(Math.max(item.luminositySol || 0, 1e-6)) - 2 * Math.log10(depth))));
        ranged("system", item, priority, depth);
      });
    }
    return found;
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

  return { update: update, shown: layer.shown, element: container, clear: layer.clear };
}
