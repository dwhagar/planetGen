// html/static/sectorframe.js
//
// The frame of the sector the Galaxy Map has opened in place (MAP.327): a wireframe round its cell that stays where the
// cell is (it never follows the pointer), and at the centre of each wall that faces another sector a labelled icon whose
// click flies to that sector. The geometry (`wallCentres`) is pure and runs under plain node (tests/js/sectorframe.test.mjs);
// `createSectorFrame` adds the three.js lines and the DOM icons.
//
// A cell is bounds {r0, r1, c0, c1, t0, t1}: radius in parsecs, colatitude and longitude in radians (a block's, as
// galaxystages.cellEdgePoints reads them). Its walls are the inner and outer radius, the two colatitudes and the two
// longitudes. A wall on the galaxy's absolute boundary (no radius under zero, no colatitude past a pole, no longitude
// wall when the cell goes the whole way round) has no neighbour and no icon.

const TWO_PI = 2 * Math.PI;
const EPSILON = 1e-9;

// How far past a wall, in parsecs, the neighbour is looked for: a quarter of a sector's edge.
export const BEYOND_EDGE_SHARE = 0.25;

function at(r, c, t) {
  const plane = r * Math.sin(c);
  return [plane * Math.cos(t), plane * Math.sin(t), r * Math.cos(c)];
}

// The unit vectors of the spherical frame at (c, t): radial, colatitude (southward) and longitude (eastward).
function frameAt(c, t) {
  return {
    r: [Math.sin(c) * Math.cos(t), Math.sin(c) * Math.sin(t), Math.cos(c)],
    c: [Math.cos(c) * Math.cos(t), Math.cos(c) * Math.sin(t), -Math.sin(c)],
    t: [-Math.sin(t), Math.cos(t), 0],
  };
}

function along(point, unit, length) {
  return [point[0] + unit[0] * length, point[1] + unit[1] * length, point[2] + unit[2] * length];
}

// The walls of a cell that may face a neighbour: `[{dir, point, beyond}]`, `point` the centre of the wall and `beyond` a
// point `beyondPc` parsecs past it, outward from the cell, where the neighbour is found. `dir` is `out`, `in` (from the
// galaxy's centre), `north` and `south` (colatitude), `west` and `east` (longitude, the lesser and greater).
export function wallCentres(b, beyondPc) {
  const reach = beyondPc === undefined ? 1 : beyondPc;
  const rMid = (b.r0 + b.r1) / 2;
  const cMid = (b.c0 + b.c1) / 2;
  const fullTurn = b.t1 - b.t0 >= TWO_PI - EPSILON;
  const tMid = (b.t0 + b.t1) / 2;
  const walls = [];
  const wall = function (dir, point, unit, sign) {
    walls.push({ dir: dir, point: point, beyond: along(point, unit, sign * reach) });
  };
  const middle = frameAt(cMid, tMid);
  wall("out", at(b.r1, cMid, tMid), middle.r, 1);
  if (b.r0 > EPSILON) wall("in", at(b.r0, cMid, tMid), middle.r, -1);
  if (b.c0 > EPSILON) wall("north", at(rMid, b.c0, tMid), frameAt(b.c0, tMid).c, -1);
  if (b.c1 < Math.PI - EPSILON) wall("south", at(rMid, b.c1, tMid), frameAt(b.c1, tMid).c, 1);
  if (!fullTurn) {
    wall("west", at(rMid, cMid, b.t0), frameAt(cMid, b.t0).t, -1);
    wall("east", at(rMid, cMid, b.t1), frameAt(cMid, b.t1).t, 1);
  }
  return walls;
}

// The arrow and the words a wall's icon carries before the neighbour's name.
export const WALL_LABELS = {
  out: { arrow: "↑", words: "outward" },
  in: { arrow: "↓", words: "inward" },
  north: { arrow: "▲", words: "north" },
  south: { arrow: "▼", words: "south" },
  west: { arrow: "◀", words: "west" },
  east: { arrow: "▶", words: "east" },
};

// A frame in `scene` for `camera` over `canvasEl`, its icons in `parent`. `options`: `THREE`, `scene`, `camera`, `canvasEl`,
// `parent`, `accent` (the line color), `edgePc`, `neighbour(wall) -> {shell, band, wedge}` or null (the sector beside the
// cell across `wall`, `{dir, point, beyond}` as `wallCentres` gives it; null where the galaxy ends), `label(sector) -> text`, `cellPoints(bounds) -> flat [x, y, z, ...]`
// (the wireframe's segments), `onMove(sector)` (an icon was clicked).
export function createSectorFrame(options) {
  const THREE = options.THREE;
  const group = new THREE.Group();
  group.renderOrder = 7;
  options.scene.add(group);
  const container = document.createElement("div");
  container.className = "map-wall-icons";
  options.parent.appendChild(container);
  const point = new THREE.Vector3();
  let icons = [];
  let line = null;
  let bounds = null;

  function clear() {
    if (line) {
      group.remove(line);
      line.geometry.dispose();
      line.material.dispose();
      line = null;
    }
    icons.forEach(function (icon) { icon.element.remove(); });
    icons = [];
  }

  // Frames the cell `cell` (bounds), or nothing for null.
  function set(cell) {
    clear();
    bounds = cell || null;
    if (!bounds) return;
    // MAP.240: the frame is an object at the middle of its edges, with offsets from it.
    const points = options.cellPoints(bounds);
    const origin = [0, 1, 2].map(function (axis) {
      let low = Infinity;
      let high = -Infinity;
      for (let i = axis; i < points.length; i += 3) {
        low = Math.min(low, points[i]);
        high = Math.max(high, points[i]);
      }
      return Number.isFinite(low) ? (low + high) / 2 : 0;
    });
    const offsets = new Float32Array(points.length);
    for (let i = 0; i < points.length; i++) offsets[i] = points[i] - origin[i % 3];
    const geometry = new THREE.BufferGeometry();
    geometry.setAttribute("position", new THREE.BufferAttribute(offsets, 3));
    line = new THREE.LineSegments(geometry, new THREE.LineBasicMaterial({
      color: options.accent, depthTest: false, depthWrite: false, transparent: true, opacity: 0.85,
    }));
    line.position.set(origin[0], origin[1], origin[2]);
    line.userData.relativeLayer = "sector frame";
    line.renderOrder = 7;
    line.frustumCulled = false;
    group.add(line);
    const seen = new Set();
    wallCentres(bounds, BEYOND_EDGE_SHARE * options.edgePc).forEach(function (wall) {
      let sector = null;
      try {
        sector = options.neighbour(wall);
      } catch (error) {
        sector = null;
      }
      if (!sector) return;
      const name = options.label(sector);
      if (seen.has(name)) return;
      seen.add(name);
      const labels = WALL_LABELS[wall.dir];
      const element = document.createElement("button");
      element.type = "button";
      element.className = "map-wall-icon";
      element.dataset.dir = wall.dir;
      element.title = "Go to " + name + ", " + labels.words;
      element.setAttribute("aria-label", "Go to " + name + ", " + labels.words);
      const arrow = document.createElement("span");
      arrow.className = "map-wall-arrow";
      arrow.textContent = labels.arrow;
      const text = document.createElement("span");
      text.className = "map-wall-name";
      text.textContent = name;
      element.appendChild(arrow);
      element.appendChild(text);
      element.addEventListener("click", function (event) {
        event.stopPropagation();
        options.onMove(sector);
      });
      // A press on an icon is not a press on the map under it (a drag, a click on a star).
      ["pointerdown", "pointerup", "dblclick"].forEach(function (type) {
        element.addEventListener(type, function (event) { event.stopPropagation(); });
      });
      element.hidden = true;
      container.appendChild(element);
      icons.push({ element: element, wall: wall, sector: sector, name: name });
    });
  }

  // Puts each icon at its wall's centre on the screen; one whose wall is behind the camera or off the map is hidden.
  function update() {
    if (!icons.length) return;
    const width = options.canvasEl.clientWidth || 1;
    const height = options.canvasEl.clientHeight || 1;
    options.camera.updateMatrixWorld();
    icons.forEach(function (icon) {
      point.set(icon.wall.point[0], icon.wall.point[1], icon.wall.point[2]).project(options.camera);
      const x = (point.x + 1) / 2 * width;
      const y = (1 - point.y) / 2 * height;
      const visible = point.z > -1 && point.z < 1 && x >= 0 && x <= width && y >= 0 && y <= height;
      icon.element.hidden = !visible;
      if (visible) {
        icon.element.style.left = x.toFixed(1) + "px";
        icon.element.style.top = y.toFixed(1) + "px";
      }
    });
  }

  return {
    set: set,
    update: update,
    group: group,
    // For the browser tests: the frame's lines, and each icon's wall, sector, text and whether it shows.
    state: function () {
      return {
        framed: !!line,
        icons: icons.map(function (icon) {
          return { dir: icon.wall.dir, name: icon.name, sector: icon.sector, shown: !icon.element.hidden, wall: icon.wall.point };
        }),
      };
    },
    positions: function () { return line ? Array.from(line.geometry.getAttribute("position").array) : []; },
  };
}
