// html/static/galaxyprisms.js
//
// The 3D Galaxy Map's sector grid and predicted-density shading
// (galaxymap3d.js), drawn as cylindrical segment prisms. Space is cut the
// way stellarObjects/galaxyGeometry.py cuts it into sectors: rings one
// sector edge wide around the galactic axis, layers one edge tall (layer 0
// straddles z = 0) and each ring split into wedges about one edge long.
// One cell is one prism: an annular wedge with curved inner and outer
// walls, flat top and bottom, and two flat radial sides.
//
// Level of detail: zoomed out, one prism stands for a group ("mega
// sector") of m x m x ~m whole sectors, m a power of three, picked per
// view so a group is still big enough on screen to see and click
// (MIN_PRISM_PX) and the prisms in view stay within PRISM_BUDGET. The
// grid at any one size is fixed in space, so panning never shifts the
// prisms, only adds and drops them at the edges; zoomed all the way in,
// one prism is one sector.
//
// Density comes from the galaxy's own analytic model
// (stellarObjects.galaxyDensity.relative_density), evaluated right here
// from the shape parameters the page embeds -- no server round trip,
// no sampled point cloud. Each prism is shaded by its mean density over
// a few sample points.
//
// No three.js import: this module only does the arithmetic and hands
// back plain typed arrays, so it runs under plain node for tests
// (src/tests/test_galaxyprisms.py).

// TODO(galaxy-map #13): the plan replacing the level of detail above
// (galaxy-megablocks/report.md, approved 2026-09-30):
//   - Block size from the screen, not a volume guess: m is the smallest
//     power of 3 with m * edge >= BLOCK_MIN_PX (4) * pcPerPixel at the
//     focus. One block then always covers at least a pixel's worth of
//     sectors.
//   - Blocks are drawn full size, making one continuous solid (#14), and
//     only exposed blocks are listed. Whether a block exists depends only on
//     its ring and layer (the density bound ignores angle), so interior
//     blocks are skipped before any wedge is looked at.
//   - A slice layer hides blocks above the cut (#14).
//   - Everything here must stay free of three.js, so it can run in a Web
//     Worker (#17).

// How many prisms one view may draw, roughly.
export var PRISM_BUDGET = 14000;
// Prisms thinner than this relative density are skipped entirely -- only
// when the shape doesn't carry the galaxy's own sector threshold
// (shape.sector_min_density, the density a sector needs to expect one
// star). With it, a prism is drawn exactly when it holds at least one
// sector the generator's skeleton would allow, so at one sector per prism
// the map's outline is the galaxy's layers.
export var PRISM_MIN_DENSITY = 0.02;
// The disk counts as this many scale heights thick when estimating how
// many prisms a view holds.
var DISK_EXTENT_SCALE_HEIGHTS = 3;
// Curved walls get one segment per this many radians of arc.
var ARC_SEGMENT_RAD = 0.12;
var MAX_ARC_SEGMENTS = 24;
// Sample points per prism: radial x azimuth x vertical.
var SAMPLES_R = 2;
var SAMPLES_THETA = 2;
var SAMPLES_Z = 3;

function sechSquared(x) {
  var ax = Math.abs(x);
  if (ax > 700) {
    return 0;
  }
  var e = Math.exp(-2 * ax);
  return (4 * e) / ((1 + e) * (1 + e));
}

// Port of galaxyDensity.relative_density, same field names as the
// GalaxyShape namedtuple.
export function relativeDensity(x, y, z, shape) {
  return densityParts(x, y, z, shape)[0];
}

// [density, the same with the arm factor held at 1]: the second is the
// density's azimuthal mean around the ring (bulge + disk, no arms), since
// the arm factor averages to 1 around any circle.
function densityParts(x, y, z, shape) {
  var rCyl = Math.hypot(x, y);
  var r3d = Math.sqrt(x * x + y * y + z * z);
  var bulge = shape.bulge_amplitude * Math.exp(-r3d / shape.bulge_scale_radius_pc);
  var disk = Math.exp(-rCyl / shape.disk_scale_length_pc) * sechSquared(z / shape.disk_scale_height_pc);
  var armFactor = 1;
  if (rCyl > 1e-9) {
    var theta = Math.atan2(y, x);
    var thetaArm = shape.spiral_reference_angle_rad
      + Math.log(rCyl / shape.spiral_reference_radius_pc) / Math.tan(shape.pitch_angle_rad);
    armFactor = 1 + shape.arm_amplitude * Math.cos(shape.arm_count * (theta - thetaArm));
  }
  return [shape.k_norm * (bulge + disk * armFactor), shape.k_norm * (bulge + disk)];
}

// The most any point in the prism could have -- cheap enough to run on
// every candidate, so empty space is skipped before it's sampled.
function densityUpperBound(r0, zMinAbs, shape) {
  var r3dMin = Math.hypot(r0, zMinAbs);
  var bulge = shape.bulge_amplitude * Math.exp(-r3dMin / shape.bulge_scale_radius_pc);
  var disk = Math.exp(-r0 / shape.disk_scale_length_pc)
    * sechSquared(zMinAbs / shape.disk_scale_height_pc)
    * (1 + Math.abs(shape.arm_amplitude));
  return shape.k_norm * (bulge + disk);
}

// --- The sector grid (mirrors stellarObjects/galaxyGeometry.py) -----------
//
// Ring i: cylindrical radius [i*e, (i+1)*e). Layer j: z in [(j-1/2)*e,
// (j+1/2)*e), so layer 0 straddles the plane. Slot k: one of
// ringSectorCount(i) equal wedges counterclockwise from +X. Every point in
// space falls in exactly one cell.

// Master wedges ring `ring` sits in: 3 at the center, doubling (6, 12,
// ... 1,536) at the first ring where each doubled wedge would hold at
// least 8 slots. Every master line runs from where it starts out to the
// edge, and is a slot boundary in every ring outward. Mirrors
// galaxyGeometry.ring_master_count.
export function ringMasterCount(ring) {
  var c = 2 * Math.PI * (ring + 0.5);
  var master = 3;
  if (!isFinite(c)) return NaN;
  while (c >= 2 * master * 8) master *= 2;
  return master;
}

// Slots in sector ring `ring`: the multiple of its master count nearest
// 2 * pi * (i + 1/2) -- 3, 9, 15, 21, 27, 36, ... -- so each slot's arc is
// within about 6% of one edge. The same on every layer, so the columns
// line up through the whole stack. Mirrors galaxyGeometry.ring_sector_count.
export function ringSectorCount(ring) {
  var master = ringMasterCount(ring);
  return Math.max(master, master * Math.round((2 * Math.PI * (ring + 0.5)) / master));
}

// The same rule for a group grid's rings (kept under its old name).
export var azimuthSegments = ringSectorCount;

// The Galaxy Map's wedge lines: in-plane lines from the core out to the
// edge, one per ring 0 slot boundary, as [{ bearingDeg, angleRad, r0 }]
// (r0: the radius the line starts at, pc). Bearings are counterclockwise
// from +X (the zero meridian, ring slot 0's leading edge).
// TODO(galaxy-map #13): with #12's master wedges, return every master line
// (3 at the core, doubling outward) with the radius its zone starts at.
export function wedgeLines() {
  var n = ringSectorCount(0);
  var lines = [];
  for (var k = 0; k < n; k++) {
    var angle = (2 * Math.PI * k) / n;
    lines.push({ bearingDeg: Math.round((360 * k) / n) % 360, angleRad: angle, r0: 0 });
  }
  return lines;
}

// The (ring, layer, slot) cell holding galaxy-frame point (x, y, z).
export function sectorAddressAt(x, y, z, edgePc) {
  var ring = Math.floor(Math.hypot(x, y) / edgePc);
  var layer = Math.floor(z / edgePc + 0.5);
  var n = ringSectorCount(ring);
  var theta = ((Math.atan2(y, x) % (2 * Math.PI)) + 2 * Math.PI) % (2 * Math.PI);
  var slot = Math.min(n - 1, Math.floor((theta * n) / (2 * Math.PI)));
  return { ring: ring, layer: layer, slot: slot };
}

// One sector's cell as bounds: {r0, r1, t0, t1, z0, z1}.
export function sectorCellBounds(ring, layer, slot, edgePc) {
  var step = (2 * Math.PI) / ringSectorCount(ring);
  return {
    r0: ring * edgePc, r1: (ring + 1) * edgePc,
    t0: slot * step, t1: (slot + 1) * step,
    z0: (layer - 0.5) * edgePc, z1: (layer + 0.5) * edgePc,
  };
}

// The 8 corners of a cell given as bounds, galaxy frame. Index is
// 4*r_bit + 2*z_bit + theta_bit, as galaxyGeometry.sector_cell_vertices_pc.
export function cellVertices(b) {
  var out = [];
  [b.r0, b.r1].forEach(function (r) {
    [b.z0, b.z1].forEach(function (z) {
      [b.t0, b.t1].forEach(function (t) {
        out.push([r * Math.cos(t), r * Math.sin(t), z]);
      });
    });
  });
  return out;
}

// A cell's center (ring centerline, middle angle, mid height) in the three
// coordinate systems the map shows: Cartesian (x, y, z), cylindrical
// (R, theta, z) and spherical (r, theta, polar angle from galactic north).
export function cellCoordinates(b) {
  var r = (b.r0 + b.r1) / 2;
  var t = (b.t0 + b.t1) / 2;
  var z = (b.z0 + b.z1) / 2;
  var x = r * Math.cos(t);
  var y = r * Math.sin(t);
  var rho = Math.sqrt(x * x + y * y + z * z);
  return {
    cartesian: [x, y, z],
    cylindrical: [r, t, z],
    spherical: [rho, t, rho > 0 ? Math.acos(z / rho) : 0],
  };
}

// TODO(galaxy-map #13): with #12's master-aligned rings, a block's wedges
// become ringMasterCount(first member ring) / 2**p for the largest p that
// keeps the wedge arc at least m edges. Every wedge side is then a real
// slot wall in every member ring, so a block is an exact set of whole
// sectors and groupSectorCount can count slots per master wedge instead
// of binning by center angle.
// Edge cases:
// - A block that spans a zone boundary must use the innermost ring's
//   (smaller) master count; the outer zone's masters are multiples of it.
// - Ring 0's block has 3 wedges at most; blocks at the axis have no inner
//   wall.
// - A block that crosses the galaxy's outline counts only the sectors the
//   skeleton allows (galaxy_layer bounds), not all of them.
// - m larger than the galaxy: cap m at the size where the whole disk is a
//   handful of blocks.

// --- Groups of sectors ("mega sectors") -------------------------------------
//
// Zoomed out, one prism stands for a block of whole sectors: m rings by m
// layers by the slots whose centers fall inside one wedge of the group
// grid, where m (sectors per side) is a power of three. Odd sizes keep a
// group's top and bottom on sector layer boundaries while group layer 0
// still straddles the plane (a power of two would split layer 0 in half),
// and group rings always start on a sector ring. The group grid uses the
// sector grid's own rules, scaled by m, so a group is a cylindrical box
// about m sectors on every side, and at m = 1 a group is one sector.

// Sector address ranges a group covers: {ringFirst, ringLast, layerFirst,
// layerLast} (inclusive). Its slots are those whose centers fall in
// [t0, t1) of each member ring (see groupSectorCount).
export function groupSectorRanges(ring, slab, m) {
  var half = (m - 1) / 2;
  return {
    ringFirst: ring * m, ringLast: ring * m + m - 1,
    layerFirst: slab * m - half, layerLast: slab * m + half,
  };
}

// How many sectors a group holds: in each member ring, the slots whose
// center angle lies in [t0, t1), times its m layers.
export function groupSectorCount(prism, m) {
  var ranges = groupSectorRanges(prism.ring, prism.slab, m);
  var total = 0;
  for (var i = ranges.ringFirst; i <= ranges.ringLast; i++) {
    var n = ringSectorCount(i);
    var step = (2 * Math.PI) / n;
    var first = Math.ceil(prism.t0 / step - 0.5 - 1e-9);
    var last = Math.ceil(prism.t1 / step - 0.5 - 1e-9) - 1;
    total += Math.max(0, last - first + 1);
  }
  return total * m;
}

// A group's prism is at least this many screen pixels across at the view's
// focus: finer than that, single groups stop being something a viewer can
// pick out or click, and the map turns to noise.
// TODO(galaxy-map #13): becomes BLOCK_MIN_PX = 4, and sectorsPerPrism's
// volume estimate goes away. Today it assumes the whole view ball is full,
// so it overestimates the thin disk and leaves 70-290 px groups.
export var MIN_PRISM_PX = 10;

function nextPowerOfThree(x) {
  var m = 1;
  while (m < x) {
    m *= 3;
  }
  return m;
}

// Sectors per group side (a power of three) for a view of radius
// viewRadius around center: large enough that one group is at least
// MIN_PRISM_PX wide on screen (pcPerPixel is the scale at the focus; 0 or
// missing skips that check), and large enough that the estimated prism
// count fits the budget. The estimate is the part of the view ball inside
// the galaxy's disk (and radius), divided by one prism's volume.
export function sectorsPerPrism(center, viewRadius, edgePc, galaxyRadius, shape, pcPerPixel) {
  var diskHalf = DISK_EXTENT_SCALE_HEIGHTS * shape.disk_scale_height_pc;
  var zLo = Math.max(center[2] - viewRadius, -diskHalf);
  var zHi = Math.min(center[2] + viewRadius, diskHalf);
  var thickness = Math.max(zHi - zLo, 0);
  // The core's bulge is round, not flat.
  thickness = Math.max(thickness, Math.min(2 * viewRadius, 2 * DISK_EXTENT_SCALE_HEIGHTS * shape.bulge_scale_radius_pc));
  var across = Math.min(viewRadius, galaxyRadius || viewRadius);
  var volume = Math.min((4 / 3) * Math.PI * Math.pow(viewRadius, 3), Math.PI * across * across * thickness);
  var budgetEdge = Math.cbrt(Math.max(volume, 0) / PRISM_BUDGET);
  var perceptualEdge = (pcPerPixel || 0) * MIN_PRISM_PX;
  return nextPowerOfThree(Math.max(1, Math.max(budgetEdge, perceptualEdge) / edgePc));
}

// The prisms for a view: sectorsPerPrism's size to start with, made
// coarser while it overflows the budget, and finer (a third the size)
// while that still fits the budget and stays MIN_PRISM_PX wide on screen
// -- the estimate is rough, since the thin outer disk holds far fewer
// prisms than its area suggests. densityCacheFor(m) returns the density
// cache for that size. Returns {sectorsPerPrism, prisms}.
export function prismsForView(center, viewRadius, edgePc, galaxyRadius, shape, densityCacheFor, pcPerPixel) {
  var m = sectorsPerPrism(center, viewRadius, edgePc, galaxyRadius, shape, pcPerPixel);
  var prisms = prismsInView(center, viewRadius, m, edgePc, shape, galaxyRadius, densityCacheFor(m));
  for (var up = 0; up < 3 && prisms.length > PRISM_BUDGET; up++) {
    m *= 3;
    prisms = prismsInView(center, viewRadius, m, edgePc, shape, galaxyRadius, densityCacheFor(m));
  }
  var minEdge = (pcPerPixel || 0) * MIN_PRISM_PX;
  for (var down = 0; down < 2 && m > 1 && (m / 3) * edgePc >= minEdge && prisms.length < PRISM_BUDGET / 27; down++) {
    var finer = prismsInView(center, viewRadius, m / 3, edgePc, shape, galaxyRadius, densityCacheFor(m / 3));
    if (finer.length > PRISM_BUDGET) {
      break;
    }
    m /= 3;
    prisms = finer;
  }
  return { sectorsPerPrism: m, prisms: prisms };
}

// Every group of m sectors a side overlapping the view ball, dense enough
// to draw. Returns [{ring, seg, slab, r0, r1, t0, t1, z0, z1, density}]:
// the group's address on the group grid and its bounds (see
// groupSectorRanges for the sectors inside). densityCache (a Map,
// optional) keeps per-prism densities between calls at the same size.
// TODO(galaxy-map #13/#14): split into
//   blockExists(ring, slab, m)
//     the density-bound test, memoized per (ring, slab);
//   surfaceBlocksInView(center, viewRadius, m, slice)
//     walks (ring, slab) pairs, skips any whose four ring/layer neighbours
//     all exist (the axis counts as filled), and only then lists that
//     pair's wedges inside the view.
// A slab above `slice` counts as empty, so the cut face is drawn.
// TODO(galaxy-map #15): the interior skip is only valid for opaque blocks.
// A translucent (unfilled) block lets its neighbours show through. So:
// - cull interior blocks only when every neighbour is filled;
// - list unfilled interior blocks while they are within the view and the
//   budget;
// - otherwise, at large m, draw the unfilled volume as a thinner shell
//   (the outline, plus the slice face).
// Each listed block gets `filled` (the count of its generated sectors,
// from the tiles' placed lists) and `total` (groupSectorCount) alongside
// `density`, so the page can set its opacity from filled / total.
// Edge cases:
// - the view ball reaching past the galaxy's edge;
// - a view centred on the axis (every wedge in view);
// - the slice exactly on a block boundary.
// The density cache is capped (it grew without bound in the mock).
export function prismsInView(center, viewRadius, m, edgePc, shape, galaxyRadius, densityCache) {
  var size = m * edgePc;
  var cx = center[0];
  var cy = center[1];
  var cz = center[2];
  var rc = Math.hypot(cx, cy);
  var thetaC = Math.atan2(cy, cx);
  var rLimit = galaxyRadius > 0 ? galaxyRadius : Infinity;
  var ringLo = Math.max(0, Math.floor((rc - viewRadius) / size));
  var ringHi = Math.floor(Math.min(rc + viewRadius, rLimit) / size);
  var slabLo = Math.round((cz - viewRadius) / size);
  var slabHi = Math.round((cz + viewRadius) / size);
  var sectorMin = shape.sector_min_density > 0 ? shape.sector_min_density : null;
  var found = [];
  for (var ring = ringLo; ring <= ringHi; ring++) {
    var r0 = ring * size;
    var r1 = r0 + size;
    var nSeg = ringSectorCount(ring);
    var dTheta = (2 * Math.PI) / nSeg;
    var segFirst = 0;
    var segCount = nSeg;
    if (rc > viewRadius) {
      var halfWidth = Math.asin(Math.min(1, viewRadius / rc)) + dTheta;
      if (halfWidth < Math.PI) {
        segFirst = Math.floor((thetaC - halfWidth) / dTheta);
        segCount = Math.min(nSeg, Math.floor((thetaC + halfWidth) / dTheta) - segFirst + 1);
      }
    }
    var rMid = r0 + size / 2;
    for (var k = 0; k < segCount; k++) {
      var seg = (((segFirst + k) % nSeg) + nSeg) % nSeg;
      var t0 = seg * dTheta;
      var t1 = t0 + dTheta;
      var tMid = t0 + dTheta / 2;
      var mx = rMid * Math.cos(tMid);
      var my = rMid * Math.sin(tMid);
      // Half the prism's diagonal, generously.
      var reach = 0.5 * Math.hypot(size, rMid * dTheta, size);
      var dxy = Math.hypot(mx - cx, my - cy);
      if (dxy > viewRadius + reach) {
        continue;
      }
      for (var slab = slabLo; slab <= slabHi; slab++) {
        var z0 = (slab - 0.5) * size;
        var z1 = z0 + size;
        var zMid = slab * size;
        if (Math.hypot(dxy, zMid - cz) > viewRadius + reach) {
          continue;
        }
        if (sectorMin !== null) {
          // The densest member sector center could be: the innermost
          // member ring's centerline, the member layer nearest the plane.
          var layers = groupSectorRanges(ring, slab, m);
          var nearLayer = layers.layerFirst <= 0 && layers.layerLast >= 0
            ? 0 : Math.min(Math.abs(layers.layerFirst), Math.abs(layers.layerLast));
          if (densityUpperBound(r0 + edgePc / 2, nearLayer * edgePc, shape) < sectorMin) {
            continue;
          }
        } else {
          var zMinAbs = z0 <= 0 && z1 >= 0 ? 0 : Math.min(Math.abs(z0), Math.abs(z1));
          if (densityUpperBound(r0, zMinAbs, shape) < PRISM_MIN_DENSITY) {
            continue;
          }
        }
        var key = ring + "/" + seg + "/" + slab;
        var sampled = densityCache ? densityCache.get(key) : undefined;
        if (sampled === undefined) {
          sampled = meanDensity(r0, r1, t0, t1, z0, z1, shape);
          if (densityCache) {
            densityCache.set(key, sampled);
          }
        }
        var density = sampled.density;
        if (sectorMin === null && density < PRISM_MIN_DENSITY) {
          continue;
        }
        found.push({ ring: ring, seg: seg, slab: slab, r0: r0, r1: r1, t0: t0, t1: t1, z0: z0, z1: z1, density: density, meanDensity: sampled.mean });
      }
    }
  }
  return found;
}

// The prism's mean density and its azimuthal mean (bulge + disk, arm
// factor 1), both over a few sample points weighted by the volume each
// stands for (r dr): the outer samples of a wide prism cover more space
// than the inner ones. The page shades by both: density / mean is the arm
// factor, which picks out the spiral.
function meanDensity(r0, r1, t0, t1, z0, z1, shape) {
  var total = 0;
  var totalMean = 0;
  var weights = 0;
  for (var i = 0; i < SAMPLES_R; i++) {
    var r = r0 + ((i + 0.5) / SAMPLES_R) * (r1 - r0);
    for (var j = 0; j < SAMPLES_THETA; j++) {
      var t = t0 + ((j + 0.5) / SAMPLES_THETA) * (t1 - t0);
      var x = r * Math.cos(t);
      var y = r * Math.sin(t);
      for (var k = 0; k < SAMPLES_Z; k++) {
        var parts = densityParts(x, y, z0 + ((k + 0.5) / SAMPLES_Z) * (z1 - z0), shape);
        total += r * parts[0];
        totalMean += r * parts[1];
        weights += r;
      }
    }
  }
  return { density: total / weights, mean: totalMean / weights };
}

// --- Geometry --------------------------------------------------------------
//
// One mesh for every prism in view: positions, per-vertex normals and a
// per-vertex prism index (the caller turns that into a color). Faces are
// wound counter-clockwise seen from outside, so front-face culling shows
// each prism's outside only.

// TODO(galaxy-map #17): build this in a Web Worker, returning the typed
// arrays as transferables. Full-size blocks share faces with their
// neighbours, so skip any face whose neighbour exists: the surface
// listing only removes whole blocks, not hidden faces. If the vertex
// count still hurts on phones, move to one InstancedMesh per wedge-arc
// count.
export function buildPrismGeometry(prisms) {
  var vertexCount = 0;
  var indexCount = 0;
  var arcsOf = new Uint16Array(prisms.length);
  prisms.forEach(function (p, n) {
    var arcs = Math.max(1, Math.min(MAX_ARC_SEGMENTS, Math.ceil((p.t1 - p.t0) / ARC_SEGMENT_RAD)));
    arcsOf[n] = arcs;
    var curved = p.r0 > 0 ? 4 : 3; // outer, [inner,] top, bottom
    vertexCount += curved * 2 * (arcs + 1) + 8;
    indexCount += curved * 6 * arcs + 12;
  });
  var positions = new Float32Array(vertexCount * 3);
  var normals = new Float32Array(vertexCount * 3);
  var owners = new Uint32Array(vertexCount);
  var uvs = new Float32Array(vertexCount * 2);
  var indices = vertexCount > 65535 ? new Uint32Array(indexCount) : new Uint16Array(indexCount);
  var nv = 0;
  var ni = 0;

  // (u, v) is the vertex's place across its own face, 0..1 each way --
  // what the shader finds the face's edges by.
  function vertex(x, y, z, nx, ny, nz, owner, u, v) {
    uvs[2 * nv] = u;
    uvs[2 * nv + 1] = v;
    var o = 3 * nv;
    positions[o] = x;
    positions[o + 1] = y;
    positions[o + 2] = z;
    normals[o] = nx;
    normals[o + 1] = ny;
    normals[o + 2] = nz;
    owners[nv] = owner;
    nv++;
  }

  // Two rows of (arcs + 1) vertices already written from `base`; joins
  // them into quads, each wound to face its own normal.
  function stitch(base, columns) {
    for (var i = 0; i < columns; i++) {
      var a = base + i;
      var b = a + 1;
      var d = base + columns + 1 + i;
      var c = d + 1;
      if (facesNormal(a, b, c, d)) {
        indices[ni++] = a; indices[ni++] = b; indices[ni++] = c;
        indices[ni++] = a; indices[ni++] = c; indices[ni++] = d;
      } else {
        indices[ni++] = a; indices[ni++] = c; indices[ni++] = b;
        indices[ni++] = a; indices[ni++] = d; indices[ni++] = c;
      }
    }
  }

  // Whether quad a-b-c-d, taken in that order, winds counter-clockwise
  // around a's normal. (c - a) x (d - b) is the quad's normal even when two
  // of its corners coincide, as they do where a wedge meets the axis.
  function facesNormal(a, b, c, d) {
    var ux = positions[3 * c] - positions[3 * a];
    var uy = positions[3 * c + 1] - positions[3 * a + 1];
    var uz = positions[3 * c + 2] - positions[3 * a + 2];
    var vx = positions[3 * d] - positions[3 * b];
    var vy = positions[3 * d + 1] - positions[3 * b + 1];
    var vz = positions[3 * d + 2] - positions[3 * b + 2];
    var nx = uy * vz - uz * vy;
    var ny = uz * vx - ux * vz;
    var nz = ux * vy - uy * vx;
    return nx * normals[3 * a] + ny * normals[3 * a + 1] + nz * normals[3 * a + 2] >= 0;
  }

  var cosT = new Float64Array(MAX_ARC_SEGMENTS + 1);
  var sinT = new Float64Array(MAX_ARC_SEGMENTS + 1);

  for (var n = 0; n < prisms.length; n++) {
    var p = prisms[n];
    var arcs = arcsOf[n];
    for (var k = 0; k <= arcs; k++) {
      var t = p.t0 + (k / arcs) * (p.t1 - p.t0);
      cosT[k] = Math.cos(t);
      sinT[k] = Math.sin(t);
    }
    var base;
    var row;
    // Outer wall, then inner wall (none for the ring at the very center).
    base = nv;
    for (row = 0; row < 2; row++) {
      for (k = 0; k <= arcs; k++) {
        vertex(p.r1 * cosT[k], p.r1 * sinT[k], row ? p.z1 : p.z0, cosT[k], sinT[k], 0, n, k / arcs, row);
      }
    }
    stitch(base, arcs);
    if (p.r0 > 0) {
      base = nv;
      for (row = 0; row < 2; row++) {
        for (k = 0; k <= arcs; k++) {
          vertex(p.r0 * cosT[k], p.r0 * sinT[k], row ? p.z1 : p.z0, -cosT[k], -sinT[k], 0, n, k / arcs, row);
        }
      }
      stitch(base, arcs);
    }
    // Top and bottom.
    for (var cap = 0; cap < 2; cap++) {
      var z = cap ? p.z1 : p.z0;
      var nz = cap ? 1 : -1;
      base = nv;
      for (row = 0; row < 2; row++) {
        var r = row ? p.r1 : p.r0;
        for (k = 0; k <= arcs; k++) {
          vertex(r * cosT[k], r * sinT[k], z, 0, 0, nz, n, k / arcs, row);
        }
      }
      stitch(base, arcs);
    }
    // The two radial sides.
    for (var side = 0; side < 2; side++) {
      var c = side ? cosT[arcs] : cosT[0];
      var s = side ? sinT[arcs] : sinT[0];
      var sign = side ? 1 : -1;
      base = nv;
      vertex(p.r0 * c, p.r0 * s, p.z0, -s * sign, c * sign, 0, n, 0, 0);
      vertex(p.r1 * c, p.r1 * s, p.z0, -s * sign, c * sign, 0, n, 1, 0);
      vertex(p.r0 * c, p.r0 * s, p.z1, -s * sign, c * sign, 0, n, 0, 1);
      vertex(p.r1 * c, p.r1 * s, p.z1, -s * sign, c * sign, 0, n, 1, 1);
      stitch(base, 1);
    }
  }

  return { positions: positions, normals: normals, uvs: uvs, owners: owners, indices: indices };
}
