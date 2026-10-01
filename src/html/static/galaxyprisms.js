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
// Level of detail: zoomed out, one prism stands for a block ("mega
// sector") of m x m x ~m whole sectors, m a power of three, picked from
// the screen scale so a block is still BLOCK_MIN_PX wide at the focus
// (blockSizeForScale). The block grid at any one size is fixed in space,
// so panning never shifts the blocks, only adds and drops them at the
// edges; zoomed all the way in, one block is one sector. Drawn full size,
// the blocks make one solid, so only its surface is listed.
//
// Density comes from the galaxy's own analytic model
// (stellarObjects.galaxyDensity.relative_density), evaluated right here
// from the shape parameters the page embeds -- no server round trip,
// no sampled point cloud. Each prism is shaded by its mean density over
// a few sample points.
//
// No three.js import: this module only does the arithmetic and hands
// back plain typed arrays, so it runs under plain node for tests
// (src/tests/test_galaxyprisms.py) and in the page's Web Worker
// (galaxyblocks.js).

// Prisms thinner than this relative density are skipped entirely -- only
// when the shape doesn't carry the galaxy's own sector threshold
// (shape.sector_min_density, the density a sector needs to expect one
// star). With it, a prism is drawn exactly when it holds at least one
// sector the generator's skeleton would allow, so at one sector per prism
// the map's outline is the galaxy's layers.
export var PRISM_MIN_DENSITY = 0.02;
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

// The Galaxy Map's wedge lines: every master line (see ringMasterCount),
// in-plane from the ring its zone starts at out to galaxyRadius, as
// [{ bearingDeg, angleRad, r0, masters }]. r0 is the radius the line starts
// at (pc) and masters the master count of the zone it starts in (3 for the
// lines through the core, then 6, 12, ...), so the page can show coarser
// lines first. Bearings are counterclockwise from +X (the zero meridian,
// every ring's slot 0 leading edge), rounded to 0.01 degree.
export function wedgeLines(edgePc, galaxyRadius) {
  var lastRing = Math.max(0, Math.ceil(galaxyRadius / edgePc));
  var lines = [];
  var previous = 0;
  for (var ring = 0; ring <= lastRing; ring++) {
    var masters = ringMasterCount(ring);
    if (masters === previous) continue;
    for (var k = 0; k < masters; k++) {
      // Lines of a coarser zone already run through here.
      if (previous && (k * previous) % masters === 0) continue;
      var angle = (2 * Math.PI * k) / masters;
      lines.push({
        bearingDeg: Math.round((36000 * k) / masters) / 100,
        angleRad: angle, r0: ring * edgePc, masters: masters,
      });
    }
    previous = masters;
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

// --- Mega-blocks -------------------------------------------------------------
//
// Zoomed out, one prism stands for a block of whole sectors, m a side,
// where m is a power of three. The block grid is the sector grid scaled by
// m: block ring I covers sector rings I*m .. I*m + m - 1, and block layer S
// covers sector layers S*m - (m-1)/2 .. S*m + (m-1)/2, so its top and bottom
// sit on sector layer walls and block layer 0 still straddles the plane
// (an even m would split sector layer 0). At m = 1 a block is one sector.
//
// Block wedges follow the master wedges of the block's innermost ring
// (ringMasterCount; a member ring in a later zone has twice as many, so
// its master lines include these). blockWedgeCount aims for
// round(2 * pi * (I + 1/2)) wedges, so a block is about m edges long, and
// takes the nearest count that is a divisor or a multiple of the master
// count, as long as that is within BLOCK_WEDGE_TOLERANCE of the aim:
// - A divisor (large m, out in the disk): every wedge side is a master
//   line, so a slot wall in every member ring, and the block is exactly
//   the sectors inside its outline.
// - A multiple: each master wedge split into equal parts. The split lines
//   don't follow slot walls.
// - Neither close enough: the aim itself, not aligned.
// A sector belongs to the block its slot's center falls in, which for an
// aligned block is simply the sectors inside it, so every sector is in
// exactly one block. Blocks hold m^3 sectors give or take 20% (a third
// at m = 3, where a block holds only a few slots per ring).

// A block is at least this many CSS pixels across at the view's focus.
// Finer than that, single blocks stop being something a viewer can pick
// out or click. The page passes its own (data.blockMinPx, from
// lib/galaxymap3d.py); this is the default.
export var BLOCK_MIN_PX = 4;
// The most blocks one view draws (default; the page passes
// data.blockBudget). Past it, blocks get three times bigger.
export var BLOCK_BUDGET = 60000;
// How far off round(2 * pi * (I + 1/2)) a master-aligned wedge count may
// be before blockWedgeCount gives up on alignment. Tighter keeps blocks
// closer to m^3 sectors but aligns fewer of them: at 1.2, 87% of block
// rings are aligned at m = 27 and 96% at m = 81.
var BLOCK_WEDGE_TOLERANCE = 1.2;
// Density caches (see blocksInView) are cleared past this many entries.
var DENSITY_CACHE_MAX = 200000;

function isPowerOfThree(m) {
  while (m > 1 && m % 3 === 0) m /= 3;
  return m === 1;
}

// Sectors per block side for a view with pcPerPixel parsecs per CSS pixel
// at its focus: the smallest power of 3 with m * edge >= minPx *
// pcPerPixel, so one block always covers at least a pixel's worth of
// sectors. Capped where the whole disk is a handful of block rings
// (galaxyRadius / 3 across), since bigger blocks only lose the shape.
export function blockSizeForScale(pcPerPixel, edgePc, minPx, galaxyRadius) {
  var need = ((minPx || BLOCK_MIN_PX) * (pcPerPixel || 0)) / edgePc;
  var cap = galaxyRadius > 0 ? galaxyRadius / (3 * edgePc) : Infinity;
  var m = 1;
  while (m < need && m * 3 <= cap) {
    m *= 3;
  }
  return m;
}

// Wedges in block ring `ring` for blocks m sectors a side (see above).
export function blockWedgeCount(ring, m) {
  if (m === 1) return ringSectorCount(ring);
  var masters = ringMasterCount(ring * m);
  var target = Math.max(3, Math.round(2 * Math.PI * (ring + 0.5)));
  var best = masters * Math.max(1, Math.round(target / masters));
  // Master counts are 3 * 2^n, so their divisors are 2^j and 3 * 2^j.
  for (var d = 1; d < masters; d *= 2) {
    [d, 3 * d].forEach(function (w) {
      if (w < masters && Math.abs(Math.log(w / target)) < Math.abs(Math.log(best / target))) {
        best = w;
      }
    });
  }
  return Math.abs(Math.log(best / target)) <= Math.log(BLOCK_WEDGE_TOLERANCE) ? best : target;
}

// Sector address ranges a block covers: {ringFirst, ringLast, layerFirst,
// layerLast} (inclusive). Its slots are those whose centers fall in its
// wedge in each member ring (see blockSlotRange).
export function blockSectorRanges(ring, slab, m) {
  var half = (m - 1) / 2;
  return {
    ringFirst: ring * m, ringLast: ring * m + m - 1,
    layerFirst: slab * m - half, layerLast: slab * m + half,
  };
}

// The slots of sector ring `sectorRing` inside wedge `seg` of block ring
// `ring`: {first, last} (inclusive; last < first when none). `wedges`
// (optional) is blockWedgeCount(ring, m), when the caller has it.
export function blockSlotRange(ring, seg, m, sectorRing, wedges) {
  wedges = wedges || blockWedgeCount(ring, m);
  var n = ringSectorCount(sectorRing);
  // Slot k's center is at (k + 1/2) / n of a turn; it is in the wedge when
  // seg / wedges <= (k + 1/2) / n < (seg + 1) / wedges. Integer math, so
  // wedge sides on slot walls never round the wrong way.
  var first = Math.ceil((2 * seg * n - wedges) / (2 * wedges));
  var last = Math.ceil((2 * (seg + 1) * n - wedges) / (2 * wedges)) - 1;
  return { first: Math.max(0, first), last: Math.min(n - 1, last) };
}

// Whether sector (ring, layer) can exist: the galaxy's density bound at
// its ring centerline and layer center reaches shape.sector_min_density
// (the skeleton's rule, stellarObjects.galaxySkeleton). Without a
// threshold in the shape, everything counts.
export function sectorAllowed(ring, layer, edgePc, shape) {
  if (!(shape.sector_min_density > 0)) return true;
  return densityUpperBound((ring + 0.5) * edgePc, Math.abs(layer) * edgePc, shape) >= shape.sector_min_density;
}

// How many sectors the block holds: in each member ring, the slots in its
// wedge, times the member layers the skeleton allows there.
export function blockSectorCount(ring, seg, slab, m, edgePc, shape) {
  var ranges = blockSectorRanges(ring, slab, m);
  var wedges = blockWedgeCount(ring, m);
  var total = 0;
  for (var i = ranges.ringFirst; i <= ranges.ringLast; i++) {
    var slots = blockSlotRange(ring, seg, m, i, wedges);
    var perLayer = Math.max(0, slots.last - slots.first + 1);
    if (!perLayer) continue;
    for (var j = ranges.layerFirst; j <= ranges.layerLast; j++) {
      if (!shape || sectorAllowed(i, j, edgePc, shape)) total += perLayer;
    }
  }
  return total;
}

// Whether a block has anything to draw. The density bound ignores angle,
// so this depends only on the block's ring and layer. With the galaxy's
// sector threshold (shape.sector_min_density) it is exactly "holds at
// least one sector the skeleton allows": the densest member is the
// innermost ring's centerline on the member layer nearest the plane.
// Without one, blocks thinner than PRISM_MIN_DENSITY everywhere are
// skipped. Nothing at or past galaxyRadius exists.
export function blockExists(ring, slab, m, edgePc, shape, galaxyRadius) {
  var size = m * edgePc;
  if (ring < 0 || (galaxyRadius > 0 && ring * size >= galaxyRadius)) return false;
  var layers = blockSectorRanges(ring, slab, m);
  var nearLayer = layers.layerFirst <= 0 && layers.layerLast >= 0
    ? 0 : Math.min(Math.abs(layers.layerFirst), Math.abs(layers.layerLast));
  if (shape.sector_min_density > 0) {
    return sectorAllowed(ring * m, nearLayer, edgePc, shape);
  }
  var z0 = (slab - 0.5) * size;
  var z1 = z0 + size;
  var zMinAbs = z0 <= 0 && z1 >= 0 ? 0 : Math.min(Math.abs(z0), Math.abs(z1));
  return densityUpperBound(ring * size, zMinAbs, shape) >= PRISM_MIN_DENSITY;
}

// The blocks for a view: blockSizeForScale's size, made three times
// coarser while the surface listing overflows the budget. densityCacheFor(m)
// returns the density cache for that size. Options:
// - pcPerPixel: parsecs per CSS pixel at the focus;
// - sliceZ: a galaxy-frame z, or null for the whole solid; every block
//   layer above the one holding it is hidden;
// - viewerZ: the camera's z (see blocksInView);
// - minPx (BLOCK_MIN_PX) and budget (BLOCK_BUDGET).
// Returns {m, slice, blocks} (slice: the block layer cut at, or null).
export function blocksForView(center, viewRadius, edgePc, galaxyRadius, shape, densityCacheFor, options) {
  options = options || {};
  var budget = options.budget > 0 ? options.budget : BLOCK_BUDGET;
  var m = blockSizeForScale(options.pcPerPixel, edgePc, options.minPx, galaxyRadius);
  var maxM = blockSizeForScale(Infinity, edgePc, options.minPx, galaxyRadius);
  for (;;) {
    var slice = options.sliceZ == null ? null : Math.round(options.sliceZ / (m * edgePc));
    var blocks = blocksInView(center, viewRadius, m, edgePc, shape, galaxyRadius, densityCacheFor(m), {
      slice: slice, viewerZ: options.viewerZ, limit: m < maxM ? budget : Infinity,
    });
    if (blocks) {
      return { m: m, slice: slice, blocks: blocks };
    }
    m *= 3;
  }
}

// Block (ring, seg, slab) m sectors a side holding galaxy-frame point
// (x, y, z): the same cell every sector whose center is there belongs to.
export function blockAddressAt(x, y, z, m, edgePc) {
  var size = m * edgePc;
  var ring = Math.floor(Math.hypot(x, y) / size);
  var wedges = blockWedgeCount(ring, m);
  var theta = ((Math.atan2(y, x) % (2 * Math.PI)) + 2 * Math.PI) % (2 * Math.PI);
  return {
    ring: ring,
    seg: Math.min(wedges - 1, Math.floor((theta * wedges) / (2 * Math.PI))),
    slab: Math.round(z / size),
  };
}

// One block as blocksInView lists it (bounds, density, meanDensity), for
// a block the listing skipped -- an interior one holding filled sectors.
// Without a shape, density and meanDensity are null.
export function blockAt(ring, seg, slab, m, edgePc, shape, densityCache) {
  var size = m * edgePc;
  var dTheta = (2 * Math.PI) / blockWedgeCount(ring, m);
  var b = {
    ring: ring, seg: seg, slab: slab,
    r0: ring * size, r1: (ring + 1) * size, t0: seg * dTheta, t1: (seg + 1) * dTheta,
    z0: (slab - 0.5) * size, z1: (slab + 0.5) * size,
  };
  if (!shape) {
    b.density = null;
    b.meanDensity = null;
    return b;
  }
  var key = ring + "/" + seg + "/" + slab;
  var sampled = densityCache ? densityCache.get(key) : undefined;
  if (sampled === undefined) {
    sampled = meanDensity(b.r0, b.r1, b.t0, b.t1, b.z0, b.z1, shape);
    if (densityCache) {
      densityCache.set(key, sampled);
    }
  }
  b.density = sampled.density;
  b.meanDensity = sampled.mean;
  return b;
}

// Blocks m sectors a side overlapping the view ball. Returns [{ring, seg,
// slab, r0, r1, t0, t1, z0, z1, density, meanDensity}]: the block's
// address on the block grid and its bounds (see blockSectorRanges for the
// sectors inside), or null once more than options.limit are found.
// Options:
// - surfaceOnly (default true): list only blocks with a missing ring or
//   layer neighbour (the axis counts as filled), since the blocks make one
//   solid and nothing inside it shows. Existence depends only on ring and
//   layer, so interior pairs are skipped before any wedge is looked at.
// - slice: a block layer; every layer above it counts as empty, so the
//   cut face is drawn.
// - viewerZ: the camera's z. A missing neighbour below a block exposes
//   only its bottom face, which faces away from a camera above that face
//   (and the same for tops), so blocks exposed only that way are dropped.
// densityCache (a Map, optional) keeps per-block densities between calls
// at the same size.
// Unfilled blocks are translucent on the page, but the listing still
// skips the interior: unfilled space is drawn as a shell (its outline and
// the slice face), and the page adds the interior blocks that hold filled
// sectors itself (blockAt), which show through that shell.
export function blocksInView(center, viewRadius, m, edgePc, shape, galaxyRadius, densityCache, options) {
  if (!isPowerOfThree(m)) throw new Error("blocks are a power of 3 sectors a side");
  options = options || {};
  var surfaceOnly = options.surfaceOnly !== false;
  var viewerZ = options.viewerZ == null ? null : options.viewerZ;
  var slice = options.slice == null ? Infinity : options.slice;
  var limit = options.limit == null ? Infinity : options.limit;
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
  var slabHi = Math.min(slice, Math.round((cz + viewRadius) / size));
  if (densityCache && densityCache.size > DENSITY_CACHE_MAX) {
    densityCache.clear();
  }
  var existsMemo = new Map();
  function exists(ring, slab) {
    if (ring < 0) return true; // the axis
    if (slab > slice) return false;
    var key = ring + "/" + slab;
    var value = existsMemo.get(key);
    if (value === undefined) {
      value = blockExists(ring, slab, m, edgePc, shape, galaxyRadius);
      existsMemo.set(key, value);
    }
    return value;
  }
  var found = [];
  for (var ring = ringLo; ring <= ringHi; ring++) {
    var r0 = ring * size;
    var r1 = r0 + size;
    var rMid = r0 + size / 2;
    var nSeg = blockWedgeCount(ring, m);
    var dTheta = (2 * Math.PI) / nSeg;
    // Half a block's diagonal, generously.
    var reach = 0.5 * Math.hypot(size, rMid * dTheta, size);
    var segFirst = 0;
    var segCount = nSeg;
    if (rc > viewRadius) {
      var halfWidth = Math.asin(Math.min(1, viewRadius / rc)) + dTheta;
      if (halfWidth < Math.PI) {
        segFirst = Math.floor((thetaC - halfWidth) / dTheta);
        segCount = Math.min(nSeg, Math.floor((thetaC + halfWidth) / dTheta) - segFirst + 1);
      }
    }
    for (var slab = slabLo; slab <= slabHi; slab++) {
      var zMid = slab * size;
      if (Math.abs(zMid - cz) > viewRadius + reach || !exists(ring, slab)) {
        continue;
      }
      var z0 = zMid - size / 2;
      var z1 = z0 + size;
      if (surfaceOnly && exists(ring - 1, slab) && exists(ring + 1, slab)
          && (exists(ring, slab - 1) || (viewerZ !== null && viewerZ >= z0))
          && (exists(ring, slab + 1) || (viewerZ !== null && viewerZ <= z1))) {
        continue;
      }
      for (var k = 0; k < segCount; k++) {
        var seg = (((segFirst + k) % nSeg) + nSeg) % nSeg;
        var t0 = seg * dTheta;
        var t1 = t0 + dTheta;
        var tMid = t0 + dTheta / 2;
        var dxy = Math.hypot(rMid * Math.cos(tMid) - cx, rMid * Math.sin(tMid) - cy);
        if (Math.hypot(dxy, zMid - cz) > viewRadius + reach) {
          continue;
        }
        var key = ring + "/" + seg + "/" + slab;
        var sampled = densityCache ? densityCache.get(key) : undefined;
        if (sampled === undefined) {
          sampled = meanDensity(r0, r1, t0, t1, z0, z1, shape);
          if (densityCache) {
            densityCache.set(key, sampled);
          }
        }
        if (found.length >= limit) {
          return null;
        }
        found.push({ ring: ring, seg: seg, slab: slab, r0: r0, r1: r1, t0: t0, t1: t1, z0: z0, z1: z1, density: sampled.density, meanDensity: sampled.mean });
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

// TODO(galaxy-map #19): full-size blocks share faces with their
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
