// html/static/galaxyprisms.js
//
// The 3D Galaxy Map's predicted-density shading (galaxymap3d.js), drawn
// as cylindrical segment prisms instead of a cloud of spheres. Space is
// cut in galactic cylindrical coordinates (r_cyl, theta, z): rings a
// whole number of sector shells wide, each ring split into equal azimuth
// wedges about as long as the ring is wide, and layers of the same height
// stacked on the galactic plane (layer 0 straddles z = 0) -- the same grid
// as the cylindrical sector plan's sectors, scaled up by a power of two.
// One cell of that grid is one prism: an annular wedge with curved inner
// and outer walls, flat top and bottom, and two flat radial sides.
//
// Level of detail: the cell size is a power-of-two number of shells
// (edge_pc each), picked per view so the prisms in view stay within
// PRISM_BUDGET -- a view of the whole galaxy uses rings hundreds of
// parsecs wide, one zoomed in on a sector uses single shells. The grid
// at any one size is fixed in space, so panning never shifts the prisms,
// only adds and drops them at the edges.
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

// How many prisms one view may draw, roughly.
export var PRISM_BUDGET = 14000;
// Prisms thinner than this relative density are skipped entirely.
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
  var rCyl = Math.hypot(x, y);
  var r3d = Math.sqrt(x * x + y * y + z * z);
  var bulge = shape.bulge_amplitude * Math.exp(-r3d / shape.bulge_scale_radius_pc);
  var diskRadial = Math.exp(-rCyl / shape.disk_scale_length_pc);
  var fz = sechSquared(z / shape.disk_scale_height_pc);
  var armFactor = 1;
  if (rCyl > 1e-9) {
    var theta = Math.atan2(y, x);
    var thetaArm = shape.spiral_reference_angle_rad
      + Math.log(rCyl / shape.spiral_reference_radius_pc) / Math.tan(shape.pitch_angle_rad);
    armFactor = 1 + shape.arm_amplitude * Math.cos(shape.arm_count * (theta - thetaArm));
  }
  return shape.k_norm * (bulge + diskRadial * fz * armFactor);
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

// Wedges in prism ring `ring`: the cylindrical sector plan's own rule for
// sector ring i, 4 * round(2 * pi * (i + 1/2) / 4) -- always a multiple of
// four, so wedge boundaries fall on the quadrant lines, and about as long
// as the ring is wide. At one shell per prism this is exactly the sector
// grid, so one prism is one sector.
export function azimuthSegments(ring) {
  return 4 * Math.max(1, Math.round((2 * Math.PI * (ring + 0.5)) / 4));
}

// Shells per prism edge (a power of two) for a view of radius viewRadius
// around center: the smallest size whose estimated prism count fits the
// budget. The estimate is the part of the view ball inside the galaxy's
// disk (and radius), divided by one prism's volume.
export function shellsPerPrism(center, viewRadius, edgePc, galaxyRadius, shape) {
  var diskHalf = DISK_EXTENT_SCALE_HEIGHTS * shape.disk_scale_height_pc;
  var zLo = Math.max(center[2] - viewRadius, -diskHalf);
  var zHi = Math.min(center[2] + viewRadius, diskHalf);
  var thickness = Math.max(zHi - zLo, 0);
  // The core's bulge is round, not flat.
  thickness = Math.max(thickness, Math.min(2 * viewRadius, 2 * DISK_EXTENT_SCALE_HEIGHTS * shape.bulge_scale_radius_pc));
  var across = Math.min(viewRadius, galaxyRadius || viewRadius);
  var volume = Math.min((4 / 3) * Math.PI * Math.pow(viewRadius, 3), Math.PI * across * across * thickness);
  var needed = Math.cbrt(Math.max(volume, 0) / PRISM_BUDGET);
  var shells = Math.max(1, needed / edgePc);
  return Math.pow(2, Math.ceil(Math.log2(shells)));
}

// The prisms for a view: shellsPerPrism's size to start with, halved
// while the result still uses under a third of the budget (the estimate
// is rough -- a view over the whole galaxy holds far fewer prisms than
// its disk area suggests, since the outer disk is too thin to draw).
// densityCacheFor(shells) returns the density cache for that size.
// Returns {shells, prisms}.
export function prismsForView(center, viewRadius, edgePc, galaxyRadius, shape, densityCacheFor) {
  var shells = shellsPerPrism(center, viewRadius, edgePc, galaxyRadius, shape);
  var prisms = prismsInView(center, viewRadius, shells * edgePc, shape, galaxyRadius, densityCacheFor(shells));
  for (var tries = 0; tries < 2 && shells > 1 && prisms.length < PRISM_BUDGET / 3; tries++) {
    var finer = prismsInView(center, viewRadius, (shells / 2) * edgePc, shape, galaxyRadius, densityCacheFor(shells / 2));
    if (finer.length > PRISM_BUDGET) {
      break;
    }
    shells /= 2;
    prisms = finer;
  }
  return { shells: shells, prisms: prisms };
}

// Every prism of the given size overlapping the view ball, dense enough
// to draw. Returns [{ring, seg, slab, r0, r1, t0, t1, z0, z1, density}].
// densityCache (a Map, optional) keeps per-prism densities between calls
// at the same size.
export function prismsInView(center, viewRadius, size, shape, galaxyRadius, densityCache) {
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
  var found = [];
  for (var ring = ringLo; ring <= ringHi; ring++) {
    var r0 = ring * size;
    var r1 = r0 + size;
    var nSeg = azimuthSegments(ring);
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
        var zMinAbs = z0 <= 0 && z1 >= 0 ? 0 : Math.min(Math.abs(z0), Math.abs(z1));
        if (densityUpperBound(r0, zMinAbs, shape) < PRISM_MIN_DENSITY) {
          continue;
        }
        var key = ring + "/" + seg + "/" + slab;
        var density = densityCache ? densityCache.get(key) : undefined;
        if (density === undefined) {
          density = meanDensity(r0, r1, t0, t1, z0, z1, shape);
          if (densityCache) {
            densityCache.set(key, density);
          }
        }
        if (density < PRISM_MIN_DENSITY) {
          continue;
        }
        found.push({ ring: ring, seg: seg, slab: slab, r0: r0, r1: r1, t0: t0, t1: t1, z0: z0, z1: z1, density: density });
      }
    }
  }
  return found;
}

function meanDensity(r0, r1, t0, t1, z0, z1, shape) {
  var total = 0;
  for (var i = 0; i < SAMPLES_R; i++) {
    var r = r0 + ((i + 0.5) / SAMPLES_R) * (r1 - r0);
    for (var j = 0; j < SAMPLES_THETA; j++) {
      var t = t0 + ((j + 0.5) / SAMPLES_THETA) * (t1 - t0);
      var x = r * Math.cos(t);
      var y = r * Math.sin(t);
      for (var k = 0; k < SAMPLES_Z; k++) {
        total += relativeDensity(x, y, z0 + ((k + 0.5) / SAMPLES_Z) * (z1 - z0), shape);
      }
    }
  }
  return total / (SAMPLES_R * SAMPLES_THETA * SAMPLES_Z);
}

// --- Geometry --------------------------------------------------------------
//
// One mesh for every prism in view: positions, per-vertex normals and a
// per-vertex prism index (the caller turns that into a color). Faces are
// wound counter-clockwise seen from outside, so front-face culling shows
// each prism's outside only.

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
