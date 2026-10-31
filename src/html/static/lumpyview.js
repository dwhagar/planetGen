// html/static/lumpyview.js
//
// The lump of a small body drawn on a plain 2D canvas (MAP.363): the object page's picture (bodyview.js) of a moon or asteroid that its own weight has not pulled round. The System Map
// and the 3D system view hang the lump (`lumpymesh.js`: the same look record, so the same shape) on a three.js sphere; the object page has no WebGL, so this module turns the same
// shape into a grid of points on the body, turns them as the body turns (`Spin`: its axis leaned by the axial tilt, turned `angle` about the axis) and fills the triangles nearest the
// viewer into the pixels of the picture, with a depth test, smooth normals and the longitude and latitude of the surface point for the painted texture to be read at.
//
// A pure module: no canvas, no `document`. `lumpMesh(look)` is the grid, `rasterise(...)` draws one turn of it into a callback, so the tests need no browser.

const VERSION_QUERY = new URL(import.meta.url).search;
const { shapeOf, surfacePoint } = await import(`./lumpymesh.js${VERSION_QUERY}`);

export const GRID_WIDTH = 48;
export const GRID_HEIGHT = 32;
/** The grid of points on the lump: 48 around and 32 from pole to pole, the segments of the three.js sphere of the other views (`lumpymesh.lumpyGeometry`), so the three views draw the same lump. */

/**
 * The lump of `look` as a grid: `{ count, positions, directions, normals, indices, extent }`. `positions`, `directions` (the unit direction on the sphere each point is a bump of, which is
 * the longitude and latitude the texture is read at) and `normals` (smooth, in the body's own frame) hold three numbers per point; `indices` the triangles; `extent` the farthest the
 * surface reaches from the middle.
 */
export function lumpMesh(look, width = GRID_WIDTH, height = GRID_HEIGHT) {
  const shape = shapeOf(look);
  const columns = width + 1;
  const rows = height + 1;
  const count = columns * rows;
  const directions = new Float32Array(count * 3);
  const positions = new Float32Array(count * 3);
  for (let j = 0; j < rows; j++) {
    const latitude = Math.PI / 2 - (Math.PI * j) / height;
    const flat = Math.cos(latitude);
    for (let i = 0; i < columns; i++) {
      const longitude = (2 * Math.PI * (i % width)) / width;
      const at = (j * columns + i) * 3;
      const dx = flat * Math.sin(longitude);
      const dy = Math.sin(latitude);
      const dz = flat * Math.cos(longitude);
      directions[at] = dx;
      directions[at + 1] = dy;
      directions[at + 2] = dz;
      const p = surfacePoint(shape, dx, dy, dz);
      positions[at] = p[0];
      positions[at + 1] = p[1];
      positions[at + 2] = p[2];
    }
  }
  const normals = new Float32Array(count * 3);
  const point = (i, j) => (Math.min(rows - 1, Math.max(0, j)) * columns + ((i % width) + width) % width) * 3;
  for (let j = 0; j < rows; j++) {
    for (let i = 0; i < columns; i++) {
      const east = point(i + 1, j);
      const west = point(i - 1, j);
      const north = point(i, j - 1);
      const south = point(i, j + 1);
      // The surface runs north to south and west to east: the outward normal is (toward east) x (toward south)... taken so that it points away from the middle.
      const ex = positions[east] - positions[west];
      const ey = positions[east + 1] - positions[west + 1];
      const ez = positions[east + 2] - positions[west + 2];
      const sx = positions[south] - positions[north];
      const sy = positions[south + 1] - positions[north + 1];
      const sz = positions[south + 2] - positions[north + 2];
      let nx = ey * sz - ez * sy;
      let ny = ez * sx - ex * sz;
      let nz = ex * sy - ey * sx;
      const at = (j * columns + i) * 3;
      let length = Math.hypot(nx, ny, nz);
      if (length < 1e-9) {          // a pole: the surface there faces along the axis
        nx = directions[at];
        ny = directions[at + 1];
        nz = directions[at + 2];
        length = 1;
      } else if (nx * directions[at] + ny * directions[at + 1] + nz * directions[at + 2] < 0) {
        nx = -nx;
        ny = -ny;
        nz = -nz;
      }
      normals[at] = nx / length;
      normals[at + 1] = ny / length;
      normals[at + 2] = nz / length;
    }
  }
  const indices = new Uint32Array(width * height * 6);
  let k = 0;
  for (let j = 0; j < height; j++) {
    for (let i = 0; i < width; i++) {
      const a = j * columns + i;
      const b = a + 1;
      const c = a + columns;
      const d = c + 1;
      indices[k++] = a; indices[k++] = c; indices[k++] = b;
      indices[k++] = b; indices[k++] = c; indices[k++] = d;
    }
  }
  return { count: count, positions: positions, directions: directions, normals: normals, indices: indices, extent: shape.extent };
}

/** The body-frame vector `(x, y, z)` turned `angle` about the axis (y) and leaned by `tiltDeg`, as the viewer sees it: `[x right, up, z toward the viewer]`. The inverse of `Spin.bodyPoint`'s lean. */
export function toView(x, y, z, cosAngle, sinAngle, cosTilt, sinTilt, out, at) {
  const turnedX = x * cosAngle + z * sinAngle;
  const turnedZ = -x * sinAngle + z * cosAngle;
  out[at] = turnedX * cosTilt + y * sinTilt;
  out[at + 1] = -turnedX * sinTilt + y * cosTilt;
  out[at + 2] = turnedZ;
}

/**
 * Draws the lump turned `angle` radians (leaned by `look.tiltDeg`) into a `size` by `size` picture: for every pixel of the nearest surface, calls
 * `paint(pixel, longitude, latitude, nx, ny, nz)` with the pixel's index (`y * size + x`), the longitude and latitude of the body point there in the body's own frame, and the surface
 * normal in the view (x right, y up, z toward the viewer). Returns the pixels it covered, in the buffer `covered` (1 for a pixel the lump fills, 0 elsewhere), which the caller passes in.
 *
 * @param {object} mesh `lumpMesh(look)`.
 * @param {object} view `{ angle, tiltDeg }`.
 * @param {number} size The picture's side, pixels.
 * @param {number} centreX The middle of the body, pixels.
 * @param {number} centreY The middle of the body, pixels (y down).
 * @param {number} scale Pixels per unit of the mesh.
 * @param {Uint8Array} covered A `size * size` buffer, cleared by this call.
 * @param {Float32Array} depth A `size * size` buffer, reset by this call.
 * @param {function} paint The pixel callback.
 */
export function rasterise(mesh, view, size, centreX, centreY, scale, covered, depth, paint) {
  covered.fill(0);
  depth.fill(-Infinity);
  const tilt = Number.isFinite(view.tiltDeg) ? (view.tiltDeg * Math.PI) / 180 : 0;
  const cosAngle = Math.cos(view.angle);
  const sinAngle = Math.sin(view.angle);
  const cosTilt = Math.cos(tilt);
  const sinTilt = Math.sin(tilt);
  const seen = new Float32Array(mesh.count * 3);
  const faces = new Float32Array(mesh.count * 3);
  for (let n = 0; n < mesh.count; n++) {
    const at = n * 3;
    toView(mesh.positions[at], mesh.positions[at + 1], mesh.positions[at + 2], cosAngle, sinAngle, cosTilt, sinTilt, seen, at);
    toView(mesh.normals[at], mesh.normals[at + 1], mesh.normals[at + 2], cosAngle, sinAngle, cosTilt, sinTilt, faces, at);
  }
  const { indices, directions } = mesh;
  for (let t = 0; t < indices.length; t += 3) {
    const a = indices[t] * 3;
    const b = indices[t + 1] * 3;
    const c = indices[t + 2] * 3;
    const ax = centreX + seen[a] * scale;
    const ay = centreY - seen[a + 1] * scale;
    const bx = centreX + seen[b] * scale;
    const by = centreY - seen[b + 1] * scale;
    const cx = centreX + seen[c] * scale;
    const cy = centreY - seen[c + 1] * scale;
    const area = (bx - ax) * (cy - ay) - (cx - ax) * (by - ay);
    if (Math.abs(area) < 1e-9) continue;
    // No face is culled: the lump is one closed star-shaped surface, so the depth test keeps the nearest of the two faces over a pixel.
    const minX = Math.max(0, Math.floor(Math.min(ax, bx, cx)));
    const maxX = Math.min(size - 1, Math.ceil(Math.max(ax, bx, cx)));
    const minY = Math.max(0, Math.floor(Math.min(ay, by, cy)));
    const maxY = Math.min(size - 1, Math.ceil(Math.max(ay, by, cy)));
    for (let y = minY; y <= maxY; y++) {
      for (let x = minX; x <= maxX; x++) {
        const px = x + 0.5;
        const py = y + 0.5;
        const wa = ((bx - px) * (cy - py) - (cx - px) * (by - py)) / area;
        const wb = ((cx - px) * (ay - py) - (ax - px) * (cy - py)) / area;
        const wc = 1 - wa - wb;
        if (wa < -1e-6 || wb < -1e-6 || wc < -1e-6) continue;
        const z = wa * seen[a + 2] + wb * seen[b + 2] + wc * seen[c + 2];
        const pixel = y * size + x;
        if (z <= depth[pixel]) continue;
        depth[pixel] = z;
        covered[pixel] = 1;
        const dx = wa * directions[a] + wb * directions[b] + wc * directions[c];
        const dy = wa * directions[a + 1] + wb * directions[b + 1] + wc * directions[c + 1];
        const dz = wa * directions[a + 2] + wb * directions[b + 2] + wc * directions[c + 2];
        let nx = wa * faces[a] + wb * faces[b] + wc * faces[c];
        let ny = wa * faces[a + 1] + wb * faces[b + 1] + wc * faces[c + 1];
        let nz = wa * faces[a + 2] + wb * faces[b + 2] + wc * faces[c + 2];
        const length = Math.hypot(nx, ny, nz) || 1;
        nx /= length;
        ny /= length;
        nz /= length;
        paint(pixel, Math.atan2(dx, dz), Math.asin(Math.max(-1, Math.min(1, dy))), nx, ny, nz);
      }
    }
  }
}
