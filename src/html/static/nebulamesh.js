// html/static/nebulamesh.js
//
// A nebula drawn from its shape (GEN.75, MAP.103) instead of a flat sprite:
// the triangle mesh `GET /galaxy/nebula/<id>/shape` serves (the API's
// `/api/nebulae/<id>/shape`), a translucent surface seen from outside and
// from inside. The maps (galaxymap3d.js, the Sector and System Maps) share
// this file. It takes `THREE` as an argument rather than importing it, so
// every page uses its own versioned copy and plain node can test the rest.

// A nebula's sprite gives way to its mesh once the nebula is this many
// pixels in radius on screen; smaller than that the mesh is only a blob.
export var NEBULA_MESH_MIN_PX = 20;

// The mesh surface's opacity, as a share of the nebula type's edge opacity:
// the surface is seen twice (both faces) and stacks with the sprites'
// glow, so it is kept thinner than the sprite's edge.
export var NEBULA_MESH_OPACITY_SHARE = 0.55;

// `shapePath` is a URL with `{id}` where the nebula's id goes.
export function shapeUrl(shapePath, id, lod) {
  return String(shapePath).replace("{id}", encodeURIComponent(id)) + "?lod=" + (lod || "low");
}

// A fetcher of nebula shapes that asks the server once per nebula and
// level of detail, and hands every later caller the same promise. A failed
// fetch is forgotten, so a later call tries again.
export function createShapeLoader(shapePath, fetchFn) {
  var requests = new Map();
  fetchFn = fetchFn || function (url) {
    return fetch(url, { headers: { Accept: "application/json" } }).then(function (response) {
      if (!response.ok) {
        throw new Error("shape " + response.status);
      }
      return response.json();
    });
  };
  return function load(id, lod) {
    var key = id + ":" + (lod || "low");
    if (!requests.has(key)) {
      requests.set(key, fetchFn(shapeUrl(shapePath, id, lod)).catch(function (error) {
        requests.delete(key);
        throw error;
      }));
    }
    return requests.get(key);
  };
}

// The flat position array for a served shape, scaled from nebula-radius
// units to `radius` (the map's own unit, e.g. parsecs).
// `yScale` (default 1) is -1 on a map whose y axis runs the other way
// from the galaxy frame's (the Sector Map's scene).
export function scaledPositions(shape, radius, yScale) {
  yScale = yScale == null ? 1 : yScale;
  var vertices = shape.vertices || [];
  var out = new Float32Array(vertices.length * 3);
  for (var i = 0; i < vertices.length; i++) {
    out[3 * i] = vertices[i][0] * radius;
    out[3 * i + 1] = vertices[i][1] * radius * yScale;
    out[3 * i + 2] = vertices[i][2] * radius;
  }
  return out;
}

export function flatFaces(shape) {
  var faces = shape.faces || [];
  var out = new Uint32Array(faces.length * 3);
  for (var i = 0; i < faces.length; i++) {
    out[3 * i] = faces[i][0];
    out[3 * i + 1] = faces[i][1];
    out[3 * i + 2] = faces[i][2];
  }
  return out;
}

// One nebula's mesh, centred on the origin and `radius` big; the caller
// positions it. `look` is a [hex color, core opacity, edge opacity] triple.
// `options`: `yScale` (see `scaledPositions`), `depthTest` (default false:
// drawn over the scene, as the Galaxy Map's clouds are; the Sector Map
// sets it so stars in front of the cloud stay in front).
export function buildNebulaMesh(THREE, shape, radius, look, options) {
  var o = options || {};
  var geometry = new THREE.BufferGeometry();
  geometry.setAttribute("position", new THREE.BufferAttribute(scaledPositions(shape, radius, o.yScale), 3));
  geometry.setIndex(new THREE.BufferAttribute(flatFaces(shape), 1));
  var material = new THREE.MeshBasicMaterial({
    color: new THREE.Color(look[0]),
    transparent: true,
    opacity: look[2] * NEBULA_MESH_OPACITY_SHARE,
    side: THREE.DoubleSide,
    depthWrite: false,
    depthTest: !!o.depthTest,
  });
  var mesh = new THREE.Mesh(geometry, material);
  mesh.userData.baseOpacity = material.opacity;
  return mesh;
}

export function disposeNebulaMesh(mesh) {
  mesh.geometry.dispose();
  mesh.material.dispose();
}
