// html/static/galaxyminimap.js
//
// MAP.75: the Galaxy Map's mini map, a second camera on the same renderer rather than a renderer of its own (docs/design/interface/
// map-ui-and-frontend-libraries.md section 5.5). It is the whole galaxy seen from outside, from the direction the main camera looks
// along, with a cross where the main view is centered and a line to where its eye is: the picture that says "you are here".
//
// One WebGL context: the main pass draws the whole canvas, then this draws a small square in a corner of it with the scissor test on,
// the depth buffer cleared inside the square, and its own scene and camera, so nothing of the main scene (its shaders read the main
// camera and the canvas size) has to know about it. Its zoom policy is MAP.61's "locked" one (mapcontrol.js): the camera stays at the
// distance that fits the galaxy and neither the wheel nor the buttons move it.

const VERSION_QUERY = new URL(import.meta.url).search;
const MC = await import(`./mapcontrol.js${VERSION_QUERY}`);

// The share of the canvas's short side the square takes, its bounds in CSS pixels, and the margin from the corner.
export const MINI_SHARE = 0.26;
export const MINI_MIN_SIDE = 96;
export const MINI_MAX_SIDE = 220;
export const MINI_MARGIN = 12;
// A canvas whose short side is under this has no room for it.
export const MINI_MIN_CANVAS = 320;
// The galaxy fills this share of the square's half-width.
const FIT_SHARE = 0.86;

// The zoom policy of the mini map: locked, so the camera stays at its fit.
export const MINI_POLICY = MC.zoomPolicy(MC.ZOOM_LOCKED);

// The square of the canvas (CSS pixels, origin at the bottom left like the renderer's) the mini map takes in the bottom right
// corner, or null when the canvas is too small for it.
export function miniRect(canvasWidth, canvasHeight) {
  const shortSide = Math.min(canvasWidth, canvasHeight);
  if (!(shortSide >= MINI_MIN_CANVAS)) return null;
  const side = Math.round(Math.max(MINI_MIN_SIDE, Math.min(MINI_MAX_SIDE, shortSide * MINI_SHARE)));
  return { x: canvasWidth - side - MINI_MARGIN, y: MINI_MARGIN, side: side };
}

// How far the mini camera stands from the galaxy's center so a sphere of `galaxyRadius` fills FIT_SHARE of the square, for a camera
// of `fovDeg` degrees (the vertical field of a square view). The locked policy returns exactly this whatever distance is asked.
export function miniFitDistance(galaxyRadius, fovDeg) {
  const half = (fovDeg * Math.PI) / 360;
  return (galaxyRadius / FIT_SHARE) / Math.sin(half);
}

// Where the mini camera stands: `distance` from the origin along the direction the main camera's orbit (theta from +x in the xy
// plane, phi from +z, radians) looks from.
export function miniEyeAt(orbit, distance) {
  const sinPhi = Math.sin(orbit.phi);
  return [distance * sinPhi * Math.cos(orbit.theta), distance * sinPhi * Math.sin(orbit.theta), distance * Math.cos(orbit.phi)];
}

// The vertices (as flat x, y, z pairs for a THREE.LineSegments) of the galaxy's outline: rings in the plane at a third, two thirds and
// the whole radius, a pair of thin rings above and below the disk's edge so its thickness reads, and the two axes of the plane.
export function outlineSegments(galaxyRadius, ringSegments) {
  const steps = ringSegments || 64;
  const out = [];
  function ring(radius, height) {
    for (let i = 0; i < steps; i++) {
      const a = (i / steps) * 2 * Math.PI;
      const b = ((i + 1) / steps) * 2 * Math.PI;
      out.push(radius * Math.cos(a), radius * Math.sin(a), height, radius * Math.cos(b), radius * Math.sin(b), height);
    }
  }
  ring(galaxyRadius / 3, 0);
  ring((2 * galaxyRadius) / 3, 0);
  ring(galaxyRadius, 0);
  ring(galaxyRadius, galaxyRadius * 0.05);
  ring(galaxyRadius, -galaxyRadius * 0.05);
  out.push(-galaxyRadius, 0, 0, galaxyRadius, 0, 0, 0, -galaxyRadius, 0, 0, galaxyRadius, 0);
  return out;
}

// The cross at the main view's center and the line to its eye, as segments: `target` and `eye` in parsecs, `arm` the cross's half-width.
export function markerSegments(target, eye, arm) {
  const [tx, ty, tz] = target;
  return [
    tx - arm, ty, tz, tx + arm, ty, tz,
    tx, ty - arm, tz, tx, ty + arm, tz,
    tx, ty, tz - arm, tx, ty, tz + arm,
    tx, ty, tz, eye[0], eye[1], eye[2],
  ];
}

// The mini map on `renderer`. `options`: galaxyRadius (pc), fovDeg, accent and outline colours (hex numbers).
export function createMiniMap(THREE, renderer, options) {
  const radius = options.galaxyRadius;
  const fovDeg = options.fovDeg || 50;
  const distance = MC.clampDistance(MINI_POLICY, 0, miniFitDistance(radius, fovDeg));
  const scene = new THREE.Scene();
  const camera = new THREE.PerspectiveCamera(fovDeg, 1, distance * 0.01, distance * 4);
  camera.up.set(0, 0, 1);

  const outlineGeometry = new THREE.BufferGeometry();
  outlineGeometry.setAttribute("position", new THREE.Float32BufferAttribute(outlineSegments(radius), 3));
  scene.add(new THREE.LineSegments(outlineGeometry, new THREE.LineBasicMaterial({ color: options.outline || 0x8a93b8, transparent: true, opacity: 0.7, depthTest: false })));

  const markerPositions = new Float32Array(8 * 3);
  const markerGeometry = new THREE.BufferGeometry();
  markerGeometry.setAttribute("position", new THREE.BufferAttribute(markerPositions, 3));
  const marker = new THREE.LineSegments(markerGeometry, new THREE.LineBasicMaterial({ color: options.accent || 0xffd24a, depthTest: false }));
  marker.frustumCulled = false;
  scene.add(marker);

  const panel = new THREE.Color(options.panel || 0x05070d);
  const saved = new THREE.Color();
  const sizeHolder = new THREE.Vector2();
  let signature = "";
  let drawn = null;

  // Draws the mini map over the canvas the main pass just drew; `orbit` and `target` are the main camera's, `eye` its position.
  function render(orbit, target, eye) {
    renderer.getSize(sizeHolder);
    const rect = miniRect(sizeHolder.x, sizeHolder.y);
    drawn = rect;
    if (!rect) return;
    // Dirty flag (section 5.5): the marker's geometry and the camera only change when the main view does.
    const next = [orbit.theta, orbit.phi, target.x, target.y, target.z, eye.x, eye.y, eye.z].join(",");
    if (next !== signature) {
      signature = next;
      const placed = miniEyeAt(orbit, distance);
      camera.position.set(placed[0], placed[1], placed[2]);
      camera.lookAt(0, 0, 0);
      camera.updateMatrixWorld();
      markerPositions.set(markerSegments([target.x, target.y, target.z], [eye.x, eye.y, eye.z], radius * 0.04));
      markerGeometry.attributes.position.needsUpdate = true;
    }
    const wasAutoClear = renderer.autoClear;
    const wasAlpha = renderer.getClearAlpha();
    renderer.getClearColor(saved);
    renderer.autoClear = false;
    renderer.setScissorTest(true);
    renderer.setViewport(rect.x, rect.y, rect.side, rect.side);
    renderer.setScissor(rect.x, rect.y, rect.side, rect.side);
    renderer.setClearColor(panel, 0.72);
    renderer.clear(true, true, false);
    renderer.render(scene, camera);
    renderer.setScissorTest(false);
    renderer.setViewport(0, 0, sizeHolder.x, sizeHolder.y);
    renderer.setClearColor(saved, wasAlpha);
    renderer.autoClear = wasAutoClear;
  }

  return {
    render: render,
    camera: camera,
    distance: distance,
    // For the browser tests: the square last drawn (null when the canvas had no room) and where the mini camera stands.
    state: function () { return { rect: drawn, eye: camera.position.toArray(), distance: distance }; },
    dispose: function () { outlineGeometry.dispose(); markerGeometry.dispose(); },
  };
}
