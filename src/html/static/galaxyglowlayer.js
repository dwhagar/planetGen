// html/static/galaxyglowlayer.js
//
// MAP.266: the unresolved glow of the Galaxy Map, the drawing half (the choosing is galaxyglow.js). One screen-facing
// quad per blob, soft and round, in the color of the light's mean temperature, drawn under the stars with no depth test.
// The quad is sized in pixels in the vertex shader, so a blob as wide as the screen is no trouble (a point sprite is
// capped by the driver). `setBlobs` takes what `blobsOf` made; `update` sets each blob's strength from the share the cut
// gives its node and from its apparent brightness (call it every frame the camera or the cut may have moved).

const VERSION_QUERY = new URL(import.meta.url).search;
const THREE = await import(`./vendor/three.module.min.js${VERSION_QUERY}`);
const { glowLook } = await import(`./galaxyglow.js${VERSION_QUERY}`);

const CORNERS = [[-1, -1], [1, -1], [1, 1], [-1, 1]];

// A star's color from its surface temperature (Tanner Helland's fit), as [r, g, b] in 0..1.
function temperatureColor(temperatureK) {
  const t = Math.min(40000, Math.max(1000, temperatureK || 5800)) / 100;
  const r = t <= 66 ? 255 : 329.698727446 * Math.pow(t - 60, -0.1332047592);
  const g = t <= 66 ? 99.4708025861 * Math.log(t) - 161.1195681661 : 288.1221695283 * Math.pow(t - 60, -0.0755148492);
  const b = t >= 66 ? 255 : t <= 19 ? 0 : 138.5177312231 * Math.log(t - 10) - 305.0447927307;
  return [r, g, b].map((v) => Math.min(255, Math.max(0, v)) / 255);
}

// A glow layer in `scene` for `camera` on `canvasEl`. `options`: `light` (true on a light background, where the glow is
// drawn dark-on-light by normal blending instead of added), `renderOrder` (default 4: under the stars' 5), `look`
// (`(blob, camera, focalPx) -> {alpha, radiusPx}`, default `glowLook`: the charted regions, MAP.316, are drawn by the same layer
// with their own), `cut` (`(alphas) -> alphas`, changed in place: the strengths of every blob, which a budget may hold down, MAP.176;
// default none) and a blob's own `color` ([r, g, b], 0 to 1) in place of the color of its temperature and its `ring` (true: a thin
// line near the rim).
export function createGlowLayer(scene, camera, canvasEl, options) {
  const o = options || {};
  const lookOf = o.look || glowLook;
  const cutOf = o.cut || null;
  const material = new THREE.ShaderMaterial({
    uniforms: { viewport: { value: new THREE.Vector2(800, 600) } },
    vertexShader: [
      "attribute vec2 corner;",
      "attribute vec3 glowColor;",
      "attribute float glowRadius;",
      "attribute float glowAlpha;",
      "attribute float glowRing;",
      "uniform vec2 viewport;",
      "varying vec2 vCorner;",
      "varying vec3 vColor;",
      "varying float vAlpha;",
      "varying float vRing;",
      "void main() {",
      "  vCorner = corner;",
      "  vRing = glowRing;",
      "  vColor = glowColor;",
      "  vAlpha = glowAlpha;",
      "  vec4 clip = projectionMatrix * modelViewMatrix * vec4(position, 1.0);",
      "  clip.xy += corner * glowRadius * 2.0 / viewport * clip.w;",
      "  gl_Position = clip;",
      "}",
    ].join("\n"),
    fragmentShader: [
      "varying vec2 vCorner;",
      "varying vec3 vColor;",
      "varying float vAlpha;",
      "varying float vRing;",
      "void main() {",
      "  float r = length(vCorner);",
      "  if (r >= 1.0) discard;",
      // A soft round falloff: bright in the middle, nothing at the rim.
      "  float soft = exp(-r * r * 3.5) * (1.0 - r);",
      // A blob with a ring (a group of nebulae, MAP.176) also has a thin line near its rim, so it reads as a group.
      "  float line = vRing * smoothstep(0.78, 0.88, r) * (1.0 - smoothstep(0.88, 0.98, r)) * 0.9;",
      "  float a = max(soft, line) * vAlpha;",
      o.light ? "  gl_FragColor = vec4(vColor, a);" : "  gl_FragColor = vec4(vColor * a, a);",
      "}",
    ].join("\n"),
    transparent: true,
    depthWrite: false,
    depthTest: false,
  });
  // Where glows overlap the brighter one shows and they are not added (a crowd of overlapping blobs added up to white
  // and hid the stars, MAP.324): a maximum blend of the premultiplied colors on a dark background, normal blending
  // dark-on-light on a light one.
  if (o.light) {
    material.blending = THREE.NormalBlending;
  } else {
    material.blending = THREE.CustomBlending;
    material.blendEquation = THREE.MaxEquation;
    material.blendSrc = THREE.OneFactor;
    material.blendDst = THREE.OneFactor;
    material.blendEquationAlpha = THREE.MaxEquation;
    material.blendSrcAlpha = THREE.OneFactor;
    material.blendDstAlpha = THREE.OneFactor;
  }
  const mesh = new THREE.Mesh(new THREE.BufferGeometry(), material);
  mesh.renderOrder = o.renderOrder === undefined ? 4 : o.renderOrder;
  mesh.frustumCulled = false;
  scene.add(mesh);

  let blobs = [];
  let frame = [0, 0, 0];
  let alphaAttribute = null;
  let radiusAttribute = null;
  let ringAttribute = null;
  let on = true;
  let drawnCount = 0;
  let strongestAlpha = 0;
  let signature = "";

  // Draws exactly `list` (`blobsOf`'s). Positions are offsets from `centre` rounded to a whole parsec, the mesh sitting
  // there (camera-relative, as the star layer's, so a float32 holds them).
  function setBlobs(list, centre) {
    blobs = list;
    signature = "";
    frame = [Math.round(centre[0]), Math.round(centre[1]), Math.round(centre[2])];
    const n = list.length;
    const positions = new Float32Array(12 * n);
    const corners = new Float32Array(8 * n);
    const colors = new Float32Array(12 * n);
    const radii = new Float32Array(4 * n);
    const alphas = new Float32Array(4 * n);
    const rings = new Float32Array(4 * n);
    const indices = new Uint32Array(6 * n);
    list.forEach((blob, i) => {
      const color = blob.color || temperatureColor(blob.temperatureK);
      for (let c = 0; c < 4; c++) {
        const v = 4 * i + c;
        positions.set([blob.x - frame[0], blob.y - frame[1], blob.z - frame[2]], 3 * v);
        corners.set(CORNERS[c], 2 * v);
        colors.set(color, 3 * v);
        rings[v] = blob.ring ? 1 : 0;
      }
      indices.set([4 * i, 4 * i + 1, 4 * i + 2, 4 * i, 4 * i + 2, 4 * i + 3], 6 * i);
    });
    const geometry = new THREE.BufferGeometry();
    geometry.setAttribute("position", new THREE.BufferAttribute(positions, 3));
    geometry.setAttribute("corner", new THREE.BufferAttribute(corners, 2));
    geometry.setAttribute("glowColor", new THREE.BufferAttribute(colors, 3));
    radiusAttribute = new THREE.BufferAttribute(radii, 1);
    alphaAttribute = new THREE.BufferAttribute(alphas, 1);
    geometry.setAttribute("glowRadius", radiusAttribute);
    geometry.setAttribute("glowAlpha", alphaAttribute);
    ringAttribute = new THREE.BufferAttribute(rings, 1);
    geometry.setAttribute("glowRing", ringAttribute);
    geometry.setIndex(new THREE.BufferAttribute(indices, 1));
    mesh.geometry.dispose();
    mesh.geometry = geometry;
    mesh.position.set(frame[0], frame[1], frame[2]);
    mesh.updateMatrixWorld();
  }

  // Whether a blob of `radiusPx` reaches the screen (its centre in front of the camera, its disc over the canvas).
  const spot = new THREE.Vector3();
  function onScreen(blob, radiusPx) {
    spot.set(blob.x, blob.y, blob.z).project(camera);
    if (!(spot.z > -1 && spot.z < 1)) return false;
    const width = Math.max(1, canvasEl.clientWidth);
    const height = Math.max(1, canvasEl.clientHeight);
    const x = (spot.x + 1) / 2 * width;
    const y = (1 - spot.y) / 2 * height;
    return x > -radiusPx && x < width + radiusPx && y > -radiusPx && y < height + radiusPx;
  }

  // Sets each blob's radius and strength: `makeWeigh()` gives the cut's weight calculator for the camera (called only when something changed; `null` for whole shares), `focalPx` the
  // camera's focal length in pixels. Nothing is done while the camera, the cut (`stamp` names it) and the layer's state are as they were.
  function update(makeWeigh, focalPx, stamp) {
    const p = camera.position;
    const key = [on ? 1 : 0, blobs.length, p.x, p.y, p.z, canvasEl.clientWidth, canvasEl.clientHeight, focalPx, stamp].join(",");
    if (key === signature) return;
    signature = key;
    material.uniforms.viewport.value.set(Math.max(1, canvasEl.clientWidth), Math.max(1, canvasEl.clientHeight));
    mesh.visible = on && blobs.length > 0;
    drawnCount = 0;
    strongestAlpha = 0;
    if (!mesh.visible || !alphaAttribute) return;
    const position = camera.position.toArray();
    const weigh = makeWeigh ? makeWeigh() : null;
    const looks = [];
    let alphas = blobs.map((blob, i) => {
      const share = weigh ? weigh(blob.block) : 1;
      const look = lookOf(blob, position, focalPx);
      looks[i] = look;
      return share > 0 ? look.alpha * share : 0;
    });
    if (cutOf) alphas = cutOf(alphas);
    blobs.forEach((blob, i) => {
      const alpha = alphas[i];
      const look = looks[i];
      if (alpha > 0.002) drawnCount += 1;
      if (alpha > strongestAlpha && onScreen(blob, look.radiusPx)) strongestAlpha = alpha;
      for (let c = 0; c < 4; c++) {
        alphaAttribute.array[4 * i + c] = alpha;
        radiusAttribute.array[4 * i + c] = look.radiusPx;
      }
    });
    alphaAttribute.needsUpdate = true;
    radiusAttribute.needsUpdate = true;
  }

  function setVisible(value) {
    on = !!value;
    signature = "";
    mesh.visible = on && blobs.length > 0;
  }

  return {
    mesh, setBlobs, update, setVisible,
    // What is drawn, for the browser tests: how many blobs have a visible strength, the strongest of those that reach the
    // screen (an alpha, 0 to GLOW_ALPHA_MAX), how many are held, and whether the layer is on.
    drawn: () => drawnCount,
    strongest: () => strongestAlpha,
    count: () => blobs.length,
    isVisible: () => mesh.visible,
    blobs: () => blobs.slice(),
  };
}
