// html/static/nebulavolumelayer.js
//
// MAP.365: the nebulae of the Galaxy Map drawn from their shape, the drawing half (nebulavolume.js chooses and decodes). Each held
// node's points (soft splats, pixel-sized in the vertex shader like the glow's quads) are added up into a half-resolution
// render target, with no depth test: where nebulae overlap their colors add. A full-screen pass then bends that sum under the
// cap (`cappedSum`, the same formula, in the shader) and adds it to the map, so a crowd of nebulae never washes to white and the
// stars stay readable. `prepare` fills the target before the map is drawn; the pass is part of the map's scene.
//
// A layer is `null` where the browser cannot render to a float target; the map then keeps the nebulae it had.

const VERSION_QUERY = new URL(import.meta.url).search;
const THREE = await import(`./vendor/three.module.min.js${VERSION_QUERY}`);
const {
  CAP, POINT_ALPHA, SPLAT_STEPS, MIN_SPLAT_PX, MAX_SPLAT_PX, positionsOf, centerOf, wanted, volumeAlpha, createVolumeStore, MAX_LOADED,
} = await import(`./nebulavolume.js${VERSION_QUERY}`);

const CORNERS = [[-1, -1], [1, -1], [1, 1], [-1, 1]];
// The target is this share of the canvas across: the nebulae are soft, and the splats are cheap to add up at a quarter of the pixels.
export const TARGET_SCALE = 0.5;
const MAX_FAMILIES = 8;

// Whether `renderer` can render to a half-float target (and blend into it).
export function canRenderFloat(renderer) {
  const has = (name) => renderer.extensions && renderer.extensions.has(name);
  return !!(renderer.capabilities && renderer.capabilities.isWebGL2 && (has("EXT_color_buffer_float") || has("EXT_color_buffer_half_float")));
}

// A layer in `scene` for `camera` on `canvasEl` drawn by `renderer`. `options`: `colours` (the families' [r, g, b], 0 to 1, in the order of the
// volume's family byte), `fetchBytes(key) -> Promise<ArrayBuffer>`, `renderOrder` (default 6: over the stars' 5) and `onChange()`.
export function createNebulaVolumeLayer(scene, camera, canvasEl, renderer, options) {
  if (!canRenderFloat(renderer)) return null;
  const o = options || {};
  const familyColours = o.colours || [];
  const familyOn = familyColours.map(() => 1);
  const accumulation = new THREE.Scene();
  let target = null;
  let targetWidth = 0;
  let targetHeight = 0;
  let on = true;
  let nodeAlpha = new Map();
  let lastDrawn = 0;
  const held = new Map();       // node key -> { mesh, volume }
  const store = createVolumeStore(o.fetchBytes, () => { if (o.onChange) o.onChange(); });

  const familyUniform = Array.from({ length: MAX_FAMILIES }, () => new THREE.Vector4(0, 0, 0, 0));
  function writeFamilies() {
    familyColours.forEach((colour, i) => { if (i < MAX_FAMILIES) familyUniform[i].set(colour[0], colour[1], colour[2], familyOn[i]); });
  }
  writeFamilies();

  function splatMaterial(volume) {
    return new THREE.ShaderMaterial({
      uniforms: {
        viewport: { value: new THREE.Vector2(800, 600) },
        splat: { value: SPLAT_STEPS * volume.spacingPc },
        focal: { value: 800 },
        alpha: { value: 0 },
        families: { value: familyUniform },
        // MAP.222: a group window cuts the shape to its sphere: the center in the mesh's frame, and the radius (0 for none).
        clipCenter: { value: new THREE.Vector3() },
        clipRadius: { value: 0 },
      },
      vertexShader: [
        "attribute vec2 corner;",
        "attribute float strength;",
        "attribute float family;",
        "uniform vec2 viewport;",
        "uniform float splat;",
        "uniform float focal;",
        "uniform float alpha;",
        `uniform vec4 families[${MAX_FAMILIES}];`,
        "uniform vec3 clipCenter;",
        "uniform float clipRadius;",
        "varying vec2 vCorner;",
        "varying vec3 vColor;",
        "void main() {",
        "  vCorner = corner;",
        "  vec4 f = families[int(family + 0.5)];",
        `  vColor = f.rgb * (strength * ${POINT_ALPHA.toFixed(5)} * alpha * f.a);`,
        "  if (clipRadius > 0.0) vColor *= 1.0 - smoothstep(0.9 * clipRadius, clipRadius, distance(position, clipCenter));",
        "  vec4 clip = projectionMatrix * modelViewMatrix * vec4(position, 1.0);",
        `  float px = clamp(splat * focal / max(clip.w, 1e-6), ${MIN_SPLAT_PX.toFixed(1)}, ${MAX_SPLAT_PX.toFixed(1)});`,
        "  clip.xy += corner * (0.5 * px) * 2.0 / viewport * clip.w;",
        "  gl_Position = clip;",
        "}",
      ].join("\n"),
      fragmentShader: [
        "varying vec2 vCorner;",
        "varying vec3 vColor;",
        "void main() {",
        "  float r = length(vCorner);",
        "  if (r >= 1.0) discard;",
        "  float soft = exp(-r * r * 3.0) * (1.0 - r);",
        "  gl_FragColor = vec4(vColor * soft, 1.0);",
        "}",
      ].join("\n"),
      transparent: true,
      depthWrite: false,
      depthTest: false,
      blending: THREE.CustomBlending,
      blendEquation: THREE.AddEquation,
      blendSrc: THREE.OneFactor,
      blendDst: THREE.OneFactor,
      blendEquationAlpha: THREE.AddEquation,
      blendSrcAlpha: THREE.OneFactor,
      blendDstAlpha: THREE.OneFactor,
    });
  }

  function meshOf(volume) {
    const n = volume.count;
    const frame = centerOf(volume).map(Math.round);
    const positions = new Float32Array(12 * n);
    const corners = new Float32Array(8 * n);
    const strengths = new Float32Array(4 * n);
    const families = new Float32Array(4 * n);
    const indices = new Uint32Array(6 * n);
    const base = positionsOf(volume, frame);
    for (let i = 0; i < n; i++) {
      for (let c = 0; c < 4; c++) {
        const v = 4 * i + c;
        positions.set(base.subarray(3 * i, 3 * i + 3), 3 * v);
        corners.set(CORNERS[c], 2 * v);
        strengths[v] = volume.strength[i];
        families[v] = volume.family[i];
      }
      indices.set([4 * i, 4 * i + 1, 4 * i + 2, 4 * i, 4 * i + 2, 4 * i + 3], 6 * i);
    }
    const geometry = new THREE.BufferGeometry();
    geometry.setAttribute("position", new THREE.BufferAttribute(positions, 3));
    geometry.setAttribute("corner", new THREE.BufferAttribute(corners, 2));
    geometry.setAttribute("strength", new THREE.BufferAttribute(strengths, 1));
    geometry.setAttribute("family", new THREE.BufferAttribute(families, 1));
    geometry.setIndex(new THREE.BufferAttribute(indices, 1));
    const mesh = new THREE.Mesh(geometry, splatMaterial(volume));
    mesh.position.set(frame[0], frame[1], frame[2]);
    mesh.frustumCulled = false;
    mesh.updateMatrixWorld();
    return mesh;
  }

  // The pass that bends the sum under the cap and adds it to the map.
  const composite = new THREE.Mesh(
    new THREE.PlaneGeometry(2, 2),
    new THREE.ShaderMaterial({
      uniforms: { sum: { value: null }, cap: { value: CAP } },
      vertexShader: "varying vec2 vUv;\nvoid main() { vUv = uv; gl_Position = vec4(position.xy, 0.0, 1.0); }",
      fragmentShader: [
        "uniform sampler2D sum;",
        "uniform float cap;",
        "varying vec2 vUv;",
        "void main() {",
        "  vec3 s = texture2D(sum, vUv).rgb;",
        "  float m = max(s.r, max(s.g, s.b));",
        "  if (m <= 0.0) discard;",
        "  float bent = cap * (1.0 - exp(-m / cap));",
        "  gl_FragColor = vec4(s * (bent / m), 1.0);",
        "}",
      ].join("\n"),
      transparent: true,
      depthWrite: false,
      depthTest: false,
      blending: THREE.CustomBlending,
      blendEquation: THREE.AddEquation,
      blendSrc: THREE.OneFactor,
      blendDst: THREE.OneFactor,
      blendEquationAlpha: THREE.AddEquation,
      blendSrcAlpha: THREE.ZeroFactor,
      blendDstAlpha: THREE.OneFactor,
    }),
  );
  composite.frustumCulled = false;
  composite.renderOrder = o.renderOrder === undefined ? 6 : o.renderOrder;
  composite.visible = false;
  scene.add(composite);
  // The same pass in a scene of its own, so a test can draw it alone and read what it makes.
  const probeScene = new THREE.Scene();
  const probe = new THREE.Mesh(composite.geometry, composite.material);
  probe.frustumCulled = false;
  probeScene.add(probe);

  function ensureTarget() {
    const width = Math.max(2, Math.floor(Math.max(1, canvasEl.clientWidth) * TARGET_SCALE));
    const height = Math.max(2, Math.floor(Math.max(1, canvasEl.clientHeight) * TARGET_SCALE));
    if (target && width === targetWidth && height === targetHeight) return;
    if (target) target.dispose();
    target = new THREE.WebGLRenderTarget(width, height, { type: THREE.HalfFloatType, format: THREE.RGBAFormat, depthBuffer: false });
    targetWidth = width;
    targetHeight = height;
    composite.material.uniforms.sum.value = target.texture;
  }

  // MAP.222: the window's sphere, `{ center: [x, y, z], radiusPc }` or null; the shape is faded out to it, in each mesh's own frame.
  let clip = null;
  function applyClip(mesh) {
    const uniforms = mesh.material.uniforms;
    if (!clip) {
      uniforms.clipRadius.value = 0;
      return;
    }
    uniforms.clipCenter.value.set(clip.center[0] - mesh.position.x, clip.center[1] - mesh.position.y, clip.center[2] - mesh.position.z);
    uniforms.clipRadius.value = clip.radiusPc;
  }
  function setClip(sphere) {
    clip = sphere ? { center: sphere.center.slice(), radiusPc: sphere.radiusPc } : null;
    held.forEach((entry) => applyClip(entry.mesh));
  }

  // Chooses which nodes to hold (`sources` from nebulavolume.sourceOf, `makeWeigh()` gives the cut's `weigh(block) -> 0..1` share, called only
  // when something is held, or is null for whole shares),
  // asks for the ones missing and sets how strongly each shows. Call it every frame the camera or the cut may have moved.
  function update(sources, makeWeigh, focalPx) {
    const position = camera.position.toArray();
    const found = on ? wanted(sources, position, focalPx) : [];
    const keep = new Set();
    nodeAlpha = new Map();
    const weigh = found.length && makeWeigh ? makeWeigh() : null;
    found.forEach((item) => {
      const key = item.source.key;
      keep.add(key);
      const state = store.request(key);
      if (state.status !== "ready" || !state.volume || !state.volume.count) return;
      let entry = held.get(key);
      if (!entry || entry.volume !== state.volume) {
        if (entry) { accumulation.remove(entry.mesh); entry.mesh.geometry.dispose(); entry.mesh.material.dispose(); }
        entry = { mesh: meshOf(state.volume), volume: state.volume };
        held.set(key, entry);
        accumulation.add(entry.mesh);
      }
      const share = weigh ? weigh(item.source.block) : 1;
      const alpha = volumeAlpha(item.radiusPx, share);
      nodeAlpha.set(key, alpha);
      const uniforms = entry.mesh.material.uniforms;
      uniforms.alpha.value = alpha;
      uniforms.focal.value = focalPx;
      uniforms.viewport.value.set(Math.max(1, canvasEl.clientWidth), Math.max(1, canvasEl.clientHeight));
      applyClip(entry.mesh);
    });
    held.forEach((entry, key) => {
      if (!keep.has(key)) {
        accumulation.remove(entry.mesh);
        entry.mesh.geometry.dispose();
        entry.mesh.material.dispose();
        held.delete(key);
      }
    });
    store.release(keep);
    held.forEach((entry, key) => { entry.mesh.visible = (nodeAlpha.get(key) || 0) > 0.002; });
    lastDrawn = Array.from(held.keys()).filter((key) => (nodeAlpha.get(key) || 0) > 0.002).length;
    composite.visible = on && lastDrawn > 0;
  }

  // Fills the half-resolution target with the sum of the splats; call it before the map is drawn.
  const savedColour = new THREE.Color();
  function prepare() {
    if (!composite.visible) return;
    ensureTarget();
    const savedAlpha = renderer.getClearAlpha();
    renderer.getClearColor(savedColour);
    const previous = renderer.getRenderTarget();
    renderer.setRenderTarget(target);
    renderer.setClearColor(0x000000, 0);
    renderer.clear();
    renderer.render(accumulation, camera);
    renderer.setRenderTarget(previous);
    renderer.setClearColor(savedColour, savedAlpha);
  }

  // A half-float number (its 16 bits) as a JavaScript number.
  function halfToFloat(bits) {
    const sign = bits & 0x8000 ? -1 : 1;
    const exponent = (bits >> 10) & 0x1f;
    const fraction = bits & 0x3ff;
    if (exponent === 0) return sign * Math.pow(2, -14) * (fraction / 1024);
    if (exponent === 31) return fraction ? NaN : sign * Infinity;
    return sign * Math.pow(2, exponent - 15) * (1 + fraction / 1024);
  }

  // For the browser tests: the target as it is now, `{ width, height, sums (Float32Array r, g, b per pixel: what the splats add up to),
  // capped (Uint8Array r, g, b, a per pixel: what the cap pass makes of it, drawn alone) }`, or null when nothing is drawn.
  function snapshot() {
    if (!composite.visible) return null;
    prepare();
    const width = targetWidth;
    const height = targetHeight;
    const half = new Uint16Array(width * height * 4);
    renderer.readRenderTargetPixels(target, 0, 0, width, height, half);
    const sums = new Float32Array(width * height * 3);
    for (let i = 0; i < width * height; i++) {
      for (let c = 0; c < 3; c++) sums[3 * i + c] = halfToFloat(half[4 * i + c]);
    }
    const alone = new THREE.WebGLRenderTarget(width, height, { type: THREE.UnsignedByteType, depthBuffer: false });
    const previous = renderer.getRenderTarget();
    const savedAlpha = renderer.getClearAlpha();
    renderer.getClearColor(savedColour);
    renderer.setRenderTarget(alone);
    renderer.setClearColor(0x000000, 0);
    renderer.clear();
    renderer.render(probeScene, camera);
    const capped = new Uint8Array(width * height * 4);
    renderer.readRenderTargetPixels(alone, 0, 0, width, height, capped);
    renderer.setRenderTarget(previous);
    renderer.setClearColor(savedColour, savedAlpha);
    alone.dispose();
    return { width, height, sums, capped };
  }

  function setFamilyVisible(index, value) {
    if (index >= 0 && index < familyOn.length) {
      familyOn[index] = value ? 1 : 0;
      writeFamilies();
    }
  }

  function setVisible(value) {
    on = !!value;
    composite.visible = on && lastDrawn > 0;
  }

  return {
    update, prepare, snapshot, setFamilyVisible, setVisible, setClip,
    // How strongly node `key`'s volume shows (0 until it has arrived): the share of its nebulae the sprites and meshes give up.
    coverage: (key) => nodeAlpha.get(key) || 0,
    // What is drawn, for the browser tests: held volumes, those showing, the points in them.
    heldCount: () => held.size,
    drawn: () => lastDrawn,
    // The held volumes as plain arrays, for the browser tests: `[{ key, level, spacingPc, ix, iy, iz, strength, family }]`.
    volumes: () => Array.from(held.entries()).map(([key, entry]) => ({
      key, level: entry.volume.level, spacingPc: entry.volume.spacingPc, ix: Array.from(entry.volume.ix), iy: Array.from(entry.volume.iy),
      iz: Array.from(entry.volume.iz), strength: Array.from(entry.volume.strength), family: Array.from(entry.volume.family),
    })),
    pointCount: () => Array.from(held.values()).reduce((sum, entry) => sum + entry.volume.count, 0),
    isVisible: () => composite.visible,
    // The raw sum of the nebulae (a texture, or null before anything was drawn), which stars and glows read to take a nebula's hue (MAP.215).
    texture: () => (target ? target.texture : null),
    capLimit: CAP,
    maxLoaded: MAX_LOADED,
  };
}
