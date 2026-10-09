// html/static/phenomenonrender.js
//
// The phenomenon page's "View" panel (planetgen/web/maps/phenomenonrender.py): a
// three.js picture of a neutron star, black hole, quasar, rogue planet or
// interstellar comet, built from the numbers in `#phenomrender`'s
// `data-view` JSON. The panel's static SVG still stays when WebGL can't
// start; with `prefers-reduced-motion` the scene draws one still frame and
// only redraws while the visitor drags it.
//
// Sizes are illustrative, not to scale: the caption under the view says
// what is real (the neutron star's slowed spin rate, for one).

const VERSION_QUERY = new URL(import.meta.url).search;
const THREE = await import(`./vendor/three.module.min.js${VERSION_QUERY}`);
const { makeGlowMaterial, makeStarSurfaceTexture } = await import(`./bodyRendering.js${VERSION_QUERY}`);

var viewport = document.getElementById("phenomrender");
var canvas = document.getElementById("phenomrender-canvas");

function readView() {
  if (!viewport) {
    return null;
  }
  try {
    return JSON.parse(viewport.getAttribute("data-view"));
  } catch (err) {
    return null;
  }
}

var reducedMotion = typeof window.matchMedia === "function"
  && window.matchMedia("(prefers-reduced-motion: reduce)").matches;

// A small seeded random source, so a still frame is the same every visit.
function makeRandom(seed) {
  var state = seed >>> 0;
  return function () {
    state = (state * 1664525 + 1013904223) >>> 0;
    return state / 4294967296;
  };
}

// A soft round sprite texture: `stops` are [offset, css color] pairs.
function radialTexture(stops) {
  var size = 128;
  var el = document.createElement("canvas");
  el.width = el.height = size;
  var ctx = el.getContext("2d");
  var grad = ctx.createRadialGradient(size / 2, size / 2, 0, size / 2, size / 2, size / 2);
  stops.forEach(function (stop) { grad.addColorStop(stop[0], stop[1]); });
  ctx.fillStyle = grad;
  ctx.fillRect(0, 0, size, size);
  var texture = new THREE.CanvasTexture(el);
  texture.colorSpace = THREE.SRGBColorSpace;
  return texture;
}

function glowSprite(stops, scale, opacity) {
  var sprite = new THREE.Sprite(new THREE.SpriteMaterial({
    map: radialTexture(stops),
    blending: THREE.AdditiveBlending,
    transparent: true,
    depthWrite: false,
    opacity: opacity != null ? opacity : 1,
  }));
  sprite.scale.setScalar(scale);
  return sprite;
}

function starfield(random) {
  var count = 1200;
  var positions = new Float32Array(count * 3);
  for (var i = 0; i < count; i++) {
    var z = random() * 2 - 1;
    var angle = random() * Math.PI * 2;
    var r = Math.sqrt(1 - z * z) * 300;
    positions[i * 3] = Math.cos(angle) * r;
    positions[i * 3 + 1] = Math.sin(angle) * r;
    positions[i * 3 + 2] = z * 300;
  }
  var geometry = new THREE.BufferGeometry();
  geometry.setAttribute("position", new THREE.BufferAttribute(positions, 3));
  return new THREE.Points(geometry, new THREE.PointsMaterial({
    color: 0xc8d4ff, size: 1.4, sizeAttenuation: false, transparent: true, opacity: 0.8,
  }));
}

// --- Accretion disk ---------------------------------------------------------
// A flat ring whose color runs from `hot` at the inner edge to `cool` at
// the outer one, streaked by spiral bands that turn faster near the
// middle (Keplerian: angular speed falls as r^-1.5).
var DISK_VERTEX = [
  "varying vec2 vPos;",
  "void main() {",
  "  vPos = position.xy;",
  "  gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);",
  "}",
].join("\n");

var DISK_FRAGMENT = [
  "uniform float time;",
  "uniform float inner;",
  "uniform float outer;",
  "uniform vec3 hot;",
  "uniform vec3 cool;",
  "varying vec2 vPos;",
  "void main() {",
  "  float r = length(vPos);",
  "  float u = clamp((r - inner) / (outer - inner), 0.0, 1.0);",
  "  float angle = atan(vPos.y, vPos.x);",
  "  float spin = time * 1.6 * pow(inner / r, 1.5);",
  "  float bands = 0.65 + 0.35 * sin(angle * 3.0 + log(r) * 9.0 - spin * 3.0)",
  "    * sin(angle * 7.0 - log(r) * 4.0 - spin * 7.0 + 1.3);",
  "  vec3 color = mix(hot, cool, pow(u, 0.6)) * bands;",
  "  float edge = smoothstep(0.0, 0.08, u) * (1.0 - smoothstep(0.55, 1.0, u));",
  "  gl_FragColor = vec4(color * (1.6 - u), edge);",
  "}",
].join("\n");

function accretionDisk(inner, outer, hot, cool) {
  var material = new THREE.ShaderMaterial({
    uniforms: {
      time: { value: 0 },
      inner: { value: inner },
      outer: { value: outer },
      hot: { value: new THREE.Color(hot) },
      cool: { value: new THREE.Color(cool) },
    },
    vertexShader: DISK_VERTEX,
    fragmentShader: DISK_FRAGMENT,
    side: THREE.DoubleSide,
    blending: THREE.AdditiveBlending,
    transparent: true,
    depthWrite: false,
  });
  var disk = new THREE.Mesh(new THREE.RingGeometry(inner, outer, 160, 12), material);
  disk.rotation.x = -Math.PI / 2;
  return disk;
}

// A beam or jet: an open cone, apex at the origin, whose light fades
// from the apex to the far end.
var BEAM_VERTEX = [
  "varying float vAlong;",
  "uniform float beamLength;",
  "void main() {",
  "  vAlong = clamp(1.0 - position.y / beamLength, 0.0, 1.0);",
  "  gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);",
  "}",
].join("\n");

var BEAM_FRAGMENT = [
  "uniform vec3 color;",
  "uniform float strength;",
  "varying float vAlong;",
  "void main() {",
  "  gl_FragColor = vec4(color, strength * pow(vAlong, 1.6));",
  "}",
].join("\n");

// A cone from the origin out along +y, `length` long.
function lightCone(length, radius, color, strength) {
  var geometry = new THREE.ConeGeometry(radius, length, 40, 1, true);
  geometry.translate(0, -length / 2, 0);
  geometry.rotateX(Math.PI);
  var material = new THREE.ShaderMaterial({
    uniforms: {
      color: { value: new THREE.Color(color) },
      strength: { value: strength },
      beamLength: { value: length },
    },
    vertexShader: BEAM_VERTEX,
    fragmentShader: BEAM_FRAGMENT,
    side: THREE.DoubleSide,
    blending: THREE.AdditiveBlending,
    transparent: true,
    depthWrite: false,
  });
  return new THREE.Mesh(geometry, material);
}

// --- Builders: each returns {object, cameraDistance, update(t, camera)} ------

function neutronStar(view) {
  var root = new THREE.Group();
  // The spin axis leans toward the viewer, so a beam 35 degrees off it
  // sweeps straight past the camera once a turn.
  var MAGNETIC_TILT = 35 * Math.PI / 180;
  root.rotation.x = Math.PI / 2 - MAGNETIC_TILT;
  var spinner = new THREE.Group();
  root.add(spinner);

  var surface = new THREE.Mesh(
    new THREE.SphereGeometry(1, 48, 32),
    new THREE.MeshBasicMaterial({ map: makeStarSurfaceTexture(THREE, "#cfe8ff") }),
  );
  spinner.add(surface);
  var glowColor = view.magnetar ? "#c9b0ff" : "#9fd0ff";
  var glow = new THREE.Mesh(new THREE.SphereGeometry(1.8, 32, 16), makeGlowMaterial(THREE, glowColor, 2.2, 1.2, 1.8));
  root.add(glow);

  var magnetic = new THREE.Group();
  magnetic.rotation.z = MAGNETIC_TILT;
  spinner.add(magnetic);
  var beams = [];
  if (view.pulsing) {
    [1, -1].forEach(function (sign) {
      var beam = lightCone(9, 1.1, "#bfe2ff", 0.55);
      if (sign < 0) {
        beam.rotation.z = Math.PI;
      }
      magnetic.add(beam);
      beams.push(beam);
    });
  }
  // Hot spots at the magnetic poles: they show the spin even with no beams.
  [1, -1].forEach(function (sign) {
    var spot = glowSprite([[0, "rgba(255,255,255,1)"], [0.3, "rgba(190,225,255,0.7)"], [1, "rgba(120,180,255,0)"]], 0.7);
    spot.position.set(0, sign * 1.02, 0);
    magnetic.add(spot);
  });
  var flash = glowSprite([[0, "rgba(255,255,255,1)"], [0.2, "rgba(200,230,255,0.8)"], [1, "rgba(120,180,255,0)"]], 3, 0);
  root.add(flash);

  var omega = Math.PI * 2 / Math.max(0.2, view.period_s || 2);
  var beamDir = new THREE.Vector3();
  var toCamera = new THREE.Vector3();
  return {
    object: root,
    cameraDistance: 9,
    update: function (t, camera) {
      spinner.rotation.y = t * omega;
      if (!beams.length) {
        return;
      }
      toCamera.copy(camera.position).normalize();
      var best = 0;
      beams.forEach(function (beam) {
        beamDir.set(0, 1, 0).transformDirection(beam.matrixWorld);
        best = Math.max(best, beamDir.dot(toCamera));
      });
      var pulse = Math.pow(Math.max(0, best), 24);
      flash.material.opacity = pulse * 0.9;
      flash.scale.setScalar(3 + pulse * 5);
    },
  };
}

function blackHole(view) {
  var root = new THREE.Group();
  var shadow = new THREE.Mesh(new THREE.SphereGeometry(1, 48, 32), new THREE.MeshBasicMaterial({ color: 0x000000 }));
  root.add(shadow);
  var disk = null;
  // A smaller hole's inner disk is hotter: stellar ones read blue-white,
  // intermediate ones yellow-white.
  var hot = view.stellar ? "#e8f0ff" : "#fff2d0";
  var cool = view.stellar ? "#ff7a2a" : "#c0401a";
  if (view.disk) {
    disk = accretionDisk(2.0, view.stellar ? 6.0 : 7.5, hot, cool);
    root.add(disk);
  }
  // The photon ring, and the disk's far side bent up over the shadow.
  var ringColor = view.disk ? "rgba(255,210,150," : "rgba(200,215,255,";
  var halo = glowSprite([[0, ringColor + "0)"], [0.36, ringColor + "0)"], [0.42, ringColor + "0.95)"],
    [0.5, ringColor + "0.25)"], [1, ringColor + "0)"]], 5.4, view.disk ? 1 : 0.6);
  root.add(halo);
  root.rotation.x = 0.22;
  root.rotation.z = 0.12 * Math.sign(view.spin || 1);
  return {
    object: root,
    cameraDistance: 13,
    update: function (t) {
      if (disk) {
        disk.material.uniforms.time.value = t;
      }
    },
  };
}

function quasar(view) {
  var root = new THREE.Group();
  var disk = accretionDisk(1.2, 6, "#ffffff", "#ff9a4a");
  root.add(disk);
  root.add(glowSprite([[0, "rgba(255,255,255,1)"], [0.15, "rgba(255,245,220,0.9)"], [0.5, "rgba(255,200,140,0.2)"], [1, "rgba(255,180,120,0)"]], 9));
  var torus = new THREE.Mesh(
    new THREE.TorusGeometry(9, 2.6, 24, 96),
    new THREE.MeshStandardMaterial({ color: 0x3a2418, roughness: 1, transparent: true, opacity: 0.85 }),
  );
  torus.rotation.x = Math.PI / 2;
  root.add(torus);
  var light = new THREE.PointLight(0xfff0d8, 60, 0, 1.4);
  root.add(light);
  var knots = [];
  if (view.jets) {
    [1, -1].forEach(function (sign) {
      var jet = lightCone(22, 1.4, "#a8c8ff", 0.5);
      if (sign < 0) {
        jet.rotation.z = Math.PI;
      }
      root.add(jet);
      for (var i = 0; i < 5; i++) {
        var knot = glowSprite([[0, "rgba(220,235,255,0.9)"], [1, "rgba(150,190,255,0)"]], 1.6);
        knot.userData = { sign: sign, phase: i / 5 };
        root.add(knot);
        knots.push(knot);
      }
    });
  }
  root.rotation.x = 0.35;
  root.rotation.z = -0.2;
  return {
    object: root,
    cameraDistance: 34,
    update: function (t) {
      disk.material.uniforms.time.value = t;
      knots.forEach(function (knot) {
        var u = (knot.userData.phase + t * 0.08) % 1;
        knot.position.set(0, knot.userData.sign * (1.5 + u * 20), 0);
        knot.material.opacity = 1 - u;
      });
    },
  };
}

function planetTexture(view, random) {
  var w = 256;
  var h = 128;
  var el = document.createElement("canvas");
  el.width = w;
  el.height = h;
  var ctx = el.getContext("2d");
  if (view.gas) {
    for (var y = 0; y < h; y++) {
      var band = Math.sin(y * 0.21) * 0.5 + Math.sin(y * 0.07 + 1.3) * 0.5;
      var l = 34 + band * 10 + random() * 3;
      ctx.fillStyle = "hsl(212, 22%, " + l + "%)";
      ctx.fillRect(0, y, w, 1);
    }
  } else {
    ctx.fillStyle = "hsl(210, 8%, 46%)";
    ctx.fillRect(0, 0, w, h);
    for (var i = 0; i < 260; i++) {
      ctx.fillStyle = "hsla(" + (200 + random() * 30) + ", 12%, " + (30 + random() * 40) + "%, 0.5)";
      ctx.beginPath();
      ctx.arc(random() * w, random() * h, 2 + random() * 9, 0, Math.PI * 2);
      ctx.fill();
    }
  }
  var texture = new THREE.CanvasTexture(el);
  texture.colorSpace = THREE.SRGBColorSpace;
  return texture;
}

// Glowing cracks (rocky) or a deep red underglow (gas) for internal heat.
function heatTexture(view, random) {
  var el = document.createElement("canvas");
  el.width = 512;
  el.height = 256;
  var ctx = el.getContext("2d");
  ctx.fillStyle = view.gas ? "#2a0805" : "#000000";
  ctx.fillRect(0, 0, 512, 256);
  if (!view.gas) {
    ctx.strokeStyle = "#ff5a1e";
    ctx.lineWidth = 2;
    ctx.lineJoin = "round";
    ctx.shadowColor = "#ff3a0a";
    ctx.shadowBlur = 4;
    for (var i = 0; i < 26; i++) {
      var x = random() * 512;
      var y = 32 + random() * 192;
      ctx.beginPath();
      ctx.moveTo(x, y);
      for (var s = 0; s < 6; s++) {
        x += (random() - 0.5) * 60;
        y += (random() - 0.5) * 32;
        ctx.lineTo(x, y);
      }
      ctx.stroke();
    }
  }
  var texture = new THREE.CanvasTexture(el);
  texture.colorSpace = THREE.SRGBColorSpace;
  return texture;
}

function roguePlanet(view, random) {
  var root = new THREE.Group();
  var material = new THREE.MeshStandardMaterial({ map: planetTexture(view, random), roughness: 0.95 });
  if (view.heat) {
    material.emissiveMap = heatTexture(view, random);
    material.emissive = new THREE.Color(0xffffff);
    material.emissiveIntensity = 0.9;
  }
  var body = new THREE.Mesh(new THREE.SphereGeometry(2, 64, 32), material);
  body.rotation.z = 0.4;
  root.add(body);
  if (view.heat) {
    root.add(new THREE.Mesh(new THREE.SphereGeometry(2.25, 32, 16), makeGlowMaterial(THREE, "#ff5a28", 3, 0.35, 1.125)));
  }
  root.add(new THREE.AmbientLight(0x8090b0, 0.35));
  var starlight = new THREE.DirectionalLight(0xc8d4ff, 0.9);
  starlight.position.set(-5, 2, 1);
  root.add(starlight);
  return {
    object: root,
    cameraDistance: 7.5,
    update: function (t) {
      body.rotation.y = t * 0.12;
    },
  };
}

function interstellarComet(view, random) {
  var root = new THREE.Group();
  var geometry = new THREE.IcosahedronGeometry(0.6, 4);
  var pos = geometry.attributes.position;
  var v = new THREE.Vector3();
  // Lumpy: stretch along x and push each vertex by a few smooth waves.
  var a = random() * 6;
  var b = random() * 6;
  for (var i = 0; i < pos.count; i++) {
    v.fromBufferAttribute(pos, i);
    var bump = 1 + 0.16 * Math.sin(v.x * 7 + a) * Math.cos(v.y * 6 + b) + 0.08 * Math.sin(v.z * 11 + a + b);
    v.multiplyScalar(bump);
    v.x *= 1.5;
    pos.setXYZ(i, v.x, v.y, v.z);
  }
  geometry.computeVertexNormals();
  var nucleus = new THREE.Mesh(geometry, new THREE.MeshStandardMaterial({ color: 0x55504a, roughness: 1 }));
  root.add(nucleus);
  root.add(new THREE.AmbientLight(0x8090b0, 0.4));
  var light = new THREE.DirectionalLight(0xffffff, 1.4);
  light.position.set(-4, 1.5, 3);
  root.add(light);

  var particles = null;
  var seeds = null;
  var TAIL = 14;
  if (view.active) {
    root.add(glowSprite([[0, "rgba(230,250,255,0.9)"], [0.25, "rgba(170,220,240,0.45)"], [1, "rgba(140,200,230,0)"]], 4.5));
    var ion = lightCone(TAIL * 1.3, 0.5, "#7fb4ff", 0.35);
    ion.rotation.z = -Math.PI / 2;
    root.add(ion);
    var count = 1500;
    seeds = new Float32Array(count * 3);
    for (var p = 0; p < count; p++) {
      seeds[p * 3] = random();
      seeds[p * 3 + 1] = random() * Math.PI * 2;
      seeds[p * 3 + 2] = Math.sqrt(random());
    }
    var g = new THREE.BufferGeometry();
    g.setAttribute("position", new THREE.BufferAttribute(new Float32Array(count * 3), 3));
    particles = new THREE.Points(g, new THREE.PointsMaterial({
      color: 0xf0e6d0, size: 0.09, transparent: true, opacity: 0.55,
      blending: THREE.AdditiveBlending, depthWrite: false,
    }));
    root.add(particles);
  }
  root.rotation.y = 0.5;
  return {
    object: root,
    cameraDistance: 13,
    update: function (t) {
      nucleus.rotation.set(t * 0.1, t * 0.23, 0);
      if (!particles) {
        return;
      }
      // Dust drifts back along +x, spreading and curving as it goes.
      var arr = particles.geometry.attributes.position.array;
      for (var p = 0; p < arr.length / 3; p++) {
        var u = (seeds[p * 3] + t * 0.05) % 1;
        var spread = 0.3 + u * 2.2;
        var r = seeds[p * 3 + 2] * spread;
        arr[p * 3] = u * TAIL;
        arr[p * 3 + 1] = Math.cos(seeds[p * 3 + 1]) * r + u * u * 2.5;
        arr[p * 3 + 2] = Math.sin(seeds[p * 3 + 1]) * r;
      }
      particles.geometry.attributes.position.needsUpdate = true;
    },
  };
}

var BUILDERS = {
  neutron_star: neutronStar,
  black_hole: blackHole,
  quasar: quasar,
  rogue_planet: roguePlanet,
  interstellar_comet: interstellarComet,
};

function start(view) {
  var builder = BUILDERS[view && view.kind];
  if (!builder || !canvas) {
    return;
  }
  var renderer;
  try {
    renderer = new THREE.WebGLRenderer({ canvas: canvas, antialias: true });
  } catch (err) {
    return; // No WebGL: the SVG still stays.
  }
  var random = makeRandom(20251001);
  var scene = new THREE.Scene();
  scene.background = new THREE.Color(0x05070c);
  scene.add(starfield(random));
  var built = builder(view, random);
  var pivot = new THREE.Group();
  pivot.add(built.object);
  scene.add(pivot);
  var camera = new THREE.PerspectiveCamera(40, 1, 0.1, 1000);
  camera.position.set(0, 0, built.cameraDistance);

  canvas.hidden = false;
  viewport.classList.add("phenomrender-live");

  function resize() {
    var w = viewport.clientWidth;
    var h = viewport.clientHeight;
    if (!w || !h) {
      return;
    }
    renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
    renderer.setSize(w, h, false);
    camera.aspect = w / h;
    camera.updateProjectionMatrix();
  }

  var timer = new THREE.Timer();  // the old clock class is deprecated since r183 (MAP.144)
  var STILL_TIME = 3.7;
  function draw() {
    timer.update();
    var t = reducedMotion ? STILL_TIME : timer.getElapsed();
    scene.updateMatrixWorld();
    built.update(t, camera);
    renderer.render(scene, camera);
  }

  // Drag turns the object; it never scrolls or zooms the page.
  var dragging = null;
  canvas.addEventListener("pointerdown", function (event) {
    dragging = { x: event.clientX, y: event.clientY };
    canvas.setPointerCapture(event.pointerId);
  });
  canvas.addEventListener("pointermove", function (event) {
    if (!dragging) {
      return;
    }
    pivot.rotation.y += (event.clientX - dragging.x) * 0.01;
    pivot.rotation.x = Math.max(-1.3, Math.min(1.3, pivot.rotation.x + (event.clientY - dragging.y) * 0.01));
    dragging = { x: event.clientX, y: event.clientY };
    if (reducedMotion) {
      draw();
    }
  });
  function endDrag() {
    dragging = null;
  }
  canvas.addEventListener("pointerup", endDrag);
  canvas.addEventListener("pointercancel", endDrag);

  if (typeof ResizeObserver === "function") {
    new ResizeObserver(function () {
      resize();
      draw();
    }).observe(viewport);
  }
  resize();
  draw();
  if (reducedMotion) {
    return;
  }

  // Only animate while the view is on screen.
  var visible = true;
  if (typeof IntersectionObserver === "function") {
    new IntersectionObserver(function (entries) {
      visible = entries[0].isIntersecting;
    }).observe(viewport);
  }
  renderer.setAnimationLoop(function () {
    if (visible && !document.hidden) {
      draw();
    }
  });
}

start(readView());
