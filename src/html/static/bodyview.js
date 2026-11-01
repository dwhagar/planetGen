// html/static/bodyview.js
//
// The rendered example on an object page (planetgen/web/lib/bodyview.py): the look record in a canvas's `data-look`
// painted by the system map's own painters (`solidpaint.js`, `starpaint.js`) and shaded as a lit sphere, with the
// atmosphere's limb glow or a star's corona drawn around it. The body turns slowly on its axis (MAP.328, spin.js). A world is
// seen from the side of its orbit (UX.244): the star to the left and a little in front, so the day and night line is a curve with more than half the disc lit, the axial tilt
// leaning the poles and the line's ends, and the surface turning past it at the display speed; a tidally locked world holds
// still with its day side toward the star. Plain
// 2D canvas, so it needs no WebGL. A giant, which the painters do not cover yet, is its tinted disc. With no script the
// figure stays hidden and the property table says everything the picture does.

const VERSION_QUERY = new URL(import.meta.url).search;
const { parseLook, parseHex } = await import(`./bodytexture.js${VERSION_QUERY}`);
const { paintSolid } = await import(`./solidpaint.js${VERSION_QUERY}`);
const { paintGiant } = await import(`./giantpaint.js${VERSION_QUERY}`);
const { glowFor, paintStar } = await import(`./starpaint.js${VERSION_QUERY}`);
const Spin = await import(`./spin.js${VERSION_QUERY}`);
const { isLumpy } = await import(`./lumpymesh.js${VERSION_QUERY}`);
const { lumpMesh, rasterise } = await import(`./lumpyview.js${VERSION_QUERY}`);

const TEXTURE = [512, 256];
// UX.250: a world is lit from the left and a little from the front: the star is PHASE_DEG round from the line of sight toward the left, so the lit part is
// more than half the disc ((1 + cos phase) / 2, 79 percent) and the day and night line is a curve bowing into the night side, not a straight line down the middle.
const PHASE_DEG = 55;
const LIGHT = [-Math.sin((PHASE_DEG * Math.PI) / 180), 0, Math.cos((PHASE_DEG * Math.PI) / 180)];

function rgba(color, alpha) {
  const c = parseHex(color);
  return "rgba(" + c[0] + "," + c[1] + "," + c[2] + "," + alpha + ")";
}

function texture(look) {
  if (look.body === "star") return paintStar(look, TEXTURE[0], TEXTURE[1]);
  if (look.body === "solid") return paintSolid(look, TEXTURE[0], TEXTURE[1]);
  if (look.body === "gas") return paintGiant(look, TEXTURE[0], TEXTURE[1]);
  return null;
}

function drawHalo(context, centreX, centreY, radius, outer, color, strength) {
  const gradient = context.createRadialGradient(centreX, centreY, radius * 0.96, centreX, centreY, outer);
  gradient.addColorStop(0, rgba(color, strength));
  gradient.addColorStop(1, rgba(color, 0));
  context.fillStyle = gradient;
  context.beginPath();
  context.arc(centreX, centreY, outer, 0, Math.PI * 2);
  context.fill();
}

// UX.250: how lit a point is from the cosine of the light on it: 0 at the day and night line and full from TERMINATOR_WIDTH of a cosine past it, so the
// line is a clear curve and most of the lit side is bright, not a slope that fades across half the disc.
const TERMINATOR_WIDTH = 0.4;
// MAP.277: how far the glow of an air reaches past the limb, in radii, when the look record does not say (the composition of the air is not stored).
const DEFAULT_GLOW_WIDTH = 0.08;
// MAP.277: the air carries the star's light past the line: `look.twilight` (0 for no air) moves the line that far onto the night side and softens it by the
// same amount, so a world with a thick air has a wide dusk and a world with none has a hard edge.
function lit(cosine, twilight) {
  const t = Math.min(1, Math.max(0, (cosine + twilight) / (TERMINATOR_WIDTH + twilight)));
  return t * t * (3 - 2 * t);
}

// UX.250: a tidally locked world keeps the left of its picture for the arrow to the star (and the word under it) and the foot for the label, so neither
// is drawn over the disc: the disc is fitted into what is left.
const LOCK_MARGIN_PX = 52;
const LOCK_FOOT_PX = 28;

// UX.245: a tidally locked world says so on the picture itself: an arrow on the left, the side the star (or, for a moon, its planet)
// is on, and "Tidally locked" under the disc.
function drawLockLabel(context, size, look, arrowY) {
  const what = look.kind === "moon" ? "planet" : "star";
  context.save();
  context.font = "600 13px sans-serif";
  context.textBaseline = "middle";
  const text = "Tidally locked";
  const width = context.measureText(text).width + 16;
  context.fillStyle = "rgba(10, 14, 26, 0.78)";
  context.fillRect((size - width) / 2, size - 26, width, 20);
  context.fillStyle = "#f2f5ff";
  context.textAlign = "center";
  context.fillText(text, size / 2, size - 16);
  const y = arrowY;
  // On a chip like the label's, so it reads on a light page too.
  context.fillStyle = "rgba(10, 14, 26, 0.78)";
  context.fillRect(1, y - 26, 50, 40);
  context.strokeStyle = "rgba(242, 245, 255, 0.95)";
  context.fillStyle = "rgba(242, 245, 255, 0.95)";
  context.lineWidth = 2;
  context.beginPath();
  context.moveTo(30, y);
  context.lineTo(8, y);
  context.stroke();
  context.beginPath();
  context.moveTo(4, y);
  context.lineTo(14, y - 5);
  context.lineTo(14, y + 5);
  context.closePath();
  context.fill();
  context.textAlign = "left";
  context.fillText(what, 4, y - 14);
  context.restore();
}

// Everything about a picture that does not change as the body turns, worked out once: for each pixel of the disc where
// it is in the canvas, its longitude in the body's own frame, the texture row it samples, its shade (a world is lit from
// the left whatever it does, so the terminator stays put while the surface turns) and how much of the soft edge it
// takes. A frame is then one pass of lookups.
function prepare(canvas, look) {
  const size = canvas.width;
  const star = look.body === "star";
  const halfWidth = star ? glowFor(look.lumSolar).scale : 1.04 + (look.glow ? look.glowWidth || DEFAULT_GLOW_WIDTH : 0);
  const pixels = texture(look);
  const locked = Spin.isLocked(look) && !star;
  // The square the disc and its halo fill: the whole canvas, or for a locked world what the arrow and the label leave.
  const margin = locked ? Math.round(LOCK_MARGIN_PX * size / 320) : 0;
  const foot = locked ? Math.round(LOCK_FOOT_PX * size / 320) : 0;
  const box = Math.min(size - margin, size - foot);
  const centreX = size - box / 2;
  const centreY = (size - foot) / 2;
  const radius = box / 2 / halfWidth;
  canvas.dataset.disc = JSON.stringify({ x: centreX, y: centreY, radius: radius });
  const base = parseHex(look.color);
  const twilight = look.twilight || 0;
  const light = Math.hypot(LIGHT[0], LIGHT[1], LIGHT[2]);
  const lumpy = isLumpy(look);       // MAP.363: a body its weight has not made round is the lump the System Map draws, not a disc
  const places = [];
  for (let y = 0; y < size && !lumpy; y++) {
    for (let x = 0; x < size; x++) {
      const nx = (x + 0.5 - centreX) / radius;
      const ny = (y + 0.5 - centreY) / radius;
      const r2 = nx * nx + ny * ny;
      if (r2 >= 1) continue;
      const nz = Math.sqrt(1 - r2);
      // A star is its own light; a world is lit from the left, with a little ambient on the night side.
      const shade = star ? 0.75 + 0.25 * nz : 0.12 + 0.88 * lit((nx * LIGHT[0] + ny * LIGHT[1] + nz * LIGHT[2]) / light, twilight);
      // Orthographic sphere to the texture's longitude and latitude in the body's own frame (north up, leaned by its tilt).
      const [lon, lat] = Spin.bodyPoint(look, nx, -ny, nz);
      places.push({
        at: (y * size + x) * 4, lon: lon, shade: shade, alpha: Math.round(255 * Math.min(1, (1 - Math.sqrt(r2)) * radius)),
        row: Math.min(TEXTURE[1] - 1, Math.floor((0.5 - lat / Math.PI) * TEXTURE[1])) * TEXTURE[0],
      });
    }
  }
  const lons = Float32Array.from(places, (p) => p.lon);
  const rows = Int32Array.from(places, (p) => p.row);
  const shades = Float32Array.from(places, (p) => p.shade);
  const offsets = Int32Array.from(places, (p) => p.at);
  const alphas = Uint8Array.from(places, (p) => p.alpha);

  // The halo is behind the sphere and does not turn.
  const halo = document.createElement("canvas");
  halo.width = halo.height = size;
  const rim = star ? look.color : look.glow;
  if (rim) drawHalo(halo.getContext("2d"), centreX, centreY, radius, box / 2, rim, star ? glowFor(look.lumSolar).strength * 0.8 : look.glowStrength);
  // The sphere is drawn on a scratch canvas so its soft edge blends with what is behind it.
  const scratch = document.createElement("canvas");
  scratch.width = scratch.height = size;
  const scratchContext = scratch.getContext("2d");
  const image = scratchContext.createImageData(size, size);
  const context = canvas.getContext("2d");

  if (lumpy) {
    // The lump is turned with the body, so the silhouette and the light on it change every frame; the picture is filled from the grid of points of `lumpyview.js`.
    const mesh = lumpMesh(look);
    const scale = radius / mesh.extent;
    const covered = new Uint8Array(size * size);
    const depth = new Float32Array(size * size);
    const view = { angle: 0, tiltDeg: look.tiltDeg };
    const lightLength = Math.hypot(LIGHT[0], LIGHT[1], LIGHT[2]);
    const data = image.data;
    const width = TEXTURE[0];
    const paint = (pixel, longitude, latitude, nx, ny, nz) => {
      let r = base[0];
      let g = base[1];
      let b = base[2];
      if (pixels) {
        let u = longitude / (2 * Math.PI) + 0.5;
        u -= Math.floor(u);
        const t = (Math.min(TEXTURE[1] - 1, Math.floor((0.5 - latitude / Math.PI) * TEXTURE[1])) * width + Math.min(width - 1, Math.floor(u * width))) * 4;
        r = pixels[t];
        g = pixels[t + 1];
        b = pixels[t + 2];
      }
      const shade = 0.12 + 0.88 * lit((nx * LIGHT[0] + ny * LIGHT[1] + nz * LIGHT[2]) / lightLength, twilight);
      const at = pixel * 4;
      data[at] = Math.min(255, r * shade);
      data[at + 1] = Math.min(255, g * shade);
      data[at + 2] = Math.min(255, b * shade);
      data[at + 3] = 255;
    };
    return function lumpFrame(angle) {
      view.angle = angle;
      data.fill(0);
      rasterise(mesh, view, size, centreX, centreY, scale, covered, depth, paint);
      scratchContext.putImageData(image, 0, 0);
      context.clearRect(0, 0, size, size);
      context.drawImage(halo, 0, 0);
      context.drawImage(scratch, 0, 0);
      if (locked) drawLockLabel(context, size, look, centreY);
    };
  }

  // Draws the body turned `angle` radians (see `Spin.angleAt`).
  return function frame(angle) {
    const data = image.data;
    const width = TEXTURE[0];
    for (let i = 0; i < offsets.length; i++) {
      let r = base[0];
      let g = base[1];
      let b = base[2];
      if (pixels) {
        // The longitude the picture shows here is the body's, less how far it has turned.
        let u = (lons[i] - angle) / (2 * Math.PI) + 0.5;
        u -= Math.floor(u);
        const t = (rows[i] + Math.min(width - 1, Math.floor(u * width))) * 4;
        r = pixels[t];
        g = pixels[t + 1];
        b = pixels[t + 2];
      }
      const shade = shades[i];
      const at = offsets[i];
      data[at] = Math.min(255, r * shade);
      data[at + 1] = Math.min(255, g * shade);
      data[at + 2] = Math.min(255, b * shade);
      data[at + 3] = alphas[i];
    }
    scratchContext.putImageData(image, 0, 0);
    context.clearRect(0, 0, size, size);
    context.drawImage(halo, 0, 0);
    context.drawImage(scratch, 0, 0);
    if (locked) drawLockLabel(context, size, look, centreY);
  };
}

// MAP.328: a body with a stored rotation period turns slowly (spin.js: one law for every view), a frame every 33 ms while
// it is on screen and Spin is on; with Spin off, or reduced motion asked for, or no stored period, it is drawn once and
// stays still.
const FRAME_MS = 33;

function animate(canvas, look, frame) {
  if (Spin.isLocked(look) && look.body !== "star") {
    // A locked world keeps its day side toward the star, on the left: the point of it facing the star is where the light comes from (to the left and front), so it holds still.
    const held = Spin.bodyPoint(look, LIGHT[0], LIGHT[1], LIGHT[2])[0];
    frame(held);
    canvas.dataset.angle = String(held);
    canvas.dataset.locked = "1";
    return;
  }
  let drawnAt = null;
  let last = -Infinity;
  let visible = true;
  if (typeof IntersectionObserver === "function") {
    new IntersectionObserver((entries) => { visible = entries[entries.length - 1].isIntersecting; }).observe(canvas);
  }
  function step(now) {
    requestAnimationFrame(step);
    if (document.hidden || !visible || now - last < FRAME_MS) return;
    last = now;
    const seconds = Spin.tick(now);
    if (seconds === drawnAt) return;
    drawnAt = seconds;
    frame(Spin.angleAt(look, seconds));
    canvas.dataset.angle = String(Spin.angleAt(look, seconds));
  }
  frame(Spin.angleAt(look, Spin.seconds()));
  canvas.dataset.angle = String(Spin.angleAt(look, Spin.seconds()));
  drawnAt = Spin.seconds();
  if (Spin.rateOf(look) > 0) requestAnimationFrame(step);
}

// The Spin switch under the picture (the same standing preference as the System Map's Menu, spin.js), for a body that has
// a period to turn by.
function spinButton(figure, look) {
  const button = figure.querySelector("[data-spin]");
  if (!button || !(Spin.rateOf(look) > 0) || (Spin.isLocked(look) && look.body !== "star")) return;
  function show() {
    const reduced = Spin.reducedMotion();
    button.setAttribute("aria-pressed", String(Spin.active()));
    button.disabled = reduced;
    button.title = reduced ? "Your device asks for reduced motion, so the picture holds still" : "Let the body turn slowly on its axis, or hold it still";
  }
  button.addEventListener("click", () => { Spin.setSwitch(!Spin.switchedOn()); show(); });
  show();
  button.hidden = false;
}

for (const canvas of document.querySelectorAll("canvas[data-look]")) {
  const look = parseLook(canvas.dataset.look);
  if (!look) continue;
  animate(canvas, look, prepare(canvas, look));
  const figure = canvas.closest("figure");
  if (figure) {
    figure.hidden = false;
    spinButton(figure, look);
  }
  canvas.dataset.drawn = "1";
}
