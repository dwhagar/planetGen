// html/static/bodyview.js
//
// The rendered example on an object page (planetgen/web/lib/bodyview.py): the look record in a canvas's `data-look`
// painted by the system map's own painters (`solidpaint.js`, `starpaint.js`) and shaded as a lit sphere, with the
// atmosphere's limb glow or a star's corona drawn around it. The body turns slowly on its axis (MAP.328, spin.js). Plain
// 2D canvas, so it needs no WebGL. A giant, which the painters do not cover yet, is its tinted disc. With no script the
// figure stays hidden and the property table says everything the picture does.

const VERSION_QUERY = new URL(import.meta.url).search;
const { parseLook, parseHex } = await import(`./bodytexture.js${VERSION_QUERY}`);
const { paintSolid } = await import(`./solidpaint.js${VERSION_QUERY}`);
const { glowFor, paintStar } = await import(`./starpaint.js${VERSION_QUERY}`);
const Spin = await import(`./spin.js${VERSION_QUERY}`);

const TEXTURE = [512, 256];
const LIGHT = [-0.45, -0.4, 0.8];

function rgba(color, alpha) {
  const c = parseHex(color);
  return "rgba(" + c[0] + "," + c[1] + "," + c[2] + "," + alpha + ")";
}

function texture(look) {
  if (look.body === "star") return paintStar(look, TEXTURE[0], TEXTURE[1]);
  if (look.body === "solid") return paintSolid(look, TEXTURE[0], TEXTURE[1]);
  return null;
}

function drawHalo(context, centre, radius, outer, color, strength) {
  const gradient = context.createRadialGradient(centre, centre, radius * 0.96, centre, centre, outer);
  gradient.addColorStop(0, rgba(color, strength));
  gradient.addColorStop(1, rgba(color, 0));
  context.fillStyle = gradient;
  context.beginPath();
  context.arc(centre, centre, outer, 0, Math.PI * 2);
  context.fill();
}

// Everything about a picture that does not change as the body turns, worked out once: for each pixel of the disc where
// it is in the canvas, its longitude in the body's own frame, the texture row it samples, its shade (a world is lit from
// the upper left whatever it does, so the terminator stays put while the surface turns) and how much of the soft edge it
// takes. A frame is then one pass of lookups.
function prepare(canvas, look) {
  const size = canvas.width;
  const centre = size / 2;
  const star = look.body === "star";
  const halfWidth = star ? glowFor(look.lumSolar).scale : look.glow ? 1.12 : 1.04;
  const radius = size / 2 / halfWidth;
  const pixels = texture(look);
  const base = parseHex(look.color);
  const light = Math.hypot(LIGHT[0], LIGHT[1], LIGHT[2]);
  const places = [];
  for (let y = 0; y < size; y++) {
    for (let x = 0; x < size; x++) {
      const nx = (x + 0.5 - centre) / radius;
      const ny = (y + 0.5 - centre) / radius;
      const r2 = nx * nx + ny * ny;
      if (r2 >= 1) continue;
      const nz = Math.sqrt(1 - r2);
      // A star is its own light; a world is lit from the upper left, with a little ambient on the night side.
      const shade = star ? 0.75 + 0.25 * nz : 0.12 + 0.88 * Math.max(0, (nx * LIGHT[0] + ny * LIGHT[1] + nz * LIGHT[2]) / light);
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
  if (rim) drawHalo(halo.getContext("2d"), centre, radius, size / 2, rim, star ? glowFor(look.lumSolar).strength * 0.8 : 0.55);
  // The sphere is drawn on a scratch canvas so its soft edge blends with what is behind it.
  const scratch = document.createElement("canvas");
  scratch.width = scratch.height = size;
  const scratchContext = scratch.getContext("2d");
  const image = scratchContext.createImageData(size, size);
  const context = canvas.getContext("2d");

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
  };
}

// MAP.328: a body with a stored rotation period turns slowly (spin.js: one law for every view), a frame every 33 ms while
// it is on screen and Spin is on; with Spin off, or reduced motion asked for, or no stored period, it is drawn once and
// stays still.
const FRAME_MS = 33;

function animate(canvas, look, frame) {
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
  if (!button || !(Spin.rateOf(look) > 0)) return;
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
