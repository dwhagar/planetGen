// html/static/bodyview.js
//
// The rendered example on an object page (planetgen/web/lib/bodyview.py): the look record in a canvas's `data-look`
// painted by the system map's own painters (`solidpaint.js`, `starpaint.js`) and shaded as a lit sphere, with the
// atmosphere's limb glow or a star's corona drawn around it. Plain 2D canvas, so it needs no WebGL. A giant, which
// the painters do not cover yet, is its tinted disc. With no script the figure stays hidden and the property table
// says everything the picture does.

const VERSION_QUERY = new URL(import.meta.url).search;
const { parseLook, parseHex } = await import(`./bodytexture.js${VERSION_QUERY}`);
const { paintSolid } = await import(`./solidpaint.js${VERSION_QUERY}`);
const { glowFor, paintStar } = await import(`./starpaint.js${VERSION_QUERY}`);

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

function halo(context, centre, radius, outer, color, strength) {
  const gradient = context.createRadialGradient(centre, centre, radius * 0.96, centre, centre, outer);
  gradient.addColorStop(0, rgba(color, strength));
  gradient.addColorStop(1, rgba(color, 0));
  context.fillStyle = gradient;
  context.beginPath();
  context.arc(centre, centre, outer, 0, Math.PI * 2);
  context.fill();
}

function draw(canvas, look) {
  const size = canvas.width;
  const centre = size / 2;
  const star = look.body === "star";
  const halfWidth = star ? glowFor(look.lumSolar).scale : look.glow ? 1.12 : 1.04;
  const radius = size / 2 / halfWidth;
  const context = canvas.getContext("2d");
  context.clearRect(0, 0, size, size);
  const rim = star ? look.color : look.glow;
  if (rim) halo(context, centre, radius, size / 2, rim, star ? glowFor(look.lumSolar).strength * 0.8 : 0.55);

  const image = context.createImageData(size, size);
  const pixels = texture(look);
  const base = parseHex(look.color);
  const light = Math.hypot(LIGHT[0], LIGHT[1], LIGHT[2]);
  for (let y = 0; y < size; y++) {
    for (let x = 0; x < size; x++) {
      const nx = (x + 0.5 - centre) / radius;
      const ny = (y + 0.5 - centre) / radius;
      const r2 = nx * nx + ny * ny;
      if (r2 >= 1) continue;
      const nz = Math.sqrt(1 - r2);
      let colour = base;
      if (pixels) {
        // Orthographic sphere to the texture's longitude and latitude, the sphere's north up and its middle at 0 degrees.
        const lon = Math.atan2(nx, nz);
        const lat = Math.asin(-ny);
        const u = (lon / (2 * Math.PI) + 0.5) % 1;
        const v = 0.5 - lat / Math.PI;
        const tx = Math.min(TEXTURE[0] - 1, Math.floor(u * TEXTURE[0]));
        const ty = Math.min(TEXTURE[1] - 1, Math.floor(v * TEXTURE[1]));
        const at = (ty * TEXTURE[0] + tx) * 4;
        colour = [pixels[at], pixels[at + 1], pixels[at + 2]];
      }
      // A star is its own light; a world is lit from the upper left, with a little ambient on the night side.
      const shade = star ? 0.75 + 0.25 * nz : 0.12 + 0.88 * Math.max(0, (nx * LIGHT[0] + ny * LIGHT[1] + nz * LIGHT[2]) / light);
      const edge = Math.min(1, (1 - Math.sqrt(r2)) * radius);
      const at = (y * size + x) * 4;
      image.data[at] = Math.min(255, colour[0] * shade);
      image.data[at + 1] = Math.min(255, colour[1] * shade);
      image.data[at + 2] = Math.min(255, colour[2] * shade);
      image.data[at + 3] = Math.round(255 * edge);
    }
  }
  // The sphere over the halo: draw it on a scratch canvas so its soft edge blends with what is behind it.
  const scratch = document.createElement("canvas");
  scratch.width = scratch.height = size;
  scratch.getContext("2d").putImageData(image, 0, 0);
  context.drawImage(scratch, 0, 0);
}

for (const canvas of document.querySelectorAll("canvas[data-look]")) {
  const look = parseLook(canvas.dataset.look);
  if (!look) continue;
  draw(canvas, look);
  const figure = canvas.closest("figure");
  if (figure) figure.hidden = false;
  canvas.dataset.drawn = "1";
}
