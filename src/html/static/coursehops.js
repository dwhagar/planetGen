// static/coursehops.js -- the Galaxy Map's NAV course: which hops are jumps through unknown space and how they are drawn (NAV.36).
//
// The route's points carry `hopUnknown` (the hop that arrives at that point is a jump through unknown space, NAV.12). The map draws
// each run of known hops as one solid line and each unknown hop on its own: red, dashed and glowing. A line cannot glow in WebGL, so
// the glow is a row of soft red sprites laid along the hop, as many as the hop is long on the screen; it pulses only where the visitor
// has not asked for reduced motion.

const MAX_BEADS = 240;
const BEAD_SPACING_PX = 6;
export const GLOW_PX = 22;

// The runs of a course: `[{from, to, unknown}]`, `from` and `to` the indexes of its first and last point; known hops that follow each
// other make one run, each unknown hop is a run of its own.
export function courseRuns(points) {
  const runs = [];
  for (let n = 1; n < points.length; n += 1) {
    const unknown = !!points[n].hopUnknown;
    const last = runs[runs.length - 1];
    if (!unknown && last && !last.unknown) {
      last.to = n;
    } else {
      runs.push({ from: n - 1, to: n, unknown: unknown });
    }
  }
  return runs;
}

// How many glow sprites a hop of `px` screen pixels takes.
export function beadCount(px) {
  if (!isFinite(px)) return MAX_BEADS;
  return Math.max(2, Math.min(MAX_BEADS, Math.ceil(px / BEAD_SPACING_PX)));
}

// The glow's brightness at time `now` (ms): steady under reduced motion, a slow pulse otherwise.
export function pulseFactor(now, reducedMotion) {
  if (reducedMotion) return 1;
  return 0.65 + 0.35 * Math.sin((now / 2400) * 2 * Math.PI);
}

// A soft round glow in `color`, as a sprite texture. (`THREE` is passed in: this file loads nothing, so it also runs under node.)
export function makeGlowTexture(THREE, color) {
  const size = 64;
  const canvasEl = document.createElement("canvas");
  canvasEl.width = canvasEl.height = size;
  const ctx = canvasEl.getContext("2d");
  const r = size / 2;
  const gradient = ctx.createRadialGradient(r, r, 0, r, r, r);
  gradient.addColorStop(0, color);
  gradient.addColorStop(0.35, color);
  gradient.addColorStop(1, "rgba(0,0,0,0)");
  ctx.fillStyle = gradient;
  ctx.fillRect(0, 0, size, size);
  const texture = new THREE.CanvasTexture(canvasEl);
  texture.colorSpace = THREE.SRGBColorSpace;
  return texture;
}

// The glow of one unknown hop from `a` to `b` (Vector3s): `group` goes into the scene; `update` lays the sprites along the hop for the
// camera and the size of a screen pixel in the sprites' units (`per`).
export function createHopGlow(THREE, a, b, texture) {
  const group = new THREE.Group();
  const sprites = [];
  const ends = [new THREE.Vector3(), new THREE.Vector3()];
  for (let n = 0; n < MAX_BEADS; n += 1) {
    const sprite = new THREE.Sprite(new THREE.SpriteMaterial({
      map: texture, transparent: true, depthWrite: false, depthTest: false, sizeAttenuation: false, opacity: 0.2,
    }));
    sprite.visible = false;
    group.add(sprite);
    sprites.push(sprite);
  }
  function update(camera, canvasHeight, now, reducedMotion) {
    const per = (2 * Math.tan(THREE.MathUtils.degToRad(camera.fov) / 2)) / (canvasHeight || 1);
    ends[0].copy(a).project(camera);
    ends[1].copy(b).project(camera);
    const behind = ends[0].z > 1 || ends[1].z > 1;
    const px = behind ? Infinity : Math.hypot((ends[1].x - ends[0].x) * (canvasHeight || 1) / 2, (ends[1].y - ends[0].y) * (canvasHeight || 1) / 2);
    const count = beadCount(px);
    const pulse = pulseFactor(now, reducedMotion);
    sprites.forEach((sprite, n) => {
      sprite.visible = n < count;
      if (!sprite.visible) return;
      sprite.position.lerpVectors(a, b, count === 1 ? 0.5 : n / (count - 1));
      sprite.scale.set(GLOW_PX * per, GLOW_PX * per, 1);
      sprite.material.opacity = 0.2 * pulse;
    });
  }
  return { group, update, sprites };
}
