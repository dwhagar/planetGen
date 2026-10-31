// html/static/lumpymesh.js
//
// The lumpy mesh of a small body (MAP.279): a moon or an asteroid that its own weight has not pulled round is drawn as an
// elongated, bumpy lump instead of a sphere. Whether a body is round is the roundness ratio the generator stores
// (`roundness`, GEN.1158 and GEN.1159: its radius over the radius at which its weight beats its strength; 1 or more is
// round), carried on the look record (planetgen/web/maps/bodylook.py). The shape is display-only noise from the look's `seed`,
// a pure function of the record, so the same body is the same lump on every page and every machine.
//
// What the stored ratio sets: how far from round the body is, `1 - roundness` from 0 (just short of round) to 1 (far from it).
// The farther, the longer the long axis over the short one and the taller the bumps. A body whose record has no ratio (a system
// stored before GEN.1159) is a lump only under `LUMPY_MAX_KM`, Boss's line of 10:26 PM EDT 2026-10-23, with `FALLBACK_ROUNDNESS`.
// The mean radius of the mesh is 1 (the volume of the unit sphere), so the marker's own radius still sets the size on screen.
//
// Stand-ins (docs/design/research/gap-filling-plan.md section 6.160): the numbers below only shape a picture.

const VERSION_QUERY = new URL(import.meta.url).search;
const { makeRandom } = await import(`./bodytexture.js${VERSION_QUERY}`);

export const LUMPY_MAX_KM = 400;
/** The radius, km, under which a body whose roundness is not stored is drawn as a lump. [Stand-in LUM1, tier T9, tweak L0.] */

export const FALLBACK_ROUNDNESS = 0.5;
/** The roundness ratio a lump without a stored one is drawn with. [Stand-in LUM2, tier T9, tweak L4.] */

export const MAX_ELONGATION = 0.8;
/** The long axis over the short one is `1 + MAX_ELONGATION` at most, for a body far from round (the longest of the small Solar System bodies are near two to one). [Stand-in LUM3, tier T6, tweak L0.] */

export const BUMP_BASE = 0.05;
export const BUMP_PER_DEFICIT = 0.2;
/** The bumps' height as a share of the radius: `BUMP_BASE` for a body just short of round plus `BUMP_PER_DEFICIT` times how far from round it is. [Stand-in LUM4, tier T9, tweak L3.] */

export const BUMP_WAVES = 7;
/** How many smooth waves are added to make the bumps. [Stand-in LUM5, tier T9, tweak L0.] */

const SEED_SALT = 0x4c554d50;

/** Whether the body is drawn as a lump: a solid body with a stored roundness under 1, or without one and under `LUMPY_MAX_KM`. */
export function isLumpy(look) {
  if (!look || look.body !== "solid") return false;
  if (typeof look.roundness === "number") return look.roundness < 1;
  return typeof look.radiusKm === "number" && look.radiusKm < LUMPY_MAX_KM;
}

/**
 * The shape of the lump of `look`: `{ axes: [a, b, c], waves: [{ dir, k, phase, weight }], amplitude, extent }`. The axes are the
 * long, middle and short scale with a geometric mean of 1, `amplitude` the bump height, and `extent` the farthest the surface
 * reaches from the centre (the picture is framed to it).
 */
export function shapeOf(look) {
  const ratio = typeof look.roundness === "number" ? look.roundness : FALLBACK_ROUNDNESS;
  const deficit = Math.min(Math.max(1 - ratio, 0), 1);
  const random = makeRandom(((look.seed >>> 0) ^ SEED_SALT) >>> 0);
  const long = 1 + MAX_ELONGATION * deficit * (0.5 + 0.5 * random());
  const middle = 1 + (long - 1) * (0.25 + 0.75 * random());
  const mean = Math.cbrt(long * middle);
  const axes = [long / mean, middle / mean, 1 / mean];
  const waves = [];
  let total = 0;
  for (let i = 0; i < BUMP_WAVES; i++) {
    const z = 2 * random() - 1;
    const turn = 2 * Math.PI * random();
    const flat = Math.sqrt(1 - z * z);
    const weight = 1 / (i + 1);
    total += weight;
    waves.push({ dir: [flat * Math.cos(turn), flat * Math.sin(turn), z], k: 2 * Math.PI * (1.2 + 2.2 * random()), phase: 2 * Math.PI * random(), weight: weight });
  }
  for (const wave of waves) wave.weight /= total;
  const amplitude = BUMP_BASE + BUMP_PER_DEFICIT * deficit;
  return { axes: axes, waves: waves, amplitude: amplitude, extent: Math.max(...axes) * (1 + amplitude) };
}

/** How far the surface is from the centre along the unit direction (x, y, z), as a multiple of the mean radius, before the axes stretch it: `1 +- amplitude`. */
export function bumpAt(shape, x, y, z) {
  let sum = 0;
  for (const wave of shape.waves) sum += wave.weight * Math.sin(wave.k * (wave.dir[0] * x + wave.dir[1] * y + wave.dir[2] * z) + wave.phase);
  return 1 + shape.amplitude * sum;
}

/** The point of the surface along the unit direction (x, y, z): `[px, py, pz]`. */
export function surfacePoint(shape, x, y, z) {
  const bump = bumpAt(shape, x, y, z);
  return [x * bump * shape.axes[0], y * bump * shape.axes[1], z * bump * shape.axes[2]];
}

/**
 * The lump of `look` as a three.js geometry on the UV sphere the system map's painted textures are mapped to (so the surface
 * paint lies on it as on a sphere), and the extent the picture is framed to: `{ geometry, extent }`.
 */
export function lumpyGeometry(THREE, look, widthSegments = 48, heightSegments = 32) {
  const shape = shapeOf(look);
  const geometry = new THREE.SphereGeometry(1, widthSegments, heightSegments);
  const position = geometry.attributes.position;
  for (let i = 0; i < position.count; i++) {
    const [px, py, pz] = surfacePoint(shape, position.getX(i), position.getY(i), position.getZ(i));
    position.setXYZ(i, px, py, pz);
  }
  position.needsUpdate = true;
  geometry.computeVertexNormals();
  geometry.computeBoundingSphere();
  return { geometry: geometry, extent: shape.extent };
}
