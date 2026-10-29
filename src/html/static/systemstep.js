// html/static/systemstep.js
//
// Stepping from the open system to a neighbouring one (MAP.121): the same six directions the step buttons and keys use
// at the other levels (north, south, east, west, out from the galaxy's centre and in toward it), read in the frame of
// the galaxy at the system's own place, go to the nearest system that lies that way. The System Map has no grid to step
// through, so "that way" is a cone about the direction. Pure: runs under plain node (tests/js/systemstep.test.mjs).

// A display tweak with no source (L4 on the stand-in scale; a step is a way to move about, not a measurement): the
// half-angle, in degrees, of the cone about a direction within which a system counts as lying that way. At 60 degrees
// the six cones overlap a little, so a system lies in at least one and a step is rarely refused where there are neighbours.
export const CONE_HALF_ANGLE_DEG = 60;

// The unit vector of `direction` at the point (x, y, z) in parsecs, in the galaxy's frame (the origin at its centre, z
// toward the north pole): out and in along the radius, north and south along the colatitude (south is the way it
// grows), east and west along the longitude (east is the way it grows). Null for an unknown direction or a point on the
// axis through the poles or at the centre, where north, east and west have no meaning.
export function directionAt(point, direction) {
  const [x, y, z] = point;
  const plane = Math.hypot(x, y);
  const r = Math.hypot(plane, z);
  if (!(r > 0)) return null;
  const radial = [x / r, y / r, z / r];
  if (direction === "out") return radial;
  if (direction === "in") return radial.map((v) => -v);
  if (!(plane > 0)) return null;
  const east = [-y / plane, x / plane, 0];
  const south = [(z * x) / (plane * r), (z * y) / (plane * r), -plane / r];
  if (direction === "east") return east;
  if (direction === "west") return east.map((v) => -v);
  if (direction === "south") return south;
  if (direction === "north") return south.map((v) => -v);
  return null;
}

// The nearest of `candidates` ({x, y, z, ...}) from `from` ({x, y, z}) that lies within `coneDeg` degrees of
// `direction`, or null. `skip(candidate)` leaves the system we stand at (and anything else not to be offered) out.
export function nearestInDirection(from, candidates, direction, skip, coneDeg) {
  const axis = directionAt([from.x, from.y, from.z], direction);
  if (!axis) return null;
  const cosCone = Math.cos(((coneDeg === undefined ? CONE_HALF_ANGLE_DEG : coneDeg) * Math.PI) / 180);
  let best = null;
  let bestDistance = Infinity;
  candidates.forEach((candidate) => {
    if (skip && skip(candidate)) return;
    const dx = candidate.x - from.x;
    const dy = candidate.y - from.y;
    const dz = candidate.z - from.z;
    const distance = Math.hypot(dx, dy, dz);
    if (!(distance > 0)) return;
    if ((dx * axis[0] + dy * axis[1] + dz * axis[2]) / distance < cosCone) return;
    if (distance < bestDistance) {
      best = candidate;
      bestDistance = distance;
    }
  });
  return best;
}
