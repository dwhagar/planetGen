// html/static/maplogruler.js
//
// MAP.336: the Galaxy Map's vertical log ruler, the pure half: the named scales along it, where a distance stands on it and what the
// distance readout says. The ruler is one logarithmic scale for the camera's distance from the whole galaxy (the top) down to the
// nearest the map goes (the bottom, a moon framed); a stop is a named scale and a Go to for it, centred on the focus. The DOM is
// mapnavpad.js's. No DOM and no three.js: tests/js/maplogruler.test.mjs holds it to this description.

// The named scales, nearest last: {key, label, pc}. The galaxy stop is the farthest the camera goes, whatever the galaxy's size
// (`galaxyFarPc`), so it stands at the top; a stop nearer than the floor, or farther than the galaxy, is left off.
export const STOPS = [
  { key: "galaxy", label: "Galaxy", pc: null },
  { key: "arm", label: "Arm and halo", pc: 1000 },
  { key: "region", label: "Region", pc: 40 },
  { key: "sector", label: "Sector", pc: 4 },
  { key: "neighbourhood", label: "Neighbourhood", pc: 1 },
  { key: "system", label: "Planetary system", pc: 1e-3 },
  { key: "planet", label: "Planet", pc: 1e-8 },
  { key: "moon", label: "Moon", pc: 1e-10 },
  { key: "surface", label: "Surface", pc: 1e-12 },
];

// The stops that fit between `farPc` (the whole galaxy in view) and `floorPc`, with each one's distance and share along the ruler
// (0 the top, 1 the bottom).
export function stopsBetween(farPc, floorPc) {
  const span = Math.log(farPc) - Math.log(floorPc);
  return STOPS.map(function (stop) {
    const pc = stop.pc === null ? farPc : stop.pc;
    return { key: stop.key, label: stop.label, pc: pc, share: (Math.log(farPc) - Math.log(pc)) / span };
  }).filter(function (stop) {
    return stop.pc <= farPc * (1 + 1e-9) && stop.pc >= floorPc * (1 - 1e-9) && stop.share >= 0 && stop.share <= 1 + 1e-9;
  });
}

// The distance at `share` along the ruler, and the share of `pc`.
export function distanceAt(share, farPc, floorPc) {
  const s = Math.max(0, Math.min(1, share));
  return Math.exp(Math.log(farPc) - s * (Math.log(farPc) - Math.log(floorPc)));
}

export function shareOf(pc, farPc, floorPc) {
  const share = (Math.log(farPc) - Math.log(pc)) / (Math.log(farPc) - Math.log(floorPc));
  return Math.max(0, Math.min(1, share));
}

// The named scale nearest `pc` on the log scale, so the readout can say where the camera is.
export function scaleNamed(pc, stops) {
  let named = stops[0];
  stops.forEach(function (stop) {
    if (Math.abs(Math.log(pc / stop.pc)) < Math.abs(Math.log(pc / named.pc))) named = stop;
  });
  return named;
}

// The readout: "Distance from focus: 2.4 pc (7.8 ly)", and with a named thing of the scene up "; nearest: Kanbee 3, 0.9 pc (2.9 ly)".
// `format(pc)` writes a distance in the reader's units, with its companion unit (unitladder's formatDistancePc); `nearest` is
// {name, distance} or null.
export function readoutText(distPc, nearest, format) {
  let text = "Distance from focus: " + format(distPc);
  if (nearest && nearest.name) text += "; nearest: " + nearest.name + ", " + format(nearest.distance);
  return text;
}
