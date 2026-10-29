// html/static/regionlanding.js
//
// MAP.336: where "Go to" a region lands. Two modes, chosen by the reader's `map.regionLanding` setting (mappreferences.js; UX.236
// builds the page that sets it): "frame" (the default) puts the camera where the region's charted extent fills the view, "densest"
// goes on to the cluster of the region's heaviest sectors (the ranking MAP.322 names labels by) and frames that. No DOM and no
// three.js: tests/js/regionlanding.test.mjs holds it.

// How many of the heaviest sectors make up one cluster, and how far apart (in sector edges, centre to centre) two of them may be
// and still be in the same cluster.
export const CLUSTER_SECTORS = 9;
export const CLUSTER_EDGES = 2.5;

// The cluster of the region's heaviest sectors, framed: {center: [x, y, z] pc, radius pc, count, massKg}, or null when the region
// has no ranking (no sector of it is in `sectors` with a mass). `sectors` is {x, y, z, sizePc, massKg} (the labels' sectors:
// the ones holding a generated sector, `slice.massKg` their total mass); `region` is {center, radius} (the sphere the region's
// extent bounds); `edgePc` is a sector's edge. The cluster is the group within CLUSTER_EDGES edges of one sector with the most
// mass among the CLUSTER_SECTORS heaviest sectors inside the region, taken about its mass-weighted middle.
export function densestCluster(sectors, region, edgePc) {
  if (!region || !region.center || !(region.radius > 0) || !(edgePc > 0)) return null;
  const reach = region.radius + edgePc;
  const inside = sectors.filter(function (s) {
    return s.massKg > 0 && Math.hypot(s.x - region.center[0], s.y - region.center[1], s.z - region.center[2]) <= reach;
  });
  if (!inside.length) return null;
  const heaviest = inside.slice().sort(function (a, b) { return b.massKg - a.massKg; }).slice(0, CLUSTER_SECTORS);
  let best = null;
  heaviest.forEach(function (seed) {
    const members = heaviest.filter(function (s) {
      return Math.hypot(s.x - seed.x, s.y - seed.y, s.z - seed.z) <= CLUSTER_EDGES * edgePc;
    });
    const massKg = members.reduce(function (sum, s) { return sum + s.massKg; }, 0);
    if (!best || massKg > best.massKg) best = { members: members, massKg: massKg };
  });
  const total = best.massKg;
  const middle = [0, 1, 2].map(function (axis) {
    const key = ["x", "y", "z"][axis];
    return best.members.reduce(function (sum, s) { return sum + s[key] * s.massKg; }, 0) / total;
  });
  // Each member's own half diagonal (a sector is a cube of `sizePc` edge) is inside the frame.
  const radius = best.members.reduce(function (far, s) {
    return Math.max(far, Math.hypot(s.x - middle[0], s.y - middle[1], s.z - middle[2]) + 0.5 * Math.sqrt(3) * (s.sizePc || edgePc));
  }, 0);
  return { center: middle, radius: radius, count: best.members.length, massKg: total };
}
