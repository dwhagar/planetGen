// html/static/groupwindow.js
//
// MAP.222: the group window of the Galaxy Map. Choosing a grouping of things (a nebula with the nebulae that overlap it; a globular
// cluster, MAP.223) flies the camera to frame it and shows that grouping alone: the map draws only what lies inside the window's bounding
// sphere, whatever sectors the group crosses, and a way back returns to where the camera was. A window is plain data: a kind, a key, a
// title, a bounding sphere and its members; each kind says how its subject becomes a window and how it is told in words, so a second
// kind plugs in with `registerKind`.
//
// No DOM and no three.js here: tests/js/groupwindow.test.mjs holds it to this description.

// The window frames the bounding sphere this many times over, and a member joins a nebula's window only when its sphere reaches the
// subject's; the window never grows past MAX_GROWTH times the subject's radius, so a chain of touching nebulae does not take the map.
export const FRAME_MARGIN = 1.15;
export const MAX_GROWTH = 3;
// A thing just outside the sphere is kept this share of the radius beyond it (so a star at the edge is not cut in half).
export const EDGE_SLACK = 0.02;

const kinds = new Map();

// Adds a kind: `spec = { label(subject) -> string, members(subject, candidates) -> list of spheres {x, y, z, radius_pc}, noun(count) -> string }`.
export function registerKind(kind, spec) {
  kinds.set(kind, spec);
}

export function kindNames() {
  return Array.from(kinds.keys());
}

function distance(a, b) {
  return Math.hypot(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}

// The smallest sphere about the middle of the members' bounding box that holds every member sphere, `{ centre: [x, y, z], radiusPc }`.
export function boundsOf(spheres) {
  if (!spheres.length) return null;
  const lo = [Infinity, Infinity, Infinity];
  const hi = [-Infinity, -Infinity, -Infinity];
  spheres.forEach((s) => {
    [s.x, s.y, s.z].forEach((v, i) => {
      lo[i] = Math.min(lo[i], v - s.radius_pc);
      hi[i] = Math.max(hi[i], v + s.radius_pc);
    });
  });
  const centre = lo.map((v, i) => (v + hi[i]) / 2);
  const radiusPc = Math.max(...spheres.map((s) => distance(centre, [s.x, s.y, s.z]) + s.radius_pc));
  return { centre, radiusPc };
}

// The window for `subject` of `kind`, built from `candidates` (what the map draws now, for the kind to choose members from), or null for
// an unknown kind or a subject with no extent.
export function createWindow(kind, subject, candidates) {
  const spec = kinds.get(kind);
  if (!spec || !subject) return null;
  let members = spec.members(subject, candidates || []);
  let bounds = boundsOf(members);
  if (!bounds) return null;
  const own = spec.sphere ? spec.sphere(subject) : subject;
  if (bounds.radiusPc > MAX_GROWTH * own.radius_pc) {
    members = [subject];
    bounds = boundsOf([own]);
  }
  return {
    kind,
    key: spec.key ? spec.key(subject) : String(subject.id),
    title: spec.label(subject),
    centre: bounds.centre,
    radiusPc: bounds.radiusPc,
    members,
    count: members.length,
    summary: spec.noun(members.length),
  };
}

// Whether the point (x, y, z) lies in the window.
export function contains(win, x, y, z) {
  if (!win) return true;
  return distance([x, y, z], win.centre) <= win.radiusPc * (1 + EDGE_SLACK);
}

// The sphere the camera frames: the window's bounds with the margin.
export function framing(win) {
  return { center: win.centre.slice(), radius: win.radiusPc * FRAME_MARGIN };
}

// The nebula kind: the subject is a cloud (`{ type, id, name, descriptor, x, y, z, radius_pc }`); its members are itself and every candidate
// cloud whose sphere reaches its own.
registerKind("nebula", {
  key: (cloud) => cloud.type + ":" + cloud.id,
  label: (cloud) => cloud.name || (cloud.type === "nebula" ? "Nebula" : "Supernova remnant"),
  sphere: (cloud) => ({ x: cloud.x, y: cloud.y, z: cloud.z, radius_pc: cloud.radius_pc }),
  members: (cloud, candidates) => {
    const own = [cloud.x, cloud.y, cloud.z];
    const found = [cloud];
    candidates.forEach((other) => {
      if (!other || (other.type === cloud.type && other.id === cloud.id)) return;
      if (distance(own, [other.x, other.y, other.z]) <= cloud.radius_pc + other.radius_pc) found.push(other);
    });
    return found;
  },
  noun: (count) => (count === 1 ? "1 nebula" : count + " nebulae"),
});

// The globular cluster kind (MAP.223): the subject is a cluster of the map's scene data (`{ ref, name, x, y, z, radius_pc }`) and its window is its own
// sphere, whatever sectors that crosses: the stars the map draws inside it are the cluster's members (the stored ones of the filled sectors and, where a
// sector is not yet filled, the plan's, drawn from the cluster's King profile), and nothing else is. A cluster has no neighbours to join it.
registerKind("cluster", {
  key: (cluster) => cluster.ref,
  label: (cluster) => cluster.name || "Globular cluster",
  sphere: (cluster) => ({ x: cluster.x, y: cluster.y, z: cluster.z, radius_pc: cluster.radius_pc }),
  members: (cluster) => [{ x: cluster.x, y: cluster.y, z: cluster.z, radius_pc: cluster.radius_pc }],
  noun: () => "a globular cluster",
});
