// html/static/systemgroups.js
//
// MAP.428: the stars of one system, shown as one system. A wide binary or a multiple draws each of its stars as a point of its own (the
// Galaxy Map's star nodes place every star), and as the view comes close the points part and each is named: nothing tells the eye they
// are one system. This works out, from the stars the map draws (`{id, system_id, name, x, y, z, luminosity_sol, ...}`), the groups of two
// or more stars of one system, the system's one name, and the link lines that join the members on screen. Pure (no DOM, no three.js):
// tests/js/systemgroups.test.mjs holds it to this description.

// The members of a group are this far apart on screen (px) or more before the link is drawn: closer they are one dot to the eye already.
export const MIN_LINK_PX = 4;
// The link is not drawn for a pair this far apart on screen or more: the stars of a group that wide are on separate parts of the screen,
// and a line across it says less than it hides.
export const MAX_LINK_PX = 900;
// The most groups linked in one frame (the brightest first).
export const MAX_LINKED_GROUPS = 120;

// The groups of `stars` of one system: [{systemId, members (brightest first), luminositySol}], for each system with two or more stars
// of the list. A star with no system (`system_id` null or undefined) is no group's member.
export function groupBySystem(stars) {
  const bySystem = new Map();
  (stars || []).forEach((star) => {
    if (!star || star.system_id == null || star.phenomenon) return;
    if (!bySystem.has(star.system_id)) bySystem.set(star.system_id, []);
    bySystem.get(star.system_id).push(star);
  });
  const groups = [];
  bySystem.forEach((members, systemId) => {
    if (members.length < 2) return;
    members.sort((a, b) => (b.luminosity_sol || 0) - (a.luminosity_sol || 0));
    groups.push({ systemId, members, luminositySol: members.reduce((sum, star) => sum + (star.luminosity_sol || 0), 0) });
  });
  groups.sort((a, b) => b.luminositySol - a.luminositySol);
  return groups;
}

// The one name of a group's system: the name of its members with the component letter off ("Ofriggor A" and "Ofriggor D" are "Ofriggor"),
// else the first name the brightest member has. "" when no member is named.
export function systemNameOf(group) {
  const letter = /^(.+\S) [A-Z]$/;
  const names = group.members.map((star) => star.name).filter(Boolean);
  if (!names.length) return "";
  const stripped = names.map((name) => {
    const found = letter.exec(name);
    return found ? found[1] : name;
  });
  // The primary may carry the system's own name already ("Ofriggor") while its companions carry letters: the shortest stripped name is the system's.
  return stripped.reduce((best, name) => (name.length < best.length ? name : best), stripped[0]);
}

// Whether the member's own name is only the system's (so the group's label stands for it): any named member of a group is named by the group.
export function labelledByGroup(star, groupsBySystem) {
  return star.system_id != null && groupsBySystem.has(star.system_id);
}

// The link lines of one group on screen: from the brightest member to each of the others, as [{from: {x, y}, to: {x, y}, id}], or [] while
// the members are closer than MIN_LINK_PX (one dot), any is further than MAX_LINK_PX, or fewer than two are on screen at all.
// `project(star)` is the star's {x, y, visible} in screen px.
export function linksOf(group, project) {
  const points = group.members.map((star) => ({ star, at: project(star) })).filter((item) => item.at && item.at.visible);
  if (points.length < 2) return [];
  const first = points[0].at;
  const lines = [];
  let widest = 0;
  points.slice(1).forEach((item) => {
    const gap = Math.hypot(item.at.x - first.x, item.at.y - first.y);
    widest = Math.max(widest, gap);
    lines.push({ id: group.systemId + ":" + item.star.id, from: { x: first.x, y: first.y }, to: { x: item.at.x, y: item.at.y }, gap });
  });
  if (widest < MIN_LINK_PX || widest > MAX_LINK_PX) return [];
  return lines;
}

// Where the group's shared name is anchored: the middle of its members (parsecs), as [x, y, z].
export function groupCenter(group) {
  const n = group.members.length;
  return [
    group.members.reduce((sum, star) => sum + star.x, 0) / n,
    group.members.reduce((sum, star) => sum + star.y, 0) / n,
    group.members.reduce((sum, star) => sum + star.z, 0) / n,
  ];
}
