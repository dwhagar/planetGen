// html/static/regionlabels.js
//
// The names of charted regions on the Galaxy Map, by zoom (MAP.314): which names the nodes the map holds offer to the label
// engine (maplabels.js) through the budget (labelbudget.js, domain "region"). A region is a 9-block with at least one generated
// sector, named from its block path (API.35); its 81-block is its family ("Fardoo Reach"). What shows depends on the size of
// the nodes the view is cut into (the node cut, galaxynodecut.js):
//
//   81 and 243   the family word with Reach, one name per node: the family holding the most generated sectors, at the
//                centroid of its charted regions' spheres (weighted by their sectors), not at the middle of the block, so a
//                block mostly empty at the galaxy's edge does not float its name in empty space;
//   27 and 9     the Region name, one per charted region, at the middle of the region's bounding sphere (a 27-block lists its
//                regions, a 9-block is one);
//   3 and below  nothing: the region is the line of the info card (MAP.315) there.
//
// A node with no charted region offers nothing, so an unfilled galaxy shows no name, and no name exists for a block that has
// not been loaded. A region offered by a 27-block is offered once however many nodes list it. Every name is `first` (the
// farthest zooms name regions before anything else) and is ranked by the share of the region charted.
//
// No DOM at import, no three.js: tests/js/regionlabels.test.mjs.

// The node sizes that name a family, and those that name a region.
export const FAMILY_LEVELS = [81, 243];
export const REGION_LEVELS = [9, 27];

// The block size (sectors across) of a node key "m.shell.band.wedge", or NaN.
export function levelOf(key) {
  return parseInt(String(key).split(".")[0], 10);
}

function item(id, text, centre, share, extent, key) {
  return {
    id: id, text: text, x: centre[0], y: centre[1], z: centre[2], ref: { kind: "region", key: key, extent: extent },
    slice: { domain: "region", shareCharted: share / 100, first: true },
  };
}

// The family a list of regions names: the one with the most generated sectors (ties by name), at the sector-weighted centroid
// of its regions' centres, and the mean share charted of its regions. `null` for a list with no named region.
export function familyOf(regions) {
  const families = new Map();
  regions.forEach(function (region) {
    if (!region.parent || !region.centre) return;
    const family = families.get(region.parent) || { name: region.parent, sectors: 0, sum: [0, 0, 0], share: 0, count: 0 };
    const weight = Math.max(region.sectors || 0, 1);
    family.sectors += weight;
    family.share += region.percent || 0;
    family.count += 1;
    for (let axis = 0; axis < 3; axis++) family.sum[axis] += weight * region.centre[axis];
    families.set(region.parent, family);
  });
  let best = null;
  families.forEach(function (family) {
    if (!best || family.sectors > best.sectors || (family.sectors === best.sectors && family.name < best.name)) best = family;
  });
  if (!best) return null;
  return { name: best.name, sectors: best.sectors, share: best.share / best.count,
    centre: best.sum.map(function (total) { return total / best.sectors; }) };
}

// The names to offer for `entries`, `[{key, node}]` (the nodes of the view; a node not loaded yet has no `node`): items of
// `{id, text, x, y, z, ref: {kind: "region", key, extent}, slice}` for the engine, `x`, `y`, `z` in parsecs. `ref.key` is the
// key of the block a click on the name flies to (the node's own, or the region's 9-block) and `ref.extent` the sphere
// `{centre, radius}` to frame there, or `null`.
export function regionLabels(entries) {
  const found = new Map();
  entries.forEach(function (entry) {
    const node = entry.node;
    if (!node) return;
    const level = levelOf(entry.key);
    if (level === 9) {
      if (!node.name || !node.extent) return;
      const share = node.charted ? node.charted.percent : 0;
      found.set("region:" + entry.key, item("region:" + entry.key, node.name, node.extent.centre, share, node.extent, entry.key));
    } else if (level === 27 && node.regions) {
      node.regions.forEach(function (region) {
        if (!region.name || !region.centre) return;
        const extent = region.radius > 0 ? { centre: region.centre, radius: region.radius } : null;
        found.set("region:" + region.key, item("region:" + region.key, region.name, region.centre, region.percent || 0, extent, region.key));
      });
    } else if (FAMILY_LEVELS.indexOf(level) >= 0 && node.regions) {
      const family = familyOf(node.regions);
      if (family) found.set("family:" + entry.key, item("family:" + entry.key, family.name, family.centre, family.share, null, entry.key));
    }
  });
  return Array.from(found.values());
}
