// html/static/crossstars.js
//
// MAP.329: the systems of the sectors round the open one, so that stepping from an open system can cross the sector's edge. The open sector's
// scene holds only its own systems; the neighbours' come from their scene JSON (the one an opened sector is built from) read for the
// positions and the system each star opens, with nothing drawn. Pure: tests/js/crossstars.test.mjs holds it to this description.

// How many neighbouring sectors are read at most, nearest first, and how far (in sector edges, from the system stood at to a sector's middle)
// a sector may be to be read: a system at the edge of its sector has the sector across it about one edge away, and one in a corner the
// diagonal ones within about 1.5.
export const CROSS_SECTORS = 8;
export const CROSS_EDGES = 1.5;

// The sectors to read: of `records` ({id, shell, band, wedge, ...}) all but the open one (`openId`), whose middle (`centreOf(record)`, [x, y, z] pc)
// is within `CROSS_EDGES` edges of `point`, nearest first, at most `CROSS_SECTORS`.
export function neighbourRecords(records, centreOf, point, edgePc, openId) {
  return (records || [])
    .filter((record) => record && record.id != null && record.id !== openId)
    .map((record) => {
      const centre = centreOf(record);
      return { record, distance: centre ? Math.hypot(centre[0] - point[0], centre[1] - point[1], centre[2] - point[2]) : Infinity };
    })
    .filter((found) => found.distance <= CROSS_EDGES * edgePc)
    .sort((a, b) => a.distance - b.distance)
    .slice(0, CROSS_SECTORS)
    .map((found) => found.record);
}

// The systems a sector's scene JSON lists, as {x, y, z, name, endpoint, href, sector} with the position in the galaxy's frame in parsecs: placed as
// the scene is when it is built (sectorscene.js `place`: the scene's own units are `halfEdgePc` over `sceneHalfPx` pixels a half edge, y is flipped
// and the origin is the sector's middle). A companion of a binary is its primary's, so the pair is one system to step to. `sector` is the
// {shell, band, wedge} the system is in.
export function starsOfScene(data, sector) {
  if (!data || !data.centerPc || !(data.halfEdgePc > 0)) return [];
  const unit = data.halfEdgePc / (data.sceneHalfPx || 160);
  const [ox, oy, oz] = data.centerPc;
  const found = [];
  const seen = new Set();
  (data.stars || []).forEach((star) => {
    if (star.companionOf != null || !star.endpoint || String(star.endpoint).indexOf("system:") !== 0 || !star.href) return;
    if (seen.has(star.endpoint)) return;
    seen.add(star.endpoint);
    found.push({
      x: ox + star.x * unit, y: oy - star.y * unit, z: oz + star.z * unit,
      name: star.name, endpoint: star.endpoint, href: star.href, sector: { shell: sector.shell, band: sector.band, wedge: sector.wedge },
    });
  });
  return found;
}
