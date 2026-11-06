// html/static/sensorview.js
//
// The Sensor view's painter table (UX.275, planetgen/web/lib/sensor.py): which painter draws the record `GET /api/objects/<ref>/sensor` returns.
// This file only names the painters that exist and the module of each; it paints nothing itself (the pictures are unchanged). The ids are
// the server's `PAINTERS` (tests/test_sensor.py checks the two lists agree).

export const SENSOR_VERSION = 1;

export const PAINTERS = {
  sphere: { module: "bodyview.js", draws: "a star, planet or moon look record as a lit sphere" },
  phenomenon: { module: "phenomenonrender.js", draws: "a neutron star, black hole, quasar, rogue planet or interstellar comet in 3D" },
  nebula: { module: "nebulaview.js", draws: "a nebula's cloud from its points" },
};

// The painter that draws a Sensor record, or null where the record has none (it says why in `reason`) or is from a newer server than this script.
export function painterOf(record) {
  if (!record || typeof record !== "object" || record.v > SENSOR_VERSION || !record.painter) return null;
  return Object.prototype.hasOwnProperty.call(PAINTERS, record.painter) ? { id: record.painter, ...PAINTERS[record.painter] } : null;
}
