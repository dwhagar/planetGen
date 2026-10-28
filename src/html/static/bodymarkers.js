// html/static/bodymarkers.js
//
// MAP.321: the symbols of the worlds themselves on the Galaxy Map's system level. Out in the galaxy a star carries the symbol of
// its best world (galaxymarkers.js, from the nodes' `interests`); once the camera is in the system, the star does not repeat it
// and each planet, moon or belt carries its own symbol, by the same reasons and the same Menu layers. This is the choosing half
// with no DOM and no three.js (tests/js/bodymarkers.test.mjs holds it): which reason a body of the system scene is, the marker
// items for the bodies, and the hand-off between the two levels, which has a gap between the distance it turns on at and the one
// it turns off at so a camera resting at the boundary does not flicker (the same idea as MAP.152's ambient hand-off).

const VERSION_QUERY = new URL(import.meta.url).search;
const { TIER_NAMES } = await import(`./poiicons.js${VERSION_QUERY}`);

// The distance from the system's centre, in radii of the system, at which the body symbols come in and at which they go (the camera
// that has just flown in to fit the system sits about four radii away, so the symbols are there from the first view).
export const HANDOFF_ON = 6;
export const HANDOFF_OFF = 10;

// The reason (an id of INTEREST_REASONS) a body of the system scene carries its symbol for, with the tier index of its habitability
// (0 Blue to 3 Red, or null), or null for a body with no symbol (a star, a comet, an unscored terrestrial body). Mirrors planetgen.db.interests.REASONS for one body: an inhabited world first, then a habitable one (Yellow
// or better), a gas giant, a barren scored world; a belt is an asteroid field.
export function reasonOfBody(body) {
  if (!body) return null;
  if (body.kind === "belt") return { reason: "asteroid-field", tier: null };
  if (body.kind !== "planet" && body.kind !== "moon") return null;
  const tier = body.tier ? TIER_NAMES.indexOf(body.tier) : -1;
  const tierIndex = tier >= 0 ? tier : null;
  if (body.inhabited) return { reason: "inhabited", tier: tierIndex };
  if (body.body_type === "g") return { reason: "gas-giant", tier: null };
  if (tierIndex !== null && tierIndex <= 2) return { reason: "habitable", tier: tierIndex };
  if (tierIndex !== null) return { reason: "uninhabitable", tier: tierIndex };
  return null;
}

// The marker items for the bodies of the open system: `kinds` the interest kinds (galaxymarkers.INTEREST_KINDS), `isOn(kindId)`
// whether its Menu layer is on and `entries` the system stage's bodies (galaxysystem.js: `{ref, kind, name, body, x, y, z}`,
// x, y and z kept up to date as the body moves). Each item is the marker layer's `{kind, item}` with `item.follow` the entry, so
// the symbol rides along, and `item.bodyRef` the body's ref, so a click selects the body.
export function bodyMarkers(kinds, isOn, entries) {
  const byReason = new Map();
  kinds.forEach(function (kind) { byReason.set(kind.interest, kind); });
  const found = [];
  (entries || []).forEach(function (entry) {
    const body = Object.assign({ kind: entry.kind }, entry.body);
    const sort = reasonOfBody(body);
    const kind = sort && byReason.get(sort.reason);
    if (!kind || !isOn(kind.id)) return;
    const item = { key: entry.ref, id: entry.ref, name: entry.name, x: entry.x, y: entry.y, z: entry.z, n: 1, interest: sort.reason, bodyRef: entry.ref, follow: entry };
    if (sort.tier !== null) item.tier = sort.tier;
    found.push({ kind: bodyKind(kind), item: item });
  });
  return found;
}

// The kind of a body's symbol: the reason's kind with the hover text of one body ("Earth, habitable world, Green").
const bodyKinds = new Map();
function bodyKind(kind) {
  if (!bodyKinds.has(kind.id)) {
    bodyKinds.set(kind.id, Object.assign({}, kind, {
      describe: function (item) {
        return item.name + ", " + (kind.bodyName || kind.interest) + (item.tier == null ? "" : ", " + TIER_NAMES[item.tier]);
      },
    }));
  }
  return bodyKinds.get(kind.id);
}

// A hand-off between the star-level symbol and the body-level ones: `update(distance, radius, live)` is true while the body
// symbols show. It turns on when the system is open (`live`) and the camera is within HANDOFF_ON radii, and stays on until the camera is
// beyond HANDOFF_OFF radii or the system closes.
export function createHandoff() {
  let on = false;
  return {
    update: function (distance, radius, live) {
      if (!live || !(radius > 0)) {
        on = false;
        return on;
      }
      const ratio = distance / radius;
      if (on) {
        if (ratio >= HANDOFF_OFF) on = false;
      } else if (ratio <= HANDOFF_ON) {
        on = true;
      }
      return on;
    },
    isOn: function () { return on; },
  };
}
