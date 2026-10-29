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
const { PLANET_CLASS_LETTERS } = await import(`./galaxymarkers.js${VERSION_QUERY}`);

// The distance from the system's centre, in radii of the system, at which the body symbols come in and at which they go (the camera
// that has just flown in to fit the system sits about four radii away, so the symbols are there from the first view).
export const HANDOFF_ON = 6;
export const HANDOFF_OFF = 10;

// The reasons (ids of INTEREST_REASONS) a body of the system scene has, strongest first, each `{reason, tier}` with the tier index of its
// habitability (0 Blue to 3 Red, or null); `[]` for a body with no symbol (a star, a comet, an unscored terrestrial body with no life and
// no class). Mirrors planetgen.db.interests.REASONS for one body, in the same order: an inhabited world, a civilization, sentient life,
// life, a habitable world (Yellow or better), a barren scored one, a gas giant, its planet class; a belt is an asteroid field (MAP.321,
// MAP.311, MAP.220, MAP.343). The life fields are the stored scores as they are (0 none to 9 abundant, or missing).
export function reasonsOfBody(body) {
  if (!body) return [];
  if (body.kind === "belt") return [{ reason: "asteroid-field", tier: null }];
  if (body.kind !== "planet" && body.kind !== "moon") return [];
  const tier = body.tier ? TIER_NAMES.indexOf(body.tier) : -1;
  const tierIndex = tier >= 0 ? tier : null;
  const found = [];
  const add = function (reason, withTier) { found.push({ reason: reason, tier: withTier ? tierIndex : null }); };
  if (body.inhabited) add("inhabited", true);
  if (body.life_stage === "technological_civilization") add("civilization", false);
  if (body.life_sentient >= 1) add("sentient", false);
  if (body.life_microbial >= 1 || body.life_flora >= 1 || body.life_fauna >= 1 || body.life_sentient >= 1) add("life", false);
  if (body.body_type === "g") {
    add("gas-giant", false);
  } else if (tierIndex !== null && tierIndex <= 2) {
    add("habitable", true);
  } else if (tierIndex !== null) {
    add("uninhabitable", true);
  }
  if (body.planet_class && PLANET_CLASS_LETTERS.indexOf(body.planet_class) >= 0) add("class-" + body.planet_class, false);
  return found;
}

// The strongest reason of a body, or null: `reasonsOfBody(body)[0]`.
export function reasonOfBody(body) {
  return reasonsOfBody(body)[0] || null;
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
    // The strongest reason whose layer is on names the body's symbol (as a star's mark is drawn once, by the first reason that is on).
    const sort = reasonsOfBody(body).find(function (one) { const kind = byReason.get(one.reason); return kind && isOn(kind.id); });
    const kind = sort && byReason.get(sort.reason);
    if (!kind) return;
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
