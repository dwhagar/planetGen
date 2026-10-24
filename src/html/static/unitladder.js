// html/static/unitladder.js
//
// The site's one unit formatter for the browser (UX.23): the mirror of
// planetgen/util/format.py `format_quantity`. A value in a ladder's SI base
// unit is shown in the largest unit on the ladder it is at least 1 of, to
// three significant figures ("384,400 km", "1.52 AU", "29.8 km/s",
// "27.3 days"); a parsec value adds one familiar unit in parentheses
// ("4.2 pc (13.7 ly)"). The ladders are the generated unitladders.js
// (scripts/build_unit_ladders.py, from planetgen/util/ladders.py), so the
// two languages share one table; tests/test_unit_ladders.py checks that the
// copies give the same text.

// Imported with this module's own `?v=` query, as systemmap.js explains.
const VERSION_QUERY = new URL(import.meta.url).search;
const { formatNumber, scientificText, threeFigures } = await import(`./numberformat.js${VERSION_QUERY}`);
const { LADDER_OF_KIND, LADDERS, TINY_NUMBER } = await import(`./unitladders.js${VERSION_QUERY}`);
const { toBase, unit } = await import(`./unitregistry.js${VERSION_QUERY}`);
const { unitSystem } = await import(`./unitpreference.js${VERSION_QUERY}`);

export { LADDERS };

// Three significant figures, or scientific notation under TINY_NUMBER, where a
// unit would otherwise run to a row of zeros (planetgen/util/format.py
// `_three_figures_or_tiny`).
function threeFiguresOrTiny(number) {
  return number !== 0 && Math.abs(number) < TINY_NUMBER ? scientificText(number) : threeFigures(number);
}

// `number` to whole units, or one decimal when under 10 in size ("15",
// "5.9", "0"): `_whole_or_tenths`.
function wholeOrTenths(number) {
  let text = formatNumber(number, Math.abs(number) < 10 ? 1 : 0);
  if (text.endsWith(".0")) {
    text = text.slice(0, -2);
  }
  return text === "-0" || text === "+0" ? "0" : text;
}

// True when a companion unit is shown beside `value`: its size is in the
// companion's range (`_companion_applies`).
function companionApplies(companion, value) {
  const size = Math.abs(value);
  if (companion.min_size !== null && size < companion.min_size) {
    return false;
  }
  return companion.max_size === null || size < companion.max_size;
}

// `value` (in the rung's base unit) in the rung's unit: three significant
// figures, or a whole number from the rung's `whole_from`, or scientific
// notation under TINY_NUMBER of the unit (UX.81); the singular form when the
// number is exactly 1. `spoken` writes the unit's name, for screen readers.
function inRung(rung, value, spoken) {
  const number = (value - rung.offset) / rung.factor;
  if (rung.whole_from !== null && Math.abs(value) >= rung.whole_from) {
    return formatNumber(number) + " " + (spoken ? rung.name : rung.symbol);
  }
  const text = threeFiguresOrTiny(number);
  if (spoken) {
    return text + " " + rung.name;
  }
  return text + " " + (rung.singular && text === "1" ? rung.singular : rung.symbol);
}

// The rung `value` is shown in: the largest one it is at least `start` of,
// with a hair of tolerance (see planetgen/util/format.py `pick_rung`).
export function pickRung(ladder, value) {
  const size = Math.abs(value);
  if (size === 0 && ladder.zero !== null) {
    return ladder.rungs.find((rung) => rung.symbol === ladder.zero);
  }
  let picked = ladder.rungs[0];
  for (const candidate of ladder.rungs) {
    if (size >= candidate.start * (1 - 1e-12)) {
      picked = candidate;
    }
  }
  return picked;
}

// The ladder as it reads under a unit preference ("auto", "metric" or
// "customary", UX.78): the mirror of planetgen/util/ladders.py `for_system`.
const systemLadders = new Map();
export function forSystem(ladder, system = "auto") {
  const key = ladder.base + "|" + system;
  if (!systemLadders.has(key)) {
    const rungs = ladder.rungs.filter((rung) => rung.systems.includes(system));
    const symbols = new Set(rungs.map((rung) => rung.symbol));
    let asides = ladder.asides.filter(([, symbol]) => symbols.has(symbol));
    if (!rungs.some((rung) => rung.with_aside)) {
      asides = [];
    }
    if (asides.length && asides[asides.length - 1][0] !== 0) {
      asides = asides.slice(0, -1).concat([[0, asides[asides.length - 1][1]]]);
    }
    systemLadders.set(key, {
      ...ladder,
      rungs,
      asides,
      zero: symbols.has(ladder.zero) ? ladder.zero : null,
      companions: ladder.companions.filter((companion) => companion.systems.includes(system)),
    });
  }
  return systemLadders.get(key);
}

// `value` (in `ladder`'s base unit) as text; an en dash for null and NaN.
// `spoken` spells the units out ("1.5 astronomical units"), for
// `aria-valuetext` and the like; `companions: false` leaves out the units
// beside the primary one (°C and °F beside K); `system` is the reader's unit
// preference, the saved one (unitpreference.js) unless given.
export function formatQuantity(ladder, value, { spoken = false, companions = true, system = unitSystem() } = {}) {
  if (value == null || isNaN(value)) {
    return "–";
  }
  if (!isFinite(value)) {
    return (value > 0 ? "inf " : "-inf ") + ladder.nonfinite;
  }
  ladder = forSystem(ladder, system);
  const rung = pickRung(ladder, value);
  let text = inRung(rung, value, spoken);
  if (rung.with_aside) {
    const size = Math.abs(value);
    const [, symbol] = ladder.asides.find(([start]) => size >= start);
    text += " (" + inRung(ladder.rungs.find((candidate) => candidate.symbol === symbol), value, spoken) + ")";
  }
  if (companions) {
    const shown = ladder.companions.filter((companion) => companionApplies(companion, value)).map((companion) => {
      const number = (value - companion.offset) / companion.factor;
      return (companion.whole ? wholeOrTenths(number) : threeFiguresOrTiny(number)) + " " + (spoken ? companion.name : companion.symbol);
    });
    if (shown.length) {
      text += " (" + shown.join(", ") + ")";
    }
  }
  return text;
}

// A planet, moon or star radius: always km in scientific notation (Boss,
// 2026-09-30), the one exception to the distance ladder; the mirror of
// planetgen/util/format.py `format_body_radius_km`.
export const formatBodyRadiusKm = (km) => (km == null || isNaN(km) ? "–" : scientificText(km) + " km");

// The ladder the registry unit `key` is shown on.
const ladderOf = (key) => LADDERS[LADDER_OF_KIND[unit(key).kind]];

// `formatQuantity` for a Quantity (quantity.js): the ladder its kind is shown on and its value in the
// registry's base unit, so a value is shown the same whatever unit it was held in ("384,400 km" for
// 384.4 Mm). The mirror of planetgen/util/format.py `format_measure`.
export function formatMeasure(quantity, options) {
  const ladder = LADDERS[LADDER_OF_KIND[quantity.kind]];
  if (ladder === undefined) throw new Error(`no ladder shows a ${quantity.kind}`);
  return formatQuantity(ladder, quantity.base(), options);
}

// `formatQuantity` for `value` of the registry unit `key`, the factor the registry's; null gives the en
// dash. What the format names below are (planetgen/util/format.py `_format_in`).
const formatIn = (key, value, options) => formatQuantity(ladderOf(key), value == null ? value : toBase(value, key), options);

export const formatDistanceM = (meters) => formatIn("m", meters);
export const formatDistanceKm = (km) => formatIn("km", km);
export const formatDistanceKmSpoken = (km) => formatIn("km", km, { spoken: true });
export const formatDistanceAu = (au) => formatIn("au", au);
export const formatDistanceLy = (ly) => formatIn("ly", ly);
export const formatDistancePc = (pc) => formatIn("pc", pc);

export const formatSpeedMs = (ms) => formatIn("m_s", ms);
export const formatSpeedKms = (kms) => formatIn("km_s", kms);

export const formatDurationSeconds = (seconds) => formatIn("s", seconds);
export const formatPeriodYears = (years) => formatIn("yr", years);

// A surface temperature with °C and °F beside it, or in K alone.
export const formatTemperatureK = (kelvin, companions = true) => formatIn("K", kelvin, { companions });

export const formatPressurePa = (pascals) => formatIn("Pa", pascals);
export const formatPressureAtm = (atm) => formatIn("atm", atm);

// "3.71 m/s² (0.378 g)" for a gravity in multiples of standard gravity.
export const formatGravityG = (g) => formatIn("g0", g);

// A mass in kilograms: "1 Earth mass", "0.299 Jupiter masses", "0.122 solar masses".
export const formatMassKg = (kg) => formatIn("kg", kg);
export const formatMassSolar = (solar) => formatIn("Msun", solar);

// A luminosity or power in watts: "47 TW", "1 L☉"; a heat flow in W/m²: "87 mW/m²".
export const formatPowerW = (watts) => formatIn("W", watts);
export const formatLuminositySol = (solar) => formatIn("Lsun", solar);
export const formatFluxWm2 = (flux) => formatIn("W_m2", flux);
