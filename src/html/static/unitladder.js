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
const { LADDERS, TINY_NUMBER } = await import(`./unitladders.js${VERSION_QUERY}`);

export { LADDERS };

const factorOf = (ladder, symbol) => ladder.rungs.find((rung) => rung.symbol === symbol).factor;

// Meters, from the table (planetgen/physics/constants.py).
export const KM_M = factorOf(LADDERS.length, "km");
export const AU_M = factorOf(LADDERS.length, "AU");
export const LIGHTYEAR_M = factorOf(LADDERS.length, "ly");
export const PARSEC_M = factorOf(LADDERS.length, "pc");
export const SECONDS_PER_YEAR = factorOf(LADDERS.time, "years");
const companionOf = (ladder, symbol) => ladder.companions.find((companion) => companion.symbol === symbol).factor;
export const STANDARD_ATMOSPHERE_PA = companionOf(LADDERS.pressure, "atm");
export const STANDARD_GRAVITY = companionOf(LADDERS.gravity, "g");
export const SOLAR_MASS_KG = factorOf(LADDERS.mass, "solar masses");
export const SOLAR_LUMINOSITY_W = factorOf(LADDERS.power, "L☉");

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
// preference.
export function formatQuantity(ladder, value, { spoken = false, companions = true, system = "auto" } = {}) {
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

const scaled = (value, factor) => (value == null ? value : value * factor);

export const formatDistanceM = (meters) => formatQuantity(LADDERS.length, meters);
export const formatDistanceKm = (km) => formatDistanceM(scaled(km, KM_M));
export const formatDistanceKmSpoken = (km) => formatQuantity(LADDERS.length, scaled(km, KM_M), { spoken: true });
export const formatDistanceAu = (au) => formatDistanceM(scaled(au, AU_M));
export const formatDistanceLy = (ly) => formatDistanceM(scaled(ly, LIGHTYEAR_M));
export const formatDistancePc = (pc) => formatDistanceM(scaled(pc, PARSEC_M));

export const formatSpeedMs = (ms) => formatQuantity(LADDERS.speed, ms);
export const formatSpeedKms = (kms) => formatSpeedMs(scaled(kms, KM_M));

export const formatDurationSeconds = (seconds) => formatQuantity(LADDERS.time, seconds);
export const formatPeriodYears = (years) => formatDurationSeconds(scaled(years, SECONDS_PER_YEAR));

// A surface temperature with °C and °F beside it, or in K alone.
export const formatTemperatureK = (kelvin, companions = true) => formatQuantity(LADDERS.temperature, kelvin, { companions });

export const formatPressurePa = (pascals) => formatQuantity(LADDERS.pressure, pascals);
export const formatPressureAtm = (atm) => formatPressurePa(scaled(atm, STANDARD_ATMOSPHERE_PA));

// "3.71 m/s² (0.378 g)" for a gravity in multiples of standard gravity.
export const formatGravityG = (g) => formatQuantity(LADDERS.gravity, scaled(g, STANDARD_GRAVITY));

// A mass in kilograms: "1 Earth mass", "0.299 Jupiter masses", "0.122 solar masses".
export const formatMassKg = (kg) => formatQuantity(LADDERS.mass, kg);
export const formatMassSolar = (solar) => formatMassKg(scaled(solar, SOLAR_MASS_KG));

// A luminosity or power in watts: "47 TW", "1 L☉"; a heat flow in W/m²: "87 mW/m²".
export const formatPowerW = (watts) => formatQuantity(LADDERS.power, watts);
export const formatLuminositySol = (solar) => formatPowerW(scaled(solar, SOLAR_LUMINOSITY_W));
export const formatFluxWm2 = (flux) => formatQuantity(LADDERS.flux, flux);
