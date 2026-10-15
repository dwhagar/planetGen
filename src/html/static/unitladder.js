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
const { formatNumber, threeFigures } = await import(`./numberformat.js${VERSION_QUERY}`);
const { LADDERS } = await import(`./unitladders.js${VERSION_QUERY}`);

export { LADDERS };

const factorOf = (ladder, symbol) => ladder.rungs.find((rung) => rung.symbol === symbol).factor;

// Meters, from the table (planetgen/physics/constants.py).
export const KM_M = factorOf(LADDERS.length, "km");
export const AU_M = factorOf(LADDERS.length, "AU");
export const LIGHTYEAR_M = factorOf(LADDERS.length, "ly");
export const PARSEC_M = factorOf(LADDERS.length, "pc");
export const SECONDS_PER_YEAR = factorOf(LADDERS.time, "years");

// `value` (in the rung's base unit) in the rung's unit: three significant
// figures, or a whole number from the rung's `whole_from`; the singular
// form when the number is exactly 1.
function inRung(rung, value) {
  const number = value / rung.factor;
  if (rung.whole_from !== null && Math.abs(value) >= rung.whole_from) {
    return formatNumber(number) + " " + rung.symbol;
  }
  const text = threeFigures(number);
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

// `value` (in `ladder`'s base unit) as text; an en dash for null and NaN.
export function formatQuantity(ladder, value) {
  if (value == null || isNaN(value)) {
    return "–";
  }
  if (!isFinite(value)) {
    return (value > 0 ? "inf " : "-inf ") + ladder.nonfinite;
  }
  const rung = pickRung(ladder, value);
  let text = inRung(rung, value);
  if (rung.with_aside) {
    const size = Math.abs(value);
    const [, symbol] = ladder.asides.find(([start]) => size >= start);
    text += " (" + inRung(ladder.rungs.find((candidate) => candidate.symbol === symbol), value) + ")";
  }
  return text;
}

const scaled = (value, factor) => (value == null ? value : value * factor);

export const formatDistanceM = (meters) => formatQuantity(LADDERS.length, meters);
export const formatDistanceKm = (km) => formatDistanceM(scaled(km, KM_M));
export const formatDistanceAu = (au) => formatDistanceM(scaled(au, AU_M));
export const formatDistanceLy = (ly) => formatDistanceM(scaled(ly, LIGHTYEAR_M));
export const formatDistancePc = (pc) => formatDistanceM(scaled(pc, PARSEC_M));

export const formatSpeedMs = (ms) => formatQuantity(LADDERS.speed, ms);
export const formatSpeedKms = (kms) => formatSpeedMs(scaled(kms, KM_M));

export const formatDurationSeconds = (seconds) => formatQuantity(LADDERS.time, seconds);
export const formatPeriodYears = (years) => formatDurationSeconds(scaled(years, SECONDS_PER_YEAR));
