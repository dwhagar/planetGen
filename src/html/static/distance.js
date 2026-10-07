// html/static/distance.js
//
// The distance ladder for the maps: the browser mirror of
// planetgen/util/format.py `format_distance_m`. A value is shown in the
// largest unit it is at least 1 of (km < AU < mpc < cpc < ly < pc < kpc <
// Mpc < Gpc); a parsec-family value adds one parenthetical, in ly when it
// is at least 0.01 ly, else AU when at least 0.01 AU, else km ("4.2 pc
// (13.7 ly)", "2.4 mpc (495 AU)"). Keep the constants and rules in step
// with the Python helper; tests/test_distance_format.py checks both.

// Imported with this module's own `?v=` query, as systemmap.js explains.
const VERSION_QUERY = new URL(import.meta.url).search;
const { formatNumber, threeFigures } = await import(`./numberformat.js${VERSION_QUERY}`);

// Boss's exact values, in meters (planetgen/physics/constants.py).
export const KM_M = 1e3;
export const AU_M = 149597870700;
export const LIGHTYEAR_M = 9460730472580800;
export const PARSEC_M = 3.085677581491367e16;

const LADDER = [
  ["km", KM_M],
  ["AU", AU_M],
  ["mpc", PARSEC_M * 1e-3],
  ["cpc", PARSEC_M * 1e-2],
  ["ly", LIGHTYEAR_M],
  ["pc", PARSEC_M],
  ["kpc", PARSEC_M * 1e3],
  ["Mpc", PARSEC_M * 1e6],
  ["Gpc", PARSEC_M * 1e9],
];
const PARSEC_UNITS = new Set(["mpc", "cpc", "pc", "kpc", "Mpc", "Gpc"]);
export const DISTANCE_PAREN_MIN_LY = 0.01;
export const DISTANCE_PAREN_MIN_AU = 0.01;

function inUnit(meters, label, unitM) {
  if (label === "km" && Math.abs(meters) >= 1e6) {
    return formatNumber(Math.round(meters / unitM)) + " km";
  }
  return threeFigures(meters / unitM) + " " + label;
}

export function formatDistanceM(meters) {
  if (meters == null || isNaN(meters)) {
    return "unknown";
  }
  const size = Math.abs(meters);
  let label = LADDER[0][0];
  let unitM = LADDER[0][1];
  for (const [candidateLabel, candidateM] of LADDER) {
    if (size >= candidateM * (1 - 1e-12)) {
      label = candidateLabel;
      unitM = candidateM;
    }
  }
  let text = inUnit(meters, label, unitM);
  if (PARSEC_UNITS.has(label)) {
    text += " (" + distanceParenthetical(meters) + ")";
  }
  return text;
}

// The familiar unit a parsec-family value carries in parentheses.
export function distanceParenthetical(meters) {
  const size = Math.abs(meters);
  if (size >= DISTANCE_PAREN_MIN_LY * LIGHTYEAR_M) {
    return inUnit(meters, "ly", LIGHTYEAR_M);
  }
  if (size >= DISTANCE_PAREN_MIN_AU * AU_M) {
    return inUnit(meters, "AU", AU_M);
  }
  return inUnit(meters, "km", KM_M);
}

export function formatDistanceKm(km) {
  return formatDistanceM(km == null ? km : km * KM_M);
}

export function formatDistanceAu(au) {
  return formatDistanceM(au == null ? au : au * AU_M);
}

export function formatDistanceLy(ly) {
  return formatDistanceM(ly == null ? ly : ly * LIGHTYEAR_M);
}

export function formatDistancePc(pc) {
  return formatDistanceM(pc == null ? pc : pc * PARSEC_M);
}
