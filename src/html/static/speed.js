// html/static/speed.js
//
// The speed ladder (UX.13): the browser mirror of stellarObjects/utils.py
// `format_speed_kms`. A speed is shown in one unit, slowest to fastest
// km/h < km/s < Mm/s < c: km/h below 1 km/s, km/s up to 1,000 km/s, Mm/s
// up to a tenth of light speed, then multiples of c ("36 km/h",
// "29.8 km/s", "4.5 Mm/s", "0.25 c"). Keep the constants and rules in step
// with the Python helper; tests/test_speed_period_format.py checks both.

// Imported with this module's own `?v=` query, as systemmap.js explains.
const VERSION_QUERY = new URL(import.meta.url).search;
const { threeFigures } = await import(`./numberformat.js${VERSION_QUERY}`);

// Exact (the SI definition), km/s (stellarObjects/physical_constants.py).
export const SPEED_OF_LIGHT_KMS = 299792.458;

// [label, km/s per unit, shown from this many km/s up], slowest first.
export const SPEED_LADDER = [
  ["km/h", 1 / 3600, 0],
  ["km/s", 1, 1],
  ["Mm/s", 1e3, 1e3],
  ["c", SPEED_OF_LIGHT_KMS, 0.1 * SPEED_OF_LIGHT_KMS],
];

export function formatSpeedKms(kms) {
  if (kms == null || isNaN(kms)) {
    return "–";
  }
  if (!isFinite(kms)) {
    return (kms > 0 ? "inf" : "-inf") + " km/s";
  }
  const size = Math.abs(kms);
  let [label, unitKms] = SPEED_LADDER[0];
  for (const [candidateLabel, candidateKms, thresholdKms] of SPEED_LADDER) {
    if (size >= thresholdKms * (1 - 1e-12)) {
      label = candidateLabel;
      unitKms = candidateKms;
    }
  }
  return threeFigures(kms / unitKms) + " " + label;
}

export function formatSpeedMs(ms) {
  return formatSpeedKms(ms == null ? ms : ms / 1e3);
}
