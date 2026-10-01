// html/static/period.js
//
// The time-period ladder (UX.14): the browser mirror of
// stellarObjects/utils.py `format_duration_seconds` and
// `format_period_years`. A period is shown in the largest unit it is at
// least 1 of, µs < ms < s < minutes < hours < days < years < ky < My < Gy,
// to three significant figures, singular when the shown value is exactly 1
// ("1 day", "27.3 days", "1.88 years", "236 My"). A year is the Julian
// 365.25 days. Keep the constants and rules in step with the Python
// helper; tests/test_speed_period_format.py checks both.

// Imported with this module's own `?v=` query, as systemmap.js explains.
const VERSION_QUERY = new URL(import.meta.url).search;
const { threeFigures } = await import(`./numberformat.js${VERSION_QUERY}`);

export const SECONDS_PER_YEAR = 365.25 * 24 * 3600;

// [plural label, singular label, seconds], shortest first.
export const PERIOD_LADDER = [
  ["µs", "µs", 1e-6],
  ["ms", "ms", 1e-3],
  ["s", "s", 1],
  ["minutes", "minute", 60],
  ["hours", "hour", 3600],
  ["days", "day", 86400],
  ["years", "year", SECONDS_PER_YEAR],
  ["ky", "ky", SECONDS_PER_YEAR * 1e3],
  ["My", "My", SECONDS_PER_YEAR * 1e6],
  ["Gy", "Gy", SECONDS_PER_YEAR * 1e9],
];

export function formatDurationSeconds(seconds) {
  if (seconds == null || isNaN(seconds)) {
    return "–";
  }
  if (!isFinite(seconds)) {
    return (seconds > 0 ? "inf" : "-inf") + " years";
  }
  if (seconds === 0) {
    return "0 s";
  }
  const size = Math.abs(seconds);
  let [plural, singular, unitS] = PERIOD_LADDER[0];
  for (const candidate of PERIOD_LADDER) {
    if (size >= candidate[2] * (1 - 1e-12)) {
      [plural, singular, unitS] = candidate;
    }
  }
  const number = threeFigures(seconds / unitS);
  return number + " " + (number === "1" ? singular : plural);
}

export function formatPeriodYears(years) {
  return formatDurationSeconds(years == null ? years : years * SECONDS_PER_YEAR);
}
