// html/static/numberformat.js
//
// The site's one number formatter for the browser (UX.20): the mirror of
// stellarObjects/utils.py `format_number` and `scientific_text`. A number
// that would show 5 or more digits before the decimal point is shown in
// scientific notation with 3 significant figures ("1.23 × 10⁶"); anything
// shorter is comma-grouped with the requested decimals. Counts and
// measurements alike go through it; IDs, page numbers and coordinates
// don't. tests/test_number_format.py checks both copies against one table.

export const SCIENTIFIC_MIN_INTEGER_DIGITS = 5;
export const SCIENTIFIC_SIGNIFICANT_FIGURES = 3;

const SUPERSCRIPTS = { "-": "⁻", "0": "⁰", "1": "¹", "2": "²", "3": "³", "4": "⁴", "5": "⁵", "6": "⁶", "7": "⁷", "8": "⁸", "9": "⁹" };

// `value` as "1.23 × 10⁶".
export function scientificText(value, significant = SCIENTIFIC_SIGNIFICANT_FIGURES) {
  const [mantissa, exponent] = value.toExponential(significant - 1).split("e");
  const exp = String(parseInt(exponent, 10)).replace(/[-0-9]/g, (c) => SUPERSCRIPTS[c]);
  return mantissa + " × 10" + exp;
}

// `value` with `minDecimals`..`maxDecimals` decimals (default whole),
// comma-grouped, or scientificText(value) past 4 whole digits.
export function formatNumber(value, maxDecimals = 0, minDecimals = maxDecimals) {
  if (value == null || !isFinite(value)) {
    return String(value);
  }
  const text = value.toLocaleString("en-US", {
    maximumFractionDigits: maxDecimals, minimumFractionDigits: minDecimals,
  });
  const whole = text.replace(/^[-+]/, "").split(".")[0].replace(/,/g, "");
  if (whole.length >= SCIENTIFIC_MIN_INTEGER_DIGITS) {
    return scientificText(value);
  }
  return text;
}
