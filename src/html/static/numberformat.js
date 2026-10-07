// html/static/numberformat.js
//
// The site's one number formatter for the browser (UX.20): the mirror of
// stellarObjects/utils.py `format_number` and `scientific_text`. A number
// that would show 7 or more digits with no decimals, or 5 or more digits
// before the decimal point with decimals (UX.36), is shown in scientific
// notation with 3 significant figures ("1.23 × 10⁶"); anything shorter is
// comma-grouped with the requested decimals. Counts and
// measurements alike go through it; IDs, page numbers and coordinates
// don't. tests/test_number_format.py checks both copies against one table.

export const SCIENTIFIC_MIN_INTEGER_DIGITS = 7;
export const SCIENTIFIC_MIN_DECIMAL_INTEGER_DIGITS = 5;
export const SCIENTIFIC_SIGNIFICANT_FIGURES = 3;

const SUPERSCRIPTS = { "-": "⁻", "0": "⁰", "1": "¹", "2": "²", "3": "³", "4": "⁴", "5": "⁵", "6": "⁶", "7": "⁷", "8": "⁸", "9": "⁹" };

// `value` as "1.23 × 10⁶".
export function scientificText(value, significant = SCIENTIFIC_SIGNIFICANT_FIGURES) {
  const [mantissa, exponent] = value.toExponential(significant - 1).split("e");
  const exp = String(parseInt(exponent, 10)).replace(/[-0-9]/g, (c) => SUPERSCRIPTS[c]);
  return mantissa + " × 10" + exp;
}

// `value` with `minDecimals`..`maxDecimals` decimals (default whole),
// comma-grouped, or scientificText(value) when it shows too many digits.
export function formatNumber(value, maxDecimals = 0, minDecimals = maxDecimals) {
  if (value == null || !isFinite(value)) {
    return String(value);
  }
  const text = value.toLocaleString("en-US", {
    maximumFractionDigits: maxDecimals, minimumFractionDigits: minDecimals,
  });
  if (showsTooManyDigits(text)) {
    return scientificText(value);
  }
  return text;
}

// True when formatted `text` should be scientific instead: 7 or more whole
// digits with no decimals shown, 5 or more with decimals (UX.36).
function showsTooManyDigits(text) {
  const [whole, decimals] = text.replace(/^[-+]/, "").split(".");
  const limit = decimals ? SCIENTIFIC_MIN_DECIMAL_INTEGER_DIGITS : SCIENTIFIC_MIN_INTEGER_DIGITS;
  return whole.replace(/,/g, "").length >= limit;
}

// Three significant figures, comma-grouped, trailing zeros dropped;
// scientific past 6 whole digits (UX.36). The mirror of
// stellarObjects/utils.py `_three_figures`; distance.js, speed.js and
// period.js share it.
export function threeFigures(value) {
  if (value === 0 || !isFinite(value)) {
    return String(value);
  }
  const magnitude = Math.floor(Math.log10(Math.abs(value)));
  const decimals = Math.max(0, 2 - magnitude);
  const text = value.toLocaleString("en-US", { maximumFractionDigits: decimals, minimumFractionDigits: 0 });
  if (showsTooManyDigits(text)) {
    return scientificText(value);
  }
  return text;
}
