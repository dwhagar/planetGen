**Headline:** The browser holds measured values in a Quantity and reads every unit factor from the registry.

### Added
- **The browser's Quantity** (UX.128). `static/quantity.js` mirrors the Python `Quantity`: a number in one registry unit, answered in any other unit of its kind, with the same bounds checks. `unitladder.js` shows one with `formatMeasure`.

### Changed
- The unit formatter's `format*` names read the registry unit they are for, and the factors it exported (`KM_M`, `AU_M`, `LIGHTYEAR_M`, `PARSEC_M`, `SECONDS_PER_YEAR`, `STANDARD_GRAVITY` and the rest) are removed.
- The light-year to parsec conversion of the galaxy map and the generate buttons, the speed steps of the orbit clock, the time estimate of the limit sliders, the duration words of the settings page and the Sun's radius in the system scale go through the registry.
