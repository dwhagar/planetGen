### Changed

- Every page, map panel and generated description now shows a number with
  5 or more digits before the decimal point in scientific notation, to 3
  significant figures: "384,400 km" reads "3.84 × 10⁵ km", "12,345
  systems" reads "1.23 × 10⁴ systems" (UX.20). Counts and measurements
  both follow it. Left as they were: IDs, seeds, years in dates, page
  numbers, ring and slot numbers, coordinates and galaxy positions, and
  the raw numbers in the API's JSON. One shared formatter does it
  (`stellarObjects.utils.format_number`, `num()` in the templates, and
  `static/numberformat.js` for the maps). No regenerate is needed: pages
  render descriptions from the stored numbers.
