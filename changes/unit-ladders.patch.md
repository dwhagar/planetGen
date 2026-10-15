### Changed

- Distances, speeds and times now share one table of units (`planetgen/util/ladders.py`) and one algorithm in Python and in the browser, so later units (mass, temperature, pressure and the rest) add a table, not a new formatter. The browser's three ladder files are replaced by `unitladder.js` and a generated `unitladders.js`. The text is the same as before except two cases that disagreed between the two languages: a half-way kilometer count ("9,999.5 km") now rounds away from zero in both, and a missing or NaN value is an en dash in the browser too.
