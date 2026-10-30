**Headline:** A measured Solar System target can no longer be a row of a generator table, and a test lists every row that rests on a Solar System body.

### Changed

- TEST.215: a possibility-table row has a `calibration_only` flag, and the registry refuses a table that holds one, so the Earth, Mars, Moon, Galilean-moon and Sirius B targets live in a test file and never in a table. No registered row had it set, so nothing moved and no table value, weight or digest changes. A new test scans every registered row's source and limit for a Solar System body and requires each to be a flagged stand-in with a replacement item (checked), a row marked pending for a decision, or a row reviewed as not resting on a body (with its reason); a new row that names a body and is in none of them fails. The pending rows are listed for the TODO thread to decide.
