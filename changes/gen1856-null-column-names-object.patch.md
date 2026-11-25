**Headline:** When the server refuses a body's row for a missing value, the error now names the body and the columns left empty, not only the column.

### Fixed
- A row refused for a NULL in a required column (`rotation_period_hours` came out so once in four runs of a two-worker run, and the error named the column and nothing else) now says the table, the required columns that were left empty, and the object's class and name, kind and place, so the next occurrence points at the body and the path that made it. The error is still the server's own (same class, same number); the generated numbers do not change.

generation-output: unchanged
generation-affects: none
