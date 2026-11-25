**Headline:** A row with a missing value for a required column is refused where it is made, and the error names the body and the columns left empty.

### Fixed
- A row refused for a NULL in a required column (`rotation_period_hours` came out so once in four runs of a two-worker run, and the error named the column and nothing else) now says the table, the required columns that were left empty, and the object's class and name, kind and place, so the next occurrence points at the body and the path that made it. The check is made as the row is built, because a sector save holds its INSERTs back and the server only refuses at the flush, after the body is gone. The error is the server's own class and number (1048); the generated numbers do not change.

generation-output: unchanged
generation-affects: none
