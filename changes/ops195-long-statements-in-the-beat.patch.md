**Headline:** A long database statement of a command shows in the beat with its label and start, and the stage events of MySQL are used where the account may read them.

### Changed
- **One helper for a long statement (OPS.195).** `queue/longsql.run(conn, statement, label)` sends the statement and names it in the beat while it runs (`db/statementmark.py`, the unit's phase). Where the database account can read MySQL's current stage events the beat also carries the work done and estimated for the statement; the account the install creates cannot, so this is checked once and stays off.

generation-output: unchanged
generation-affects: none
