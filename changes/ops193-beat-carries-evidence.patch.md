**Headline:** The beat of a running job now carries proof that it is working, so a long database statement no longer looks like a job that stopped.

### Changed
- **Evidence in the beat (OPS.193).** Every ten seconds a unit's beat also writes the processor seconds it has used, the bytes it has read and written (`/proc/self/io`), and, while a database statement is out, the statement's label and the time it was sent. The database layer marks the start and end of every statement it sends (`db/statementmark.py`, reused from the `store` connection), because MySQL's own progress figures are not readable by the account the install creates. Nothing reads the new fields yet: OPS.194 teaches the judge to call a unit that moves them busy instead of quiet.

generation-output: unchanged
generation-affects: none
