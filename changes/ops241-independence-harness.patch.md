**Headline:** Every queue command now declares the tables it reads and writes, and a SQL tap test holds those declarations to what the commands really touch.

### Changed
- New `queue/touches.py` holds the reads and writes of all 66 command kinds and end hooks, the conflict rule (a write another command reads or writes, a shared file or cache), `independence.json` (the proven pairs, empty today) and `plan_edges`, the planner's edges. `planrows.record_plan` takes its `after:` edges from it; with no pair listed each step still waits for the one before it, so a plan is exactly as before.
- New `queue/sqltap.py` reads the statements a command sends and names the tables each reads and writes. `tests/test_command_touches.py` runs a toy galaxy and the commands around it under the tap and fails a command that touches an undeclared table, a registered command with no declaration, and a command the run never exercised that is not listed with the reason. Output of generation is unchanged.
