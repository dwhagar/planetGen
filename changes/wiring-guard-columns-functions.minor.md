**Headline:** The wiring guard now also finds a stored column nothing reads and a public function only tests call, and today's findings are named in the register.

### Added
- **Two new kinds in `scripts/wiring_guard.py` (TEST.520).** `column`: a column of the tables that hold objects that no module outside the plumbing (schema, models, column units, store, storage checksum, checks, export and repair code, the generic Details sheet) reads, by its own name, by the name the store gives it (`object_rows.attr`, a dictionary or an attribute assigned from the row) or by an f-string that builds it. A physics or generation module that loads the value is a reader, so a column read downstream is not listed. `function`: a public module-level function of a reachable module of `physics`, `generation` or `galaxy` that no other module, no other function of its module and no `"module:function"` string names.
- **`src/tests/test_wiring_guard.py` fails on a new one**, on a register line for something that is read now, and on a named item that is not open, as for the other kinds; seven new tests run the two kinds on a small package.
- **Today's findings in `docs/design/research/wiring-register.md`.** Seven columns (`binary_epoch_unix` and `binary_next_update_due` are by design, the orbit update's own clock; `cluster_id` waits for GEN.359, `dynamical_state` for GEN.291, `formed_at` and `origin_event_id` for GEN.346, `host_node_id` for GEN.332) and 157 functions in 90 modules, each named to TEST.520 until it has an owner or is deleted. Report: `reports/wiring-guard-test520-2026-10-31.md`.

No generator output changes.

generation-output: unchanged

TODO items resolved: none (TEST.520 stays open until each of the 157 function lines has an owner item or is deleted).
