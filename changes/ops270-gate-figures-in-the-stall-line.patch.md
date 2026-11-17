**Headline:** The stall line of a stage that still has units to start says what the memory gate holds, and the new `record-phenomenon-scatter` command has a memory figure.

### Added
- When no unit of a stage has finished for 30 minutes and units are still to start, the log line now ends with the gate's figures: the claims held by the bulk lanes (what they claimed and what their processes hold), the cores they use of the cores the gate allows, the memory the host has available of its total, what the gate keeps back, and the load a core. A stage that runs fewer units than it has workers shows what limits it (OPS.270).

### Fixed
- `record-phenomenon-scatter` had no memory figure, so `test_every_command_with_units_has_a_memory_figure` failed on main; it claims 700 MB like the other one-unit commands.
