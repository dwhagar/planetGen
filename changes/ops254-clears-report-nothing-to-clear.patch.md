**Headline:** The clear of the earlier bright stars still shows its step on the Queue page when a new galaxy has nothing to clear.

### Fixed
- **OPS.254 (bug), output unchanged.** After the clears of a new galaxy were made cheap, the clear unit reported nothing, so its bar never showed a first step and the log had no line for it (`test_unit_reports` went red on main). A clear that finds nothing stored now says so as its one phase ("No bright stars are stored, nothing to clear") and the bar fills.
generation-output: unchanged
generation-affects: none
