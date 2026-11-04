**Headline:** The "Hours outside in a day" line of a body's survival equipment now goes through the unit ladder, so a short stay reads in minutes.

### Changed
- **UX.258.** The describer wrote the hours a person can spend outside in a day with a raw format (`3.9 hours`). It now writes the quantity through the unit registry like every other describer fact (`_held("h", ...)`), so the figure is the same and carries its unit as the registry shows it: `4 hours`, `3.96 hours`, and `59.2 minutes` where the stay is under an hour. The wording and the "set by" cause are unchanged. Only the describer text changes; no stored value does.
generation-output: unchanged
