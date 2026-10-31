**Headline:** The 23 flavour lines a system can carry are rewritten in plain words: short sentences, one idea each, and no claim about a world the data could contradict.

### Changed
- **System flavour lines (GEN.1212).** `tuning.SYSTEM_FLAVOR` keeps its 23 lines in the same order, so the draw lands on the same line as before, but each line is now one or two short sentences ("Sensors show a derelict chain of deep-space platforms. They fire strong energy beams into the void for no known reason."). A line that called a world "now-dead" no longer does. A test reads every line after "Sensors show" against the style guide. A system's stored `system_flavor_text` is the new wording.

No change to which systems get a line or which line they get.

generation-output: changed

The stored `system_flavor_text` wording is the only output that changes (epoch affects: fill). It rides a merge train, which takes the epoch.

TODO items resolved: GEN.1212 (system flavour lines; the stored planet description copy is decided separately with Core 5).
