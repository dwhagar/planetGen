**Headline:** A job entering the queue can now declare where it may be cut, whether the engine may choose the cut, and what it depends on; the engine that reads this is in and logs every decision, but no stage uses it yet

### Changed
- **The declaration (OPS.161).** `queue/breakpoints.py` holds what a job declares: its suggested split, the positions a cut may fall at, spans of items that must stay whole, a split flag (`byte-identical`: the engine may choose any legal cut; `suggested`: the job's own cut is used exactly, which is also the default for a job that declares nothing), and its dependencies (`after`, `reads`, `writes`). `CommandSpec.breaks` names a command's declaration; a `whole` command declares itself as one unit; a command that declares nothing keeps its typed split and runs in plan order, as now.
- **The engine (OPS.161).** `queue/splitter.py` chooses the cut of a job flagged byte identical from its profile and the host's free memory and cores: units of about a fifth of the unit ceiling, fewer items if the unit's memory would not fit its share, between a quarter and four times the suggested unit (the typed figure stays the default and the clamp), never inside a span that must stay whole, never more than `MAX_UNITS` units. The decision (flag, old and new units and items, samples, the reason) is written to the stage's plan in the run record when the stage opens and to the run's log. `queue/splitprofile.py` is a stub that knows nothing until the profile table (OPS.155) exists, so every job still takes its suggested split.
- No stage declares a split yet and no output changes; PERF.261 moves the stages onto this one at a time, behind PERF.260's byte-identical test.

generation-output: unchanged
generation-affects: none
