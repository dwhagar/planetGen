**Headline:** A run is started by the supervisor, not by the process that asks for it.

### Changed
- **OPS.211, output unchanged.** `planetgen generate`, the Generate page and the API no longer take the machine themselves. A run leaves its ask on Redis (a place in the waiting line for a run that waits, a request of its own for a check or maintenance run, which is refused at once with the holder named when the machine is taken) and the supervisor's `RunStarts` duty starts it: the first of the line when the machine is free. The process that asked follows the run it was given (`runrecord.start_run(owner=...)` names it), so a run is never started by a process that is gone a moment later. `runline.start_when_free` carries the ask and waits for the supervisor's answer (`runline.verdict`); when no supervisor answers within a minute it fails naming the supervisor (`runline.NotStarted`) instead of waiting for ever.
generation-output: unchanged
generation-affects: none
