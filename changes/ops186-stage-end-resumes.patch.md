**Headline:** A run no longer waits for ever when the unit that finished a stage is killed before it queues the next one.

### Changed
- **OPS.186, output unchanged.** The end of a stage (its conclusion, its row and the queuing of the next opening) now records how far it got and can be finished by the retry of the unit that began it: `units.repair_unended_stages` finds a finished stage that did not advance and ends it, concluding and queuing each only once, under a lease so two attempts do not end it together. Queuing the opening is tried four times before it fails the run, and a unit whose end is unfinished with no retry left fails the run, with the stage named, instead of "the run goes on". The supervisor's reconciler (OPS.210) calls the same function.
- **OPS.187, output unchanged.** An opening that has cut its stage marks it done (`runrecord.opened_key`); a failure RQ reports for that job after that no longer fails the run, and its retry does not cut the stage again.
generation-output: unchanged
generation-affects: none
