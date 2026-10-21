**Headline:** A Generate job's log no longer says the stage twice: the banner above a one-stage step names it, so the two lines under it that repeated the stage are gone.

### Changed
- The log of a step of a web job that has one stage printed `=== Stage 4 of 13: Plan the galaxy ===` (the banner) and then `1 stage (stages 4 to 4 of 13): 4. Plan the galaxy.` and `Stage 4 of 13: Plan the galaxy.`, which say the same. Under a banner a step of one stage that is not skipped now prints neither. A skipped stage still says why, a step of several stages still lists its stages and announces each, a command-line run with no banner keeps its lines, and the Generate page's stage text and the progress file are unchanged.

TODO items resolved: PERF.145 / Added: none
