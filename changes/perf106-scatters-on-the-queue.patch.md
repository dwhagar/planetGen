**Headline:** The luminosity pass and the dimmer band of bright stars run as commands on the work queue, like the massive-star pass.

### Changed
- **The bright-star luminosity pass and the band are commands (PERF.106).** `scatter-luminosity` and `scatter-band` each cut their pass into units of consecutive shells (`planetgen commands` lists them), queue them through the wrapper as a run of their own and follow the units with the same progress bar as before. A unit runs alone from its JSON arguments (`planetgen command-run`).
- **The units are the same at any worker count, and the stars are the same as before.** The old per-pass task queue (`run_plan._scatter_shells`, its progress channel and drain thread) is gone with no replacement; the pass no longer cuts its work by `--workers`. On the test galaxy the stars stored after the scatter and after a band are identical to the ones the previous version drew, at one worker and at four.
- **The luminosity pass reads the sectors the massive pass marked when its stage starts**, between the two passes, as it did when it was told them by the starter.
- **`docs/design/redis-and-the-work-queue.md` describes it.**

TODO items resolved: none (slice one of PERF.106). Added: none.
