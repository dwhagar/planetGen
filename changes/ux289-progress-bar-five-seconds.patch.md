**Headline:** A progress bar now appears for any step that takes more than 5 seconds, and a random start's first sector has one.

### Changed
- **The bar threshold is 5 seconds (UX.289).** `tuning.PROGRESS_BAR_SECONDS` was 15 seconds; a step that is predicted to take over 5 seconds, or has no recorded speed yet, draws its bar at once, and a step that was thought quick but runs past 5 seconds gets its bar at that moment, started at the work done so far and counting the rest. The docs of `steps` and `stage_progress` say 5 seconds.
- **A random start's sector has a bar.** Stage 1 of `planetgen galaxy` ("Generate the starting sector") took about 90 seconds on a small galaxy with nothing drawn; it now shows a "Starting sector" bar (or "Starting sectors" with `--neighborhoods`) like the single-address mode already did.
- **Tests.** `test_steps.py` checks the 5 second threshold, a prediction just under and just over it, a quick step that runs past 5 seconds (with a fake timer) and a quick step that leaves no bar.
