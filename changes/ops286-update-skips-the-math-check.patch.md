**Headline:** `update.sh` no longer runs the math check as a step, so an update is about 14 s (warm) to 49 s (cold) shorter (OPS.286).

### Changed

- The update has 8 stages instead of 9: "Checking the generator's math" is gone and the later stages are renumbered. The check took about 1 s but ran as an engine run, costing 14 s on a warm box and about 49 s on a cold one, and its result was only a warning.
- The population pass offer no longer waits on the check, and the closing message no longer repeats its warning. `planetgen check-math`, the admin page's button, the web app's startup check and the test gate are unchanged.
- INSTALL.md and the update's header comment say so.
