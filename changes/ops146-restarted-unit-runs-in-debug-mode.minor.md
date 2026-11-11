**Headline:** A unit the heartbeat stopped as stuck and queued again now runs in debug mode and writes a debug file of its own, with the first attempt's stack at the top, so the second attempt can be understood; on by default, bounded and rotated (OPS.146, heartbeat part 9).

### Added
- **`queue.debug_restarted_units`** (config.json, on by default, with the other queue keys; on the Settings page as "Debug file for restarted units"). With it off nothing extra is written.
- **A debug file for each restarted unit**: `unit-debug-<job id>-attempt<N>.log` in the log folder (the job id is the run, stage and unit; `N` is 2, the restart). It is separate from `error.log` and the run's log. The follower makes it when it queues the unit again (`queue/unitdebug.py`), with why the unit was stopped and the stack of every thread of the first attempt at the top, and the unit's record names it; the horse of the restarted attempt then runs its work under `log.debug_to_file`, which appends everything it logs at debug level: phases, counter changes at every beat, waits for the gate, the sky lock and MySQL with how long they took (`unitbeat`), SQL statements with their timing, and the memory it holds at every beat. The random draws are left out (they are many times the rest).
- **Bounded and rotated.** The program rotates the file at 20 MB and keeps three older copies; the name ends in `.log`, so the folder's one logrotate rule covers it too; the newest 40 such files are kept.
- **The run's log gets one line** naming the file ("unit N is run in debug mode on its restart: see ...").
- `log.debug_to_file(path, areas, rotate, create_mode)`: a second user of the debug machinery `log.relayed` already used (they share `_debugging`); it adds to a debug run's output and takes nothing from it.

### Changed
- `unitbeat.UnitBeat` says in debug what it does (phase begins, waits for and waited, each beat with the counter and memory); these lines go nowhere unless a debug sink is on.
- No star, planet or stored number changes.

generation-output: unchanged
generation-affects: none
