**Headline:** A run is now called stale only when every unit that is running is dead or stuck, so one long unit whose count moves no longer shows a "Clear the stale run" offer after two minutes (OPS.141, heartbeat part 8).

### Changed
- **`runrecord.is_stale` reads the unit states.** `read_run` now carries, for each stage, the states the follower last published for its running units (`stage["unit_states"]`, from `unitwatch`, OPS.138). While any unit is judged, the run is stale exactly when every one of them is `dead` or `quiet`; a unit that works, is slow, loops, or waits on a lock keeps the run alive however long it goes without finishing a unit or logging an event. Before this, a run was stale after `max(120 s, 3 x its longest unit so far)` without a finished unit or an event, which one 30 minute sky cell always tripped.
- With no unit judged (none has started, or an older process follows the run) the old rule stands, so a run with queued units and no worker is still called stale after the same time. A run whose follower is gone is still found by `holder_is_dead`, which the Clear offer reads too.
- No star, planet or stored number changes.

generation-output: unchanged
generation-affects: none
