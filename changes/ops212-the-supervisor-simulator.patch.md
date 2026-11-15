**Headline:** The supervisor is tested under a fake clock with killed supervisors, stopped starters, killed stage ends and refused Redis commands.
 
### Changed
- **OPS.212, output unchanged.** `tests/supervisor_sim.py` draws, from a seed, asks for runs, units that finish, stage ends that are killed, starters that stop, supervisors that are killed and replaced, and Redis commands that are refused, and checks that every ask is answered once with one run, the line is served in order, a finished stage is ended exactly once, a starter's stop is said once, and nothing waits for ever once the faults stop. It follows the stage 0 skeleton (`enginesim.py`: a seed replays exactly, nothing sleeps). It found three faults in OPS.211, fixed there: an ask lost when Redis failed in its start, a second run possible when a start was carried out again, and a run whose starter stopped in its first turn freed as dead.
generation-output: unchanged
generation-affects: none
