**Headline:** A cache warm-up now goes on when its terminal is closed, and `planetgen warm --follow` watches a running one from another session.

### Changed
- A lost terminal (an SSH session that dropped, a hangup) no longer ends the cache warm-up as canceled (OPS.95 had done so). The command stops writing to the terminal, notes in the run's log that the terminal was closed, and the warm-up carries on to its end. SIGTERM and Ctrl+C still end it.

### Added
- `planetgen warm --follow` draws the progress of the warm-up that is running, started by `update.sh` or by hand, from its run record in Redis, the same bars and phase lines the starting session shows. It starts nothing; ending it with Ctrl+C or a hangup leaves the warm-up running. It says so when none is running.
