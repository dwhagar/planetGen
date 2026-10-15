### Added
- The admin Stats page shows Redis's health: whether it answers, whether its append-only file is on (off means a restart loses queued work), its memory policy, the memory used, the last save and how many workers are alive. A Redis that does not answer is reported with the reason, and the check never waits more than a second for it.
