**Headline:** A command the update runs (the cache warm-up, the math check) keeps the supervisor it asked for, so it no longer waits a minute and fails.

### Fixed
- OPS.298: the supervisor ended itself after one turn whenever an update was running, even when the update's own command had started it. `check-math` and `warm` inside `update.sh` then found no supervisor for a minute and the math check stopped with "no supervisor has been alive for 1m". The supervisor now stops only for an update it is not part of.

generation-output: unchanged
