**Headline:** a unit that has just left a long wait is no longer stopped as stuck, and the heartbeat has one test that runs every rule on a fake clock.

### Changed

- OPS.143: `tests/test_heartbeat_rules.py` drives the whole chain (a horse that beats, the follower's watch, the stuck action, the stale rule, the job rows) over hours of fake time: work, slow, a long wait, quiet and stopped, the false kill (a moving counter at ten times its estimate for six hours is flagged and never touched), dead, stuck twice, the stale rule for every mix of states, and the API, web and tile jobs.
- Fixed, found by that test: the time a unit spent waiting on the gate, the sky lock or MySQL counted as time its counter stood still, so a unit that left a wait of more than five minutes was judged quiet at once and stopped after two polls. The horse now writes the completed waits since the counter last changed (`waited`) and the judge subtracts them.
- Fixed, found by that test: a unit queued again by the pass that found its stop had taken effect was looked up in the verdicts the watch had just forgotten, which raised a `KeyError` in the follower; it is skipped now.
