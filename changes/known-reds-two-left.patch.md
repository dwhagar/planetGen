**Headline:** The known-red list no longer names test_work_units and test_gen_reach_groups, which pass on main.

### Changed

- `scripts/data/known_reds.json`: TEST.642 (test_work_units) and TEST.650 (test_gen_reach_groups) are removed; both pass on main after Train 21, so the gate treats a failure in them as new.

generation-affects: none
