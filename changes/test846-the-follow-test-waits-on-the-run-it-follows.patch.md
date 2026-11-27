**Headline:** The follow-the-warm-up test waits for the warm-up's run to be written and for the follower to be watching, not for a clock

### Changed

- TEST.846: `test_follow_draws_the_warm_up_that_another_process_runs_and_starts_nothing` waited for any maintenance run and then half a second, and held the family for 30 s; on a loaded machine the follower could come before the run had its title and steps (it then says "No cache warm-up is running" and ends) or the held family let go before anyone watched. The test now waits until the run is the warm-up's with its steps written, releases the family only when the follower has attached (`wrapper.follow` called), and holds it until then with no short timeout. No retry.

generation-output: unchanged
