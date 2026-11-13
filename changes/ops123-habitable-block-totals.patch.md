**Headline:** Making the blocks' habitable-world totals follow the sectors' counts is a step of its own after the habitability re-score, no longer the hidden end of its last range (OPS.123 slice 2).

### Changed
- **`make-habitable-block-totals` is its own command and stage.** The `conclude` of `rescore-habitability` is gone; its effect (`block_stats.rebuild`, MAP.143, PERF.202) is a one-unit stage, "Make the habitable block totals" (`habitability-blocks`), after the re-score's ranges. A `planetgen rescore-habitability` run and the Regenerate habitability job count one more stage. The step runs whether or not any range was stale, as the `conclude` did. Reused: `block_stats.rebuild`, the one-unit whole expansion and the stage pattern of PRs 6349 and 6356.

generation-output: unchanged
generation-affects: none
