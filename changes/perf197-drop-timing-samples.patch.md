**Headline:** The generator no longer records single-object timing samples or a timing histogram that nothing read.

### Changed

- The `object_timing_samples` table and the `histogram` column of `object_timing_stats` are dropped (migration 0195; a table or column already gone is skipped). The sample list, the 1-in-1,000 sampling, the histogram buckets and the `PLANETGEN_OBJECT_TIMING_FULL` switch are removed from `generation/object_timing.py`. The per-kind count, total, fastest and slowest, and the Stats page's Seconds per object table, are unchanged. generation-output: unchanged
