**Headline:** `planetgen sectors --block` fills one block of the Galaxy Map's ladder, and an admin API route starts the same fill as a job.

### Changed
- **Fill a block of the ladder (GEN.310).** `planetgen sectors --block <m.shell.band.wedge>` (and `planetgen galaxy --block`) generates every not-yet-generated sector the stored outline holds inside that block (243, 81, 27, 9 or 3), inner shells first, with the same size estimate, progress, cancel, `--limit` and `--yes` as a shell or a span. The old refusal ("blocks are being rebuilt on shells") is gone.
- **The Galaxy Map's fill route (GEN.310).** `POST /api/galaxy/blocks/<key>/generate` (admin) starts the same fill as a queued job with the bright-star backfill and settle after it, or answers an estimate with `"estimate_only": true`; it answers 413 for a block of more than 2 million sectors.
