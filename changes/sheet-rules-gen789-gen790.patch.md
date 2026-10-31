**Headline:** The star property sheet now says which stars have a widest stable orbit (the stars of a wide pair and the further members of a system, not a single star or a close pair), and the comet rule is confirmed.

### Changed
- **GEN.789.** The sheet rule for `wide_binary_a_crit_km` now matches what the generator stores: a wide pair's two stars and every member beyond the pair have the limit; a single star has no companion and the stars of a close pair are orbited as one body, so they have none. The open gap that excused the column is gone, and the fill audit reads the system's pairing when it checks a star.
- **GEN.790.** No change needed: the comet rule (`mean_anomaly_deg`, `orbital_period_years`, `period_class` apply to elliptical comets only) already matches the generator, which fills all three for every elliptical comet and none for a parabolic one (checked on 300 comets).
generation-output: unchanged
