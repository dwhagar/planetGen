**Headline:** Stars now have companions at every separation from 0.05 AU outward, from one distribution over orbital period, instead of a close pair or one 50 AU and wider with nothing between.

### Changed
- A star's companion chance is shared between three ranges by Moe and Di Stefano's (2017, Table 13) companions per decade of period for its mass: close (0.05 to 0.25 AU, a P-type pair), inner (0.25 to 10 AU) and outer (10 AU and out). The coin flip between close and wide (`WIDE_BINARY_DEFAULT_CHANCE`) is gone. For a G dwarf about 5.7% of its companions are inside 0.25 AU, 24% from 0.25 to 10 AU and the rest beyond; massive stars have more close ones.
- A wide pair now starts at 0.25 AU (`WIDE_BINARY_SEPARATION_MIN_AU`, was 50 AU) and its axis is drawn from the same lines up to the place's survival limit, so pairs of a few AU exist and their planets and comets are limited to the fraction of an AU the Holman and Wiegert limit leaves. A forced comet with no room around a tight pair is reported as an unmet requirement (`StarSystem.unmet_requirements`).
- The metallicity trend of GEN.103 now moves the odds of the close and inner ranges together (Moe, Kratter and Badenes 2019 measured periods under 10^4 days, 10 AU); the outer range stays flat in [Fe/H]. A place that keeps no pair wider than 0.25 AU loses the wide chance; one that keeps part of the range loses that part.
- The usual shares on the Generate page and in the prevalence checks are remeasured over 20,000 systems: systems with two stars 0.260 (was 0.303; the metallicity boost now applies to the inside-10-AU part only, not to the whole close half), binaries that are wide 0.928 (was 0.499, the close class being 0.05 to 0.25 AU only), systems with comets 0.321 (was 0.334), habitable world 0.257 (was 0.242), asteroid belt 0.607 (was 0.590).
- New `generation/companion_classes.py`, `tuning.COMPANIONS_PER_DECADE_BY_PRIMARY_MASS` and `tuning.COMPANION_LOG_PERIODS_DAYS`; `docs/design/multistar-and-compact-systems.md` 4 and `placement-by-population.md` 3.1 to 3.3 describe the model and its levels (a G dwarf inside 10 AU is about 0.6 of the pooled line's level; the cause is open).

generation-output: changed

TODO items resolved: GEN.448 / Added: the inside-10-AU level against the pooled line (0.13 against 0.20 at [Fe/H] 0) and the O and B flatness ratios (1.32, 1.43) for TEST.196 and Boss
