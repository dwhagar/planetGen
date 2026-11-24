**Headline:** A young star now holds no planets and no belts while it still has its own disk, and holds them once the disk is gone, instead of none under a fixed 10 million years.

### Changed

- The fixed cut of 10 Myr (`PLANET_MIN_STAR_AGE_GY`, gone) is the star's own disk (GEN.1552): `StarSystem._holds_its_disk` is true while the star's `disk_kind` is class 0, I or II (which its age and its own dispersal time, GEN.1551, gave it), and then `_apply_star_history`, the trojan and belt-promotion passes leave the star with no planets and no belts; a request for a class or a `SLOTS` list overrides it as before. A star whose disk is gone before 10 Myr (about 6 in 10 stars of 3 Myr) now stands with its bodies where it was empty.
- `pre_main_sequence.dispersal_yr` is the star's own time again, no longer held at 10 Myr, so a few stars hold a disk past 10 Myr as the survey curve has it (about 2 in 100 for a solar-mass star).
- `tests/audit/prevalence_disk_age.py` (the disk-and-planet census by star age, built ahead under GEN.1553) comes with it.

generation-output: changed
generation-affects: fill
