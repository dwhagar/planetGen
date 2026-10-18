**Headline:** A disk star's age is drawn from the galaxy's star formation rate inside its age window, so old and middle-aged stars cluster where the galaxy formed most, not evenly across the window.

### Changed
- A young, intermediate or old disk star's age follows the star formation history's rate inside its window (`galaxy.age_law`, from `galaxy.disk_history.population_age_laws`) instead of an even spread. The scatter's bright stars, the mass pass, both backfills and a sector's own dim systems all use the same law, so their counts and ages agree. The bulge, halo, clusters and satellites keep the even spread.
- The local B-star share between the arms falls from about 0.069% to 0.050% (observed 0.04%), so `test_the_local_b_share_reaches_the_observed_one` is no longer an expected failure (it allows 1.5 times the observed share).
- `docs/design/disk-star-formation-history.md` describes the law.

generation-output: changed

TODO items resolved: GEN.546, GEN.547 / Added: none
