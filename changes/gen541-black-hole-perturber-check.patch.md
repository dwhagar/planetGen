**Headline:** The black holes' share of the wide-pair perturber mass is checked against four papers: the generator's black hole count and masses lie inside the published ranges, so the remnants stay in.

### Changed
- `docs/design/placement-by-population.md` records the check of GEN.541: the share of core-collapse stars that leave a hole (0.29, against 0.26 to 0.45 in Sukhbold et al. 2016), the mean hole mass (12.5 against 9 to 14 there and 11.7 to 14.2 in Olejak et al. 2020) and the number for a galaxy of 5e10 solar masses (1.6e8 to 2.5e8 against Olejak et al.'s 1.3e8 and the older 1e8 to 1e9). By the widest-range rule the draw stays; the perturber mass, 1.44 solar masses (stars alone 0.55), would be 0.75 to 1.34 for the other published count and mass shapes, with the holes carrying 0.2 to 0.6 of the second moment.
- The note's two stale figures (0.594 and 1.46, from before the Hurley, Pols and Tout lifetimes) are 0.553 and 1.44.
- `tests/test_place_terms.py` records the three ranges that hold and the two perturber masses.

TODO items resolved: GEN.541 / Added: none
