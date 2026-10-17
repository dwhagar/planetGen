**Headline:** Mid and early B stars now have companions at the frequencies Moe and Di Stefano measured: 76% and 84%, with more triples and quadruples.

### Changed
- B stars of 5 to 9 solar masses at birth (mid B) have a companion with chance 0.76 and those of 9 to 16 (early B) 0.84, from Table 13's single-star fractions of 0.24 and 0.16, where every B star had 0.65. B stars under 5 solar masses keep 0.65, which sits inside the 0.59 +-0.11 of Table 13's bin for A and late B.
- The star counts of those two classes follow Table 13's single, binary and triple-or-more fractions (mid B 0.24, 0.36, 0.40; early B 0.16, 0.32, 0.52; was 0.25, 0.40, 0.35 for all B), the rest split across three to seven stars in the B row's proportions. Mean star count 2.37 and 2.63 against 2.28.
- `tuning.MULTIPLICITY_BY_PRIMARY_MASS_BIN` names Table 13's rows (companion frequency, close frequency, single, binary and multiple fractions by mass bin), `tuning.B_SUBCLASS_LOWER_MASS_SOL` the B bins, and `multistar.companion_class` picks the class from the birth mass.
- `docs/design/multistar-and-compact-systems.md` 3.3 records the rows and the check against the table.

generation-output: changed

TODO items resolved: GEN.449 / Added: Table 13's close-binary frequency by mass (0.15 to 1.0) is named in `tuning` but not yet drawn from, for TEST.196
