**Headline:** The sky backfill draws fewer faint stars: every mass ring now stops at the distance where one in 200 of its stars still shows, instead of the rings of the Sun's mass and up all reaching the sphere of a 50-luminosity star, so about 1.4 times fewer stars are drawn and about 1 percent of the stars that would show at magnitude 7.5 are not drawn (GEN.1486 slice 1).

### Changed
- **One reach rule for every ring (GEN.1486).** `sky_tiers.build` gives each mass ring the reach `ring_reach_pc` gives it (the farthest distance at which `RING_MIN_YIELD` of the ring's stars still show, never past a 50 Lsun star), where the rings of 1 solar mass and up had no rule and reached the sphere of a 50 Lsun star. `RING_MIN_YIELD` is 0.005 (it was 0.01 for the rings under 1 solar mass): Boss's choice, "skip the faintest, lose 1%", of the measured 0.001 to 0.01 (reports/gen1486-ring-drawn-vs-shown-2026-11-12.md). Stand-in SKT5 in GEN.725 is updated (tier P, tweak level L3).
- **A ring with no star under the cap is left out.** The ring of 8 solar masses and up holds no living star dimmer than 50 Lsun (all of them belong to the luminosity ledger), so its sample is empty (`ring_luminosities` no longer raises) and it reaches nowhere. The 5 to 8 ring holds few dim stars and now reaches about 42 pc instead of 238.
- At the observer 8,000 pc from the centre with a 7.5 sky, the rings visit 1.4 million cells where they visited 3.5 million (2.5 times fewer); the rings of 1 to 2 and 2 to 5 solar masses reach about 207 and 229 pc instead of 238.

generation-output: changed
generation-affects: scatter
