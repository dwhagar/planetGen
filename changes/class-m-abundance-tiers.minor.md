**Headline:** A galaxy now carries a Class M abundance setting: Science (the default, which changes nothing), Sparse, Scattered, Plentiful or Abundant.

### Added
- **Class M abundance tiers (Boss, 11:19 PM EST 2026-11-13).** Known science cannot say how common an Earth is, and the generator's own laws give under one Class M world in ten thousand Sun-like systems. A galaxy's maker now chooses how far past the science to go: the `--class-m` option of `planetgen plan` and the Class M abundance field of the Generate page (Census dial, beside Rogue planets), stored with the galaxy as the Rogues level is (no migration). The targets are 1 Class M world per 1,000, 100, 30 and 10 Sun-like systems for Sparse, Scattered, Plentiful and Abundant.
- **Three dials and nothing else** (`generation/abundance_tier.py`, `tuning.ABUNDANCE_TIER_DIALS`): the share of rocky planets in the ecosphere given Earth's water class, the share of them that captured no hydrogen envelope (the late core of GEN.1323), and the share of systems given the habitable-world request (GEN.1130). No class limit moves, no class is named and no count is set: every class still comes from the finished conditions.
- **The meter.** `scripts/class_census.py` draws Sun-like systems under a tier and prints the funnel from candidates to Class M and Blue with the rate per 1,000 systems beside the tier's target.
- **Design:** `docs/design/habitability/class-m-abundance-tiers.md`. The dial values are first-pass stand-ins (tier T9, tweak L4; registry rows 6.288:AT1 to AT4) and are set again by census after GEN.1694 to GEN.1698.

Reuse: the Rogues level's lever row, storage and page field (`galaxy/dials.py`, `galaxy_shape_pins`), `FillContext.apply_population`, `Planet(late_core=True)`, `prevalence.resolve`, `util/request_odds.py`, `physics/delivered_water.py` and `scripts/class_census.py`.

generation-output: unchanged
generation-affects: none
