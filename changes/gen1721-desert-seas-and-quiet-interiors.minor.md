**Headline:** A Class M abundance tier now gives its desert worlds small living seas and gives half of its candidates quiet volcanoes.

### Added
- **Desert seas (GEN.1721, second slice; Boss, card of 9:49 PM EST 2026-11-14: "both, small seas first").** Under a tier other than Science, a desert world (the dry part of the water dial) draws its water inside the part of the dry class that leaves a few small seas (`tuning.ABUNDANCE_DESERT_WATER_RANGE`), and its small sea is an ordinary salt sea that life can start in, not the concentrated brine so little water makes on any other world (`abundance_tier.desert_sea`, `hydrosphere.salts_concentrated`). On every other world a brine seep is what it was. Registry row 6.289:AT7.
- **Quiet interiors (Boss, 9:43 PM EST 2026-11-14: "tweak the volcanic activity a bit").** A fifth dial, `quiet_share` (0 on Sparse, 0.5 on the other tiers): the share of the rocky planets in the ecosphere drawn with a quiet interior, whose volcanic gases take a third of the oxygen a body of the same heat loses (`tuning.ABUNDANCE_QUIET_SINK_FACTOR`, the low end of Earth's own published range). One draw per candidate under a tier. Registry row 6.289:AT6.
- **Design:** `docs/design/habitability/class-m-abundance-tiers.md` sections 3, 6 and 8.

Science is untouched: every dial is 0 there, no draw is made, no planet holds a new field and the golden systems are byte for byte the same (no generator epoch, no migration). Desert worlds stay rare until the land plants of GEN.1729: the oxygen law counts the open sea alone, so a world with under 20 percent sea reaches breathable oxygen only when it is warm, bright and quiet inside.

Reuse: the envelope dial's draw point in `Planet._draw_body`, the water dial's class and pin (`delivered_water.dialed_class`, `body._water_dialed`), `hydrosphere.ocean_class` and `biogenic_oxygen.reductant_flux` unchanged for every other body; the oxygen source law is not touched.

generation-affects: none
