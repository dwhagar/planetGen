**Headline:** A Class M abundance tier now draws the carbon cycle of its Earth-shaped planets from the part of the published ranges that keeps their CO2 where life can use it.

### Added
- **The carbon-cycle dial (GEN.1720; Boss, 11:58 PM EST 2026-11-14: "recommendation approved").** A sixth dial, `cycle_share` (0.5 on Sparse, 1 on the other tiers): the share of rocky planets in the ecosphere that draw the exponent of the heat flow in the volcanic outgassing (published 1 to 2) from 1 to 1.25 and the factor of the weathering law's thermodynamic limit (published a tenth to ten) from 0.1 to 1. The balance of the carbon cycle goes as the outgassing to about the fifth power, so over the whole ranges most Earth-shaped planets hold under 2 Pa or over 100 kPa of CO2 and few the 2 Pa to 1 kPa the life classes ask for; in the narrower ranges a planet with a heat flow of 0.6 to 2 times the Earth's holds CO2 inside that band at an Earth-like temperature. The dial moves the draw of two constants and conditions nothing on the CO2 that follows. One draw per candidate under a tier. Registry row 6.289:AT8; `tuning.ABUNDANCE_TIER_DIALS` rows are now six values.
- **Attribution:** `reports/gen1720-co2-attribution-2026-11-14.md` (GEN.1713 is a minor part; the snow line fix and the reshuffled draws of GEN.1646 moved the population and the warm end).

Science is untouched: its dial is 0, no draw is made, no body holds a new field. The nudge is rolled on a stream of its own (seeded from the run stream's state, which it does not move), so a tier's body draws the two constants from the ambient stream as Science does and nothing else of it (the air recipe, the sulfur lifetimes, the seafloor law) is reshuffled (the first form took one ambient draw and moved them; found by the Class M tiers meter, Edabihion VIII).

Reuse: `continental_weathering.draw_cycle` and `mixture_base.bracket` (the constants are drawn and recorded as they were), `abundance_tier.candidate`, the dial table and its helpers of GEN.1704.

generation-output: changed (a tier other than Science only)
generation-affects: none
