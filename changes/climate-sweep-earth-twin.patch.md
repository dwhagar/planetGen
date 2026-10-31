**Headline:** The climate sensitivity test now measures the real Earth, not a random hot planet.

### Changed
- **TEST.509: the climate sweep twin is the Earth.** `test_climate_sweep` used a drawn body at 313 K with a hydrogen and carbon dioxide air, whose response said nothing about the Earth's. The twin now holds the Earth's dry air (nitrogen, oxygen, argon, 420 ppm carbon dioxide at 1 bar), gravity, radius and ocean through a stand-in for the mixture draw, with the absorbed flux set so the fixed point stands at the measured 288 K (the model is 6.5 K cold at the Earth's own flux, having no cloud). The response there is 0.40 K per W/m2, under the published lower bound of 0.54; the strict expected failure now names the budget (Planck with water vapour and lapse rate 0.52 by IPCC AR6 table 7.10, 0.58 with the surface albedo and cloud feedbacks of GEN.972 and GEN.974). Test only, no world changes.
- generation-output: unchanged
