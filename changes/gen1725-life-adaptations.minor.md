**Headline:** Life on a cold or hot world can now carry a cold or heat adaptation: a table of how life widens where it can live (ice, rock, crust, salt, hot springs) is in the data, and the two new traits are drawn and shown.

### Added
- GEN.1725: `physics/data/life_adaptations.json` and `physics/life_adaptations.py` hold 36 adaptations of life to cold, heat, dryness, salt and dim light from the GEN.1724 research (the climate, the life that carries it, the axis and the limit at the tolerant and extreme grades, the stage the life must have reached, the habitat it opens, source, tier and tweak). At the extreme grade for an oxygen-making class the table reaches the measured edges: 253 K (lichen, -20 C), 346 K (cyanobacteria), 0.04 umol photons m-2 s-1 under the Arctic pack, 75 hours of liquid water a year, a root 68 m deep. Stand-ins AD1 to AD7 in section 6.294 of the gap plan; design note `docs/design/habitability/life-adaptations.md`.
- Two adaptation traits, cold and heat, join GEN.876's three: a body whose mean surface temperature is past the ordinary limit of an oxygen-making life (under 273 K or over 313 K) draws the trait with the same law as the others (chance 0.8 times the score, extreme with chance 0.5 times the score squared), stores it in `life_traits` (`cold:3`) and the object page shows "cold (extreme)". They come last in the draw order, so every body's radiation, sulfide and carbon dioxide traits are what they were.

### Changed
- Nothing reads the new traits yet: the scores of the habitability index and the life classes are unchanged. The refuge gate of GEN.1723 will read `adaptation_traits.limit`; widening the class 1 gate waits for it.

generation-output: changed
generation-affects: fill
