**Headline:** A planet's or moon's strength now follows how hot it was when its own weight shaped it: bodies that took their shape warm (heated by the way they were assembled, by radioactive decay or by tides) are weaker and turn round at smaller sizes, and cold icy bodies are stronger than before.

### Changed

- GEN.1201: planets and moons store `shape_temperature_k` (migration 0244, `planets` and `moons`, DOUBLE, K, NULL for a giant and for a row made before): the temperature `physics/shape_temperature.py` works out from the surface temperature, the binding energy of the body's assembly, its stored radioactive content (K, Th, U and 26Al) and its tide's flux, capped at the rock's liquidus. `body_strength.values` gives it to the roundness engine, so `roundness_ratio` and `yield_strength_pa` are no longer the cold strength: ice creeps near its melting point (a tenth of the cold strength) and is stronger below the reference 158 K, and dry olivine creeps below its cold strength from about 1,000 K.
- Effect on a seeded sample of 120 systems (360 planets and moons): 5 bodies change between round and not round (1 moon becomes round, 4 cold ice-rich moons around far-out planets become not round); planets stay round, with a yield strength far below the cold rock (a planet is held at its liquidus). The strength column of a planet is now the strength at the temperature it took its shape at, not the cold strength.
- Known limits, in `gap-filling-plan.md` section 6.162: the module counts no latent heat or convective loss, so with the canonical 26Al every rocky body over about 100 km reaches the ceiling; comets, rogue planets and belt members keep the cold strength.

generation-output: changed
generation-affects: fill
