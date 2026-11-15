**Headline:** The oxygen a world's life makes now settles at a pressure that does not grow with the world's gravity, so heavy worlds stop piling up at the 35.5 kPa fire limit.

### Changed

- GEN.1733 (every galaxy): the oxidative weathering that balances the oxygen source now goes with the oxygen's partial pressure to the power 1.0, linear (`biogenic_oxygen.OXIDATIVE_WEATHERING_EXPONENT`; the sources run from the square root of the COPSE model [R5280] and the fit of [R5281] to a linear form, and the linear end gives the habitable classes' oxygen window of 16 to 35.5 kPa a width of a factor 2.2 in the surplus instead of 1.5), where it went with the moles of oxygen in the air. A world's balance pressure is Earth's 21.2 kPa times its surplus per area over Earth's (`balance_pa`); the fire limit still caps it. The measure (reports/class-m-tiers-build/oxygen-35kpa-gen1733-2026-11-15.md): 46 of 64 oxygen-making ocean worlds of the Abundant run sat at the limit, all but 4 of them over 1.2 g, because the same oxygen per area gives a pressure in proportion to gravity.
- A world with Earth's surplus per area now holds 21.2 kPa at any gravity (it held about 22 kPa per g). A world whose surplus is below Earth's holds less oxygen than before and one above it holds more until the fire limit; a refuge world's tiny surplus (GEN.1723) gives a tiny pressure.
- Measured on 400 Abundant systems (reports/class-m-tiers-build/meter-gen1733-oxygen-fix-2026-11-15.md): of 60 oxygen-making sea worlds, inside the 16 to 35.4 kPa window 24 (11 before the change), at the 35.5 kPa limit 29 (41 before), under 16 kPa 7 (8 before).

generation-output: changed
generation-affects: fill
