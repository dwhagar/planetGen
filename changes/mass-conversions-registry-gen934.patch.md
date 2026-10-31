**Headline:** Mass conversions between kilograms, solar masses, Earth masses and Jupiter masses now go through the unit registry in 71 modules instead of being written by hand, with no change to any generated number.

### Changed
- **Hand-written mass conversions replaced (GEN.934).** A multiply or divide by `SOLAR_MASS_TO_KG`, `EARTH_MASS_TO_KG` or `JUPITER_MASS_TO_KG` is now a `unit_registry.converter` call (`_kg_to_msun`, `_mearth_to_kg` and the like) wherever the swap gives the same last bit. The pinned galaxy is identical.
- **Sites kept by hand.** 16 products such as `G * m * C` and `C / m` stay as written, because a converter call would regroup the multiplication; they are listed in `tests/test_mass_conversions.py`.

### Added
- **`tests/test_mass_conversions.py`**: fails on a new hand-written mass conversion or a stale entry in the kept list, and checks the six converters equal the hand multiply and divide to the last bit.

generation-output: unchanged
