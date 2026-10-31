**Headline:** The last hand-written mass conversions that went in since the mass-conversion test was written now use the unit registry's converters, and the output is unchanged.

### Changed
- TEST.510: `describe/sheet.py` (2), `generation/member_record.py`, `generation/phenomena/rogue.py` (3), `physics/belt_mass.py` (3), `physics/binding_energy.py`, `physics/planet_architecture.py` and `tuning.py` (the giant core's Jupiter-mass threshold) call `unit_registry.converter` instead of multiplying or dividing by a mass constant. Each converter does the same one multiply or divide, so every number is the same bit for bit (the constant on the left of a divide in `rogue.py` is `_mearth_to_kg(1.0)`, which is `1.0 * C`; the `tuning.py` product is `(x * E) / J` as it was). `physics/belt_mass.py` comes off the kept-sites list in `tests/test_mass_conversions.py`; nothing was added to it and the test is not loosened.

generation-output: unchanged
