**Headline:** The dust model converts the solar mass through the unit registry, and the mass-conversion test passes again.

### Changed

- `physics/dust_model.py` (VIEW.43) multiplied and divided by `constants.SOLAR_MASS_TO_KG` by hand in two constants, so `test_mass_conversions.py` was red on main (TEST.639). Both now use `_MSUN_KG`, the registry's `converter("Msun", "kg")(1.0)`, which is the same number (checked equal to the constant), so the dust model's output is bit-identical.

generation-output: unchanged
