**Headline:** The greenhouse closure that checks a planet's stored warming against a physical opacity is in place (GEN.443, greenhouse kernel).

### Changed

- New `planetgen.physics.greenhouse`: the Robinson and Catling grey-atmosphere closure (surface temperature from an infrared optical depth, and the inverse), the greenhouse factor from the effective temperature, and the check of a stored multiplier. It reproduces Earth (optical depth 1.82, factor 0.625) and a Venus-like world. Nothing calls it yet and no stored value changes; the composition side of the check waits on the infrared opacity data and the stored multiplier.
