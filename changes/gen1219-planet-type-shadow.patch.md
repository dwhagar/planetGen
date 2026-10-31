**Headline:** A read-only report now sorts every generated planet and moon into terrestrial, transitional or gas by the share of its mass in an envelope.

### Added
- **The body type function and its shadow report (GEN.1219).** `planetgen.physics.planet_type` gives the type from `envelope_mass_fraction` (terrestrial under 0.5 percent, transitional to 50 percent, gas above), and `src/tests/planet_type_cli.py` prints how a generated sample splits, the bodies near each cut and the bodies whose radius disagrees with their type. Nothing stored changes; the stored column comes with the next epoch.
