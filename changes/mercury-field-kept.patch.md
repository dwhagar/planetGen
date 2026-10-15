### Fixed

- **Mercury-like planets keep their magnetic field (GEN.213).** A small, iron-rich planet (0.05 to 0.1 Earth masses, like Mercury) was always given no field, because the dynamo-lifetime test ran before the rule that gives it Mercury's 3e19 A m2 dipole. The rule now runs first, so such a planet has a field at any age. Fields stored before this release change when a planet is regenerated.
