**Headline:** The tuning targets now state which wrong values they can refuse, and the thorium rows share one scale.

### Changed
- **The thick-disc thorium row o2 is on the solar-twin scale of row o1 (GEN.664).** Its centre moves from -0.063 to 0.061 dex (the Hypatia thick-minus-thin offset of -0.019 added to +0.08), its target with it. The recipes are not read by the generator yet, so no generated value changes.

### Added
- **A power check for every tuning target (GEN.665).** `tuning_harness.power` tells whether a factor-2 error is refused; two weak targets are named, and the white-dwarf redox and mass-radius notes record why the width 0.22 and the weak cut stay.
