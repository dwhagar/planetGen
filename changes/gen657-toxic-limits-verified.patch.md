**Headline:** The toxic gas limits behind the survival equipment line are now read from the NIOSH Pocket Guide instead of recalled, and a test holds the gas filter grade to the physics.

### Changed
- **The toxic rows are verified (GEN.657).** Each of the ten rows in `protection_limits.TOXIC_GASES` now cites its NIOSH Pocket Guide entry; the nitrogen dioxide IDLH is 13 ppm (NIOSH 2017-202), where the design draft had recalled 20 ppm. No world reads nitrogen dioxide yet, so no generated value changes. `tests/test_equipment.py` holds each published IDLH and checks that the gas filter grade follows the ratio to the long-term limit, that a gas past its IDLH asks for own air, and that a corrosive one asks for the chemical-barrier suit.

TODO items resolved: none (GEN.657 stays open: reading the stored ammonia, hydrogen cyanide and the other toxic gases into the equipment rules changes stored values and waits for the generation lock). Added: none.
