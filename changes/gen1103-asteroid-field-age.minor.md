**Headline:** An asteroid field now has an age and a speed for its pieces, and its radius follows from the two instead of being drawn, so every field can be brought to another epoch by the same formula.

### Added
- **Age, speed, shape, origin, history and largest piece of a field (GEN.1103).** Schema 221 adds six columns to `asteroid_fields` (`rms_speed_kms`, `birth_age_yr`, `shape`, `origin`, `history`, `largest_body_kg`), with a migration that adds each one only if it is missing. A field the Rogues top-up plans keeps its planned speed, age, shape, origin and largest piece, and gets the sentence of where it seems to have come from (no system is claimed). The field's page shows them.

### Changed
- **A field's radius is a result of its age and speed (GEN.1103).** A field made by hand takes a speed log-uniform from 0.1 to 15 km/s and an age uniform from the age of the smallest field to the age at which a field retires, and its rms radius is what `field_dispersal.rms_radius_m` gives for the two, so the radius a field is stored with is the radius its age gives. The size digit of a top-up field's class now follows its planned radius (it followed the radius drawn before the plan replaced it). The speed range is stand-in row AF1 in section 6.202, and the registry text of the asteroid-field kind says a hand-made field is given a speed and an age.

generation-output: changed

TODO items resolved: none (GEN.1103 part 1: the columns, the radius from the age and the registry text; the field density leaving 0 and the calls of the field mass, rate and composition modules stay open under GEN.1103).
