**Headline:** An asteroid field now ages with the in-universe clock: its radius grows as its pieces spread, and one that passes a light year is retired from the maps.

### Added
- **Bring-to-epoch for asteroid fields (GEN.1107).** The orbit update that moves a field through the galaxy now adds the years it advanced to the field's age (`birth_age_yr`), and the radius and the size digit of the class follow from the age by the formula that made them (`asteroid_field.aged`, `radius_ly_at`). A field whose radius passes one light year is retired, as Boss decided: the row stays, with its Object ID, members and any facility, and the map reads and the widest-radius reads leave it out. Retired is a formula of the row (`radius_ly > 1`), so no column was added and there is no schema change.

generation-output: unchanged

TODO items resolved: GEN.1107.
