**Headline:** a spelling fix only: "centre" in the obliquity-variation code and its test is now "center".

### Changed

- `test_american_spelling` was red on main because the obliquity-variation module (GEN.1454 part 1) and its test used "centre" (including the field `cos_centre`). Renamed with the repo's own sweep (`scripts/american_spelling_sweep.py apply`); nothing else in the tree was found. No change in behavior.
