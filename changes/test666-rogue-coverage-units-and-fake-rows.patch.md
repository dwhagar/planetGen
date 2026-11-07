**Headline:** A free-floating brown dwarf rebuilt from its saved form keeps the 25 temperature coverage fields (empty, as it was built), so the rebuilt object has every column its row stores.

### Changed

- **A rogue planet's saved form carries the temperature coverage fields (TEST.666).** `rogue_habitability.ROGUE_FIELDS` lists the 25 `temp_*` coverage fields (GEN.1265) that `rogue_planets` stores, but `ROGUE_UNITS`, which feeds `RoguePlanet.FIELD_UNITS` and with it the saved form, left them out, so `to_dict` dropped them and `from_dict` rebuilt a brown dwarf without the attributes (`AttributeError` on `temp_livable_share` in `test_a_rebuilt_brown_dwarf_keeps_its_empty_score`). They now carry the units of the same columns on `planets` (kelvin, degrees, days, share, flag). Values do not change (empty on a brown dwarf); the saved form of every rogue planet gains the 25 keys, which changes the structure digest of the golden battery.
- **The fake database of `test_object_rows` serves the irregular-moon rows (TEST.666).** The fake connection answered only the tables that existed when the test was written; it now answers `irregular_sets` and `irregular_moons` from the rows `store.insert_irregular_moons` writes, so a rebuilt planet is compared with the planet that was stored.

generation-output: changed
