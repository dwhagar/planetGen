### Changed
- **Planets and moons are named after their star.** Only the system draws
  a generated name now. Planets are numbered in orbit order with roman
  numerals (`Voranthis I`, `Voranthis II`), moons add a letter
  (`Voranthis IIa`, `Voranthis IIb`), and asteroid belts take no number.
  A binary's two stars each put their own word after the system name
  (`Voranthis Kelmoor`, `Voranthis Ostra`) instead of `Voranthis` and
  `Voranthis B`; a close pair's planets use the system name and a wide
  pair's use their own star's (`Voranthis Kelmoor I`). See
  `src/stellarObjects/bodyNames.py`.
- **Planet and moon names are no longer searched for duplicates.** They
  derive from the system name, which is already unique, so saving a
  system skips a registry lookup and write per planet and moon, and the
  "Kin"/"Ami" companion suffixes are gone. When a system is renamed
  (`PATCH /api/systems/<id>`, or the Alpha/Beta and "Little" decorations
  a name collision applies), every star, planet and moon still carrying
  its name follows. Schema v34 drops `body_name_registry`; run
  `update.sh`. Existing rows keep their old names until the galaxy is
  regenerated.
- **`PATCH /api/systems/<id>` refuses a name already in use** by any
  sector, system, star, planet or moon, with a `409`.

### Added
- **Rename stars, planets and moons through the API.**
  `PATCH /api/stars/<id>`, `/api/planets/<id>` and `/api/moons/<id>` take
  `{"name": str}` and need an admin login, like the other writes.
  Renaming a single star renames its system; renaming a binary's star
  carries the planets and moons named after it. See `docs/api.md`,
  "Renaming".
