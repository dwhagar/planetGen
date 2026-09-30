### Changed

- **Schema v37: interstellar objects at real-world rates.** Each sector
  now draws its phenomena per star from the research densities in
  `docs/design/interstellar-object-rates.md`
  (`program_constants.PHENOMENON_DENSITY_PC3`, with a
  `PHENOMENON_RATE_SCALE` dial per type). Isolated asteroid fields are no
  longer generated (they disperse), and planetary nebulae now appear.
- Rogue planets are drawn from four mass bins (terrestrial, sub-Neptune,
  Saturn-class, Jupiter-mass), so terrestrial rogues are now the most
  common kind instead of almost never appearing.

### Added

- Free-floating brown dwarfs, stored as rogue planets with
  `rogue_planets.mass_bin = 'brown-dwarf'`; the migration fills `mass_bin`
  for existing rogues from their mass.
- Runaway and hypervelocity stars: `star_systems.runaway_class` and
  `runaway_speed_kms`, flagged on ordinary systems (hypervelocity stars
  grow rarer with distance from the galactic center).
- `queryDb.sector_detail` reports each sector's star count and its
  estimated count of interstellar comets and planetesimals.
