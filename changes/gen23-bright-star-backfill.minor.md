### Added
- **Bright stars around every generated sector (GEN.23).** As soon as a
  galaxy sector is generated, anywhere and by any route (the command
  line, the Generate page, the Sector page, or a Galaxy Map visit), every
  sector block (3x3x3 sectors) within 100 ly of it gets every star from
  100 L_sun up to what was already placed there. Each block remembers how
  dim it has gone (new `bright_star_blocks` table, schema v48), so a
  block is only drawn once, and sectors that are already filled are
  never touched. The new stars show on the Galaxy Map like the plan's
  bright stars, and later sectors in those blocks build their systems
  around them. A later `generate.py plan --bright-stars-down-to` band
  leaves those blocks out and tops each one up below its own level, so
  no star is ever drawn twice.

### Changed
- **The default generate-around sphere is 12 pc (about 39 ly), not
  100 ly.** "About 10 pc, rounded up" to the next whole 4 pc sector. It
  applies to `generate.py galaxy`'s random start, the Sector page's
  neighborhood button, `POST /api/sectors/<id>/generate-neighborhood`
  with no radius, and the Galaxy Map's neighborhood dialog. A radius you
  give is still used as is.

Run `update.sh` (or `update.ps1`) after updating: it migrates the
database to schema v48, adding one empty table. No regeneration is
needed; the backfill starts with the next sector generated.
