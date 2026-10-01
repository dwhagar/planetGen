### Added

- **Every bright star is placed before its sector is filled (schema
  v43).** `generate.py plan` now ends by drawing every star of 500 solar
  luminosities or more across the whole galaxy and storing each at a
  fixed point in its sector, in a new `bright_stars` table, while the
  sectors themselves stay unfilled. Filling a sector builds a full system
  around each of its bright stars first and draws the rest from dimmer
  stars, so a sector's expected count is unchanged. New plan options:
  `--bright-star-min-luminosity`, `--no-bright-stars`,
  `--bright-stars-only` and `--force`.
- **Star ages follow where a sector sits.** Each system in a galaxy
  sector draws its star from the young, intermediate, old or bulge
  population in proportion to their density there, so O and B stars and
  supergiants gather in the spiral arms near the plane.
