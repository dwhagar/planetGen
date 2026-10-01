### Changed

- A new galaxy's galaxy-wide bright-star scatter now places every star of 1,000 solar luminosities and up (was 500). A galaxy already scattered keeps its level and its stars.
- The bright-star backfill around each generated sector is tiered by distance (GEN.30): down to 100 solar luminosities within 10 ly, 250 within 25 ly, 500 within 50 ly and 750 out to 100 ly. A block takes the tier of its nearest sector, and a block a nearer sector reaches later is topped up with only the band it lacks.
- The Generate page has a "Bright stars from (solar luminosities)" field on New galaxy, Plan and Rebuild the bright stars, passed to the scatter as `--bright-star-min-luminosity`.
