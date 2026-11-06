**Headline:** A giant's stored albedo and greenhouse pair is tested to give its temperature back at every mass, and the scale-height inversion says plainly that it does not apply to a giant.

### Fixed
- `climate_draws.from_stored` no longer runs the barometric law on a giant's row (every generated giant came out with an albedo under 0, or as an "edited" row): it answers that a giant's scale height is not the no-air temperature's, and the row is backfilled from its stored temperature against the starlight alone, which reproduces it.
- A test draws giants of ten masses from 5 to 4,000 Earth masses and checks that the stored pair (with the giant's own heat inside its greenhouse quantity, GEN.1030) and the backfill both reproduce the temperature to a part in a billion.

generation-output: unchanged
