### Added

- The mass limit (GEN.183) is now picked from the presets 8, 10, 12, 14, 16, 18 and 20 solar masses, 20 by default: a slider in the Generate page's new galaxy and plan forms, and `--phenomenon-min-mass` on the command line, which refuses any other value. Every star, neutron star and black hole at or above it is placed across the whole galaxy; lighter ones are drawn when their sector is made.

### Changed

- A scatter or a phenomena-only re-scatter run without `--phenomenon-min-mass` now keeps the mass limit already stored with the galaxy instead of going back to 20.
