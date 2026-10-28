**Headline:** A greenhouse model on real line data is built (not yet used by the generator): doubling carbon dioxide warms Earth by about a kelvin with its vapour held, and Titan's haze brings it to its measured surface.

### Added
- **Downloaded infrared line data (GEN.715).** `update.sh` now also fetches the ExoMol correlated-k tables of water vapour, carbon dioxide and methane (about 1.3 GB) and the HITRAN collision-induced absorption of 11 gas pairs (about 0.2 GB), reduces them to 6.5 MB of band tables in the data cache and deletes the downloads. Generation requires them like the other datasets; `--disregard-missing-data` takes the labelled stand-ins instead.
- **A one-dimensional radiative-convective model (GEN.715, GEN.717, GEN.718).** `physics/radiative_convective.py` (the line opacity is `physics/line_opacity.py`) finds the surface temperature at which an air's adiabat, topped by a radiative stratosphere, sends the absorbed sunlight back to space, band by band and g point by g point, with Manabe and Wetherald's humidity, a water continuum, collision-induced pairs, and an optional share of the sunlight that a haze takes high up.

### Changed
- **The model is built but not wired.** `radiative_convective.warming` returns the surface temperature, the greenhouse factor G, the tropopause and the whole profile; `atmosphere_state.greenhouse_warming` keeps the grey closure of GEN.703 that the surface-air iteration calls several times per planet, because the new model takes one to several seconds a call. Checked, not tuned: a doubling of carbon dioxide with the vapour held warms 1.0 K (GEN.717), Titan with 60 percent of its sunlight absorbed in its haze sits at 93 K against the measured 94 K (GEN.718), clear-sky Earth at 276 K, Mars at 215 K, Venus at 711 K. Nothing the generator stores changes.

### Fixed
- **The moist adiabat took the wrong slope across the saturation corner.** A step that started unsaturated and crossed into saturation used the saturated slope for its whole length, so the top of the column came out up to 10 K too warm just below the corner and then fell as the surface warmed.
