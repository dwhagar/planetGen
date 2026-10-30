**Headline:** The generator can now read an alpha-enhanced star off a solar-scaled track and tie a disc star's metallicity to its age.

### Added
- **Alpha-enhanced equivalent metallicity and the age-metallicity relation (GEN.744, lock-free step).** `physics/alpha_tracks.py` gives the [Fe/H] a high-alpha star reads on the solar-scaled MIST grid (Salaris et al. 1993: about +0.29 dex at [alpha/Fe] +0.4) and `physics/age_metallicity.py` the disc's mean [Fe/H] by age on the thin and thick branches (0.025 and 0.15 dex per Gyr, Haywood et al. 2013). Nothing in generation calls them yet; the wiring rides the star epoch. Stand-ins AL1, AR1 to AR3 registered.
