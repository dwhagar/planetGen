**Headline:** A main-sequence star now takes its brightness, colour and size from the MIST isochrones at its mass, age and composition, and it brightens as it ages.

### Changed

- **One main-sequence table (GEN.652).** The generator's own mass-luminosity pieces and mass-radius power law are gone. A main-sequence star's luminosity, effective temperature and radius are the MIST row of its initial mass, its age and its birth [Fe/H] (`physics/main_sequence.py`), so a star is dimmer when young and brighter near the end of its main sequence, and a metal-poor star is hotter than a metal-rich one of the same mass. The pre-main-sequence tracks contract onto the star's own MIST arrival row, a subgiant starts from its last row, and a supergiant brightens from it. The sampling tables for bright stars find the age a star first clears a brightness by reading the same table. A mass below 0.1 solar masses, or an age past MIST's main sequence, takes the row at the table's edge or a power-law extension (stand-in).
- The composition note of a main-sequence star reports the same row the star's record reads.

generation-output: changed
