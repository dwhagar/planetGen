**Headline:** The stellar mass function now lives in one module, so the next change to it is made in one place.

### Changed

- **One initial mass function module (GEN.926).** The mass function was written out in six places (the draw and the count in the stellar evolution module, the moments in the system mass module, the stars per solar mass in the H II module, the mass a population returns in the disk history, the mass grid in the star population module). It is now `physics/imf.py`: the draw of a star's mass, the number or mass of stars in a range, the density and the segments a table is built on. Nothing generated changes (the golden output is the same); the old names are gone.
