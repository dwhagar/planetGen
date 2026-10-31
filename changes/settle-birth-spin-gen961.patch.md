**Headline:** A planet or giant made again at its star's final age now wears down the spin it was born with, instead of drawing a new one.

### Changed
- **The final-age settle (`refinish_rocky`, `settle_giant`) drew a fresh birth spin (GEN.961, GEN.959).** The tide history of the moons and the planet's spin starts from the birth spin kept on the body, so a body made again now runs `generate_spin(..., birth=_birth_spin)` and keeps its birth rate and tilt; a body restored without one draws as before. The settle stream draws one birth spin fewer per body, so the draws after it (precession angle, then the field, dose and colour that follow) move.
- **A moon the planet's own history replaced was left at the age it was made at (GEN.992, GEN.959).** The moons' history that runs with the planet's spin can swap a moon for a new body, and the settle then skipped it, so a moon kept an air, class and spin from an age several Gy too young (seed 38 of the final-age sweep, an 'Alramam IIIb' made at 1.7 Gy under a 13.8 Gy star). The moons are now read after the planet is made again.
- generation-output: changed
