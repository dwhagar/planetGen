**Headline:** The System Map now draws a planet's irregular moons on their own tilted, stretched orbits, many going backward, with an info panel and a link to each body's page.

### Changed

- GEN.1086 (slice 2): a planet that caught irregular moons has a scene of its own on the System Map (`planet-<id>-irregular`), opened from the planet when it has no regular moons, or by an "Irregular moons" button beside the moon system scene (and a "Moons" button back). The moons are tens of millions of kilometers out, so they have a scale of their own and the regular moons are not squeezed into a speck. Each body is drawn at its saved place on its saved orbit (eccentric, inclined, forward or backward); clicking one shows its type, radius, mass, distance, period, makeup and the window it was caught through, with a View button to its page.
- `GET /api/systems/<id>`: each planet has `irregular_moons` (an empty list when it keeps none); a new key, the API version stays. The planet's own page links to the irregular scene when it keeps only irregular moons.
- Not yet: the 3D view, and the clearance ledger taking the irregular zone as occupied space (the rest of GEN.1086).

generation-output: unchanged
