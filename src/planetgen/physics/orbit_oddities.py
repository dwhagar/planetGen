# planetgen/physics/orbit_oddities.py

"""
The oddities of the orbit lever: a highly eccentric orbit (W3) and a polar or retrograde orbit round a hot star (W4) (GEN.1049, docs/design/systems/weird-worlds.md)
=================================================================================================================================================================

Two rows of `physics/data/oddities.json` act on the orbit of a planet (lever L2) and need no module of their own beyond the angle they set:

* **W3, the highly eccentric orbit.** The eccentricity of a planet is drawn by `orbit_elements.draw_planet_eccentricity` (the law of the Kepler multis for a crowded system, a
  Rayleigh law for a single planet, up to `MAX_ECCENTRICITY`) and worn down by the star's tide for the age of the star, so the tail of high eccentricities is already
  the generator's. The row names the planets that come out at `e` 0.5 or more (HD 80606 b is 0.93, HD 20782 b 0.95) so that the page says why their year swings by hundreds of kelvin;
  nothing on the body changes and no draw is made. Whether the tail holds the share the literature gives (about a tenth of the giants beyond 1 AU above `e` 0.5, WR1) is
  the census's to say (TEST.492), and the law is retuned there if it does not.
* **W4, the polar orbit.** About a third of hot Jupiters round hot stars are misaligned with the star's equator and the true angles of the ones measured sit in a polar peak between
  80 and 125 degrees (WR2). The generator keeps no star spin axis, so the angle of the orbit is the one stored for the plane: the effect sets `orbital_inclination_deg` to a
  draw from that peak on the row's own stream. The three rows are the three host bands of the rate (cool hosts under 6,250 K, 6,250 to 7,000 K, over 7,000 K) on a giant at less than
  10 stellar radii.

Stand-ins are in `docs/design/research/gap-filling-plan.md` section 6.92 (rows RC6 and RC7).
"""

from planetgen.physics import oddities

POLAR_PEAK_DEG = (80.0, 125.0)
"""tuple: The range the inclination of a misaligned hot Jupiter's orbit is drawn from, uniformly (WR2: 18 of the 19 misaligned systems with a measured true angle lie in it). Source tier T1, tweak level L1."""


@oddities.register("polar_orbit")
def polar_orbit(body, row):
    """The effect of oddity W4: the orbit's plane is turned to a polar one (`POLAR_PEAK_DEG`), from the row's own stream."""
    low, high = POLAR_PEAK_DEG
    body.orbital_inclination_deg = oddities.stream(f"{row.id}:angle", oddities.body_key(body)).uniform(low, high)
    return None


@oddities.register("high_eccentricity")
def high_eccentricity(body, row):
    """The effect of oddity W3: nothing on the body changes (the eccentricity is the generator's own; the row only names it)."""
    return None
