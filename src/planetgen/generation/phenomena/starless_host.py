# generation/phenomena/starless_host.py

"""
The host a rogue planet is built against (GEN.1111, TEST.498): a rogue planet is the same `Planet` as a bound one, with no star.

The planet chain reads its star for the mass (the tide and the Hill sphere), the luminosity (the light it gets), the age and the donor's snow line
(whether it formed beyond it, so what it is made of). A `StarlessHost` is those four things for a body that left its star: the luminosity is zero, the
age is the rogue's own (`rogue_age.from_record`: the donor's age, less the time since the ejection), and everything else is the donor's own, so the
snow line and the formation the body was born with are the ones the architecture lost it from (`planet_architecture.Loss`). The body sits `FAR_AU`
from the host's mass, where its tide and its light are the ambient floor; the host's mass stays so that the chain's formulas are defined, and the
distance takes it out of every one of them.
"""

from planetgen.generation.planet import Planet

FAR_AU = 1.0e5
"""float: Where the body sits from the donor's mass, AU. At this distance the donor's tide on the body, and its light, are below anything the chain
resolves, and its orbital period is far longer than the age of the galaxy. STAND-IN (GEN.725): the chain's own placement for a body with no host."""


class StarlessHost:
    """
    The star a rogue planet is built against: the donor star it was ejected from, with no light and the rogue's own age.

    Args:
        donor: The star the body was ejected from (a `Star`, or a `BinaryStarProxy`).
        age_gy (float): The rogue's age, Gy (`rogue_age.from_record(donor.age, delay).age`).

    Every attribute the host does not set itself (the mass, the snow line, the spectral type, the metallicity) is read off the donor.
    """

    starless = True
    """bool: Read by the quantity builders (`mixture_subject`): this host gives the body an age and nothing else of a star's."""

    def __init__(self, donor, age_gy):
        self.donor = donor
        self.age = age_gy
        self.luminosity = 0.0

    def __getattr__(self, name):
        # Only called for what the host does not hold itself; `donor` missing (a copy before `__init__`) must not recurse.
        if name == "donor":
            raise AttributeError(name)
        return getattr(self.donor, name)


def build_ejected_planet(system_config, donor, age_gy, mass_kg, formation_au=None, moon_count=None):
    """
    A `Planet` of `mass_kg` built through the planet code with no star: the donor's own mass and snow line, no light, `age_gy` of age.

    Args:
        system_config (SystemConfig): The shared settings.
        donor: The star it was ejected from.
        age_gy (float): The rogue's age, Gy.
        mass_kg (float): The body's mass (one body of a `Loss`, whose `mass_kg` is the total of its `count`).
        formation_au (float, optional): Where it formed (`Loss.formation_au`); omitted, it is built where `FAR_AU` puts it.
        moon_count (int, optional): How many moons to make (`None` lets the chain draw them).

    Returns:
        Planet: The body. Its `star` is the `StarlessHost`.
    """
    host = StarlessHost(donor, age_gy)
    return Planet(system_config, host, host.habitable_zone, FAR_AU, mass=mass_kg, moon_count=moon_count, formation_distance_au=formation_au)
