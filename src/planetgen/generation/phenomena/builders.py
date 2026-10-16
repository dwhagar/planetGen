# planetgen/generation/phenomena/builders.py

"""
The sector-roll builders of the phenomena registry (GEN.250): one function
per kind that a sector rolls, each taking `(config, galactic_center_dist_ly)`
and returning the new object. `registry.Kind.builder` names them; only black
holes and neutron stars use the distance (their Hill sphere and galactic
orbit).
"""

from planetgen.generation.phenomena.asteroid_field import AsteroidField
from planetgen.generation.phenomena.compact_remnant import BlackHole, NeutronStar
from planetgen.generation.phenomena.nebula import Nebula
from planetgen.generation.phenomena.quasar import Quasar
from planetgen.generation.phenomena.rogue import InterstellarComet, RoguePlanet
from planetgen.generation.phenomena.supernova_remnant import SupernovaRemnant


def build_black_hole(config, galactic_center_dist_ly):
    return BlackHole(config, galactic_center_dist_ly=galactic_center_dist_ly)


def build_neutron_star(config, galactic_center_dist_ly):
    return NeutronStar(config, galactic_center_dist_ly=galactic_center_dist_ly)


def build_nebula(config, galactic_center_dist_ly):
    return Nebula(config)


def build_planetary_nebula(config, galactic_center_dist_ly):
    return Nebula(config, nebula_type="planetary")


def build_molecular_cloud(config, galactic_center_dist_ly):
    return Nebula(config, nebula_type="dark")


def build_supernova_remnant(config, galactic_center_dist_ly):
    return SupernovaRemnant(config)


def build_rogue_planet(config, galactic_center_dist_ly):
    return RoguePlanet(config)


def build_brown_dwarf(config, galactic_center_dist_ly):
    return RoguePlanet(config, mass_bin="brown-dwarf")


def build_comet(config, galactic_center_dist_ly):
    return InterstellarComet(config)


def build_asteroid_field(config, galactic_center_dist_ly):
    return AsteroidField(config)


def build_quasar(config, galactic_center_dist_ly):
    return Quasar(config)
