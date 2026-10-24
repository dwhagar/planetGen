# tests/member_systems.py

"""
A system of three stars for the tests of GEN.630, built by hand: a binary the generator
made with one more star orbiting the pair as a whole. The tree is two nodes, node 0 the pair (as a pair is stored,
now the slot-0 child of the new root) and node 1 the root, with the pair in slot 0 and the new star in slot 1.
"""

import dataclasses

import pytest

from planetgen.db import orbit_trees
from planetgen.generation.config import SystemConfig
from planetgen.generation.star import Star
from planetgen.generation.system import StarSystem
from planetgen.physics import constants
from planetgen.physics.system_tree import OrbitNode
from planetgen.util import draw

EPOCH_UNIX = 1_760_000_000


def binary(wide, planets=False):
    """A seeded binary system of the asked kind and the `SystemConfig` that made it."""
    for _ in range(40):
        cfg = SystemConfig()
        cfg.PLANETS = planets
        cfg.BINARY_SYSTEM = True
        cfg.WIDE_BINARY = wide
        system = StarSystem(system_config=cfg)
        if system.binary_type == ("wide" if wide else "close"):
            return system, cfg
    pytest.fail("no binary of the asked kind in 40 draws")


def with_member(system, cfg, mass_sol=0.4, separation_au=400.0, seed=630):
    """`system` with a third star orbiting its pair; returns it."""
    pair = system.primary_star
    with draw.bound(seed):
        member = Star(cfg, name=f"{system.name} C", initial_mass_sol=mass_sol, age_gy=pair.age,
                      galactic_center_dist_ly=getattr(pair, "galactic_center_dist_ly", None),
                      galactic_orbital_phase_deg=getattr(pair, "galactic_orbital_phase_deg", None))
    pair_node = dataclasses.replace(orbit_trees.binary_nodes(system, EPOCH_UNIX)[0], parent_node_no=1, slot=0, epoch_unix=float(EPOCH_UNIX))
    total = sum(star.mass for star in system.stars) / constants.SOLAR_MASS_TO_KG + member.mass / constants.SOLAR_MASS_TO_KG
    root = OrbitNode(node_no=1, parent_node_no=None, slot=None, mass_solar=total,
                     separation_km=separation_au * constants.AU_TO_KM, eccentricity=0.3, inclination_deg=20.0,
                     ascending_node_deg=40.0, arg_periapsis_deg=10.0, mean_anomaly_deg=75.0,
                     epoch_unix=float(EPOCH_UNIX), period_years=4000.0, kind="wide", interaction="detached", flags=0)
    system.orbit_nodes = (pair_node, root)
    system.member_stars = (member,)
    system.stars = list(system.stars) + [member]
    system.star_slots = ((0, 0), (0, 1), (1, 1))
    return system
