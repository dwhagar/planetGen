# stellarObjects/__init__.py

"""
Stellar Objects Package
=======================

This package contains the core classes for generating procedural star systems,
including stars, planets, and moons.

This package exports the following classes:
    - Star: Represents a single star and its properties.
    - Planet: Represents a single planet or moon and its properties.
    - StarSystem: Represents a full star system, including a central star and a list of planets.
    - SpaceSector: Represents a named collection of star systems laid out in 3D space.
    - SectorSystemEntry: One star system's placement (and generation recipe) within a SpaceSector.
"""

from .planetData import Planet
from .spaceSector import SectorSystemEntry, SpaceSector
from .starData import Star
from .systemData import StarSystem

__all__ = ['Planet', 'SectorSystemEntry', 'SpaceSector', 'Star', 'StarSystem']