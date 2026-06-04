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
"""

from .planetData import Planet
from .starData import Star
from .systemData import StarSystem

__all__ = ['Planet', 'Star', 'StarSystem']
