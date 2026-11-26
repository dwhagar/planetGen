"""
A copy of a planet or moon to try an air on (moved here from `generation.life` so the physics of the life's air can use it, Train 62).

The copy has its own attributes, its own dicts, lists and sets, and its own build record, so storing an air on it leaves the body as it was. The star, the moons and every container of
bodies are shared with the body (GEN.1734): each body holds the siblings of its host (`_tide_siblings`, the moons it was born with), so a deep copy of one walks every body of its system in
one chain and runs the stack out in a system of a few hundred bodies.
"""

import copy


def trial_copy(planet):
    """A copy of `planet` to try an air on: its attributes copied, with the containers (the air, the spectra, the layers) copied too, so building the air on
    it leaves the body as it was. The star and the moons are shared (nothing in the air pass changes them). Its build record is copied the same way."""
    clone = copy.copy(planet)       # the clone has a record of its own (`Planet.__copy__`)
    _copy_containers(planet, clone, skip=("star", "moons"))
    _copy_containers(planet.build, clone.build)
    return clone


def _copy_containers(source, clone, skip=()):
    """Gives `clone` its own copy of every dict, list and set `source` holds, so a pass on the clone leaves the source's as they were."""
    for name, value in vars(source).items():
        if name not in skip and isinstance(value, (dict, list, set)):
            # GEN.1734: a list of bodies (the moons the planet was born with, `build.born_moons`) is shared like `moons`: each body holds the siblings of its host, so a deep copy of
            # one walks every body of the system in one chain (a recursion of 160 deep systems) and copies bodies the air pass never changes
            setattr(clone, name, copy.copy(value) if _holds_bodies(value) else copy.deepcopy(value))


def _holds_bodies(container):
    """Whether a list, set or dict has a planet, moon or other body in it (an object with a `body_type`)."""
    items = container.values() if isinstance(container, dict) else container
    return any(getattr(item, "body_type", None) is not None for item in items)
