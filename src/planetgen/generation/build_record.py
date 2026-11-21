# planetgen/generation/build_record.py

"""
The build record (GEN.1651)
===========================

What a planet or moon carries from one pass of its build to a later one and never stores: the numbers a pass
leaves for the pass that runs the body's history again (the birth spin, the birth eccentricity), the draws an air
began with, the light an air was made on, the capture of an envelope. Each of these was once set as a private
attribute on the body (`planet.build.birth_spin`) and read back with `getattr(planet, "build").birth_spin`; they are
fields of one `BuildRecord` now, `body.build`, with a default for every field, so a reader writes `body.build.birth_spin`
and a body restored from the database (which has built nothing) reads the same defaults a new one does.

The record is not stored (`SERIALIZABLE_FIELDS` names what is), is copied with the body (`copy.copy`, `copy.deepcopy`,
`life._trial_copy`) and leaves `co2_sink` out of a pickle (a sink holds a closure over the body).

Nothing here draws a random number.
"""

import copy
import dataclasses


@dataclasses.dataclass
class BuildRecord:
    """
    The build-time state of one body. Every field is `None` (or `0.0` / `False` where said) until a pass sets it.

    Attributes:
        air_replay (tuple or None): Where the air's draws began (the stream state, the distance override, the host eccentricity and what the late
            accretion settled), kept for `planets.remake_air` once the system's belts are in (GEN.1493).
        air_eccentricity (float or None): The orbit's eccentricity the air was made on (`planets.orbit_moved`, GEN.1402).
        air_light_au (float or None): The distance the air, the temperature, the dose and the spectrum were made on; a later settled orbit moves
            `star_distance_au` and leaves this one.
        surface_w (float or None): The heat flow at the surface, kept for `planets.refresh_at_final_age` (GEN.839).
        retention_cap_pa (float or None): The retention cap of the first air solve (GEN.1340); the life stage keeps its oxygen inside it.
        co2_sink (callable or None): The carbon sink of the air kept (GEN.1309 A); the life stage balances the CO2 again at the temperature the oxygen
            rise leaves. A closure over the body, so a pickle leaves it out.
        birth_spin (tuple or None): `(rate, obliquity, capture draws)` the body was born with, before the tides and the Cassini state changed them.
        age_drawn (float or None): The star's age the air, class and heat were made at; `0.0` marks them stale (`planets.stale_air`).
        orbit_birth (float or None): The eccentricity the body was born on, before the tides and the secular chain settled it.
        venus_spin (object or None): The `VenusSpin` of a rocky planet, for `planets.generate_spin` to put on the rotation.
        resonance_e (float or None): The eccentricity a resonance holds the orbit on.
        frame_axis (object or None): The axis of the frame a moon's state is drawn in.
        born_moons (list or None): The moons the planet was born with (a list of bodies, shared like `moons` when a body is copied).
        bare_core (dict or None): The planet as it was drawn, before an envelope (`planets.apply_envelope`, GEN.1323), so applied again it never compounds.
        captured_envelope (float): The envelope mass fraction the planet captured.
        late_core (bool): A core that finished growing after the gas disk was gone and captured no envelope (GEN.1323, GEN.1510).
        star_distance_au (float or None): The distance from the star the body's exposure was worked out on (a moon's planet's).
        deferred_climate (tuple or None): `(moon_count, distance_override)` of a body built with `defer_climate` whose air and moons `finish_climate` still owes (GEN.1282).
        kozai (object or None): The `Kozai` outcome of a moon that cycles.
        late_accreted_percent (float or None): The percent of the body's mass the late accretion added (`planets.remake_air`).
        plan_eccentricity (float or None): The eccentricity the planet architecture drew for the body's slot and held to its neighbours (`planet_architecture.Seed.eccentricity`, GEN.1797), the birth value
            `planets.generate_orbit_shape` takes in place of a draw (cut at what the star's Roche limit allows); `None` for a body whose orbit is drawn there (a slot the author named, a moon, a body made on its own).
        moon_history_unfinished (bool): No integrator finished the planet's moon history (`moon_orbits.IntegrationFailed`, GEN.1816): the moons that stay are as born, read by `moon_orbits.inspiral_times_s` alone.
        air_input (tuple or None): `(setting, start_k, entry_state)` of the final air solve (`air_final.make`): what `air_history.temperatures` solves the air again from at each age (GEN.1745).
    """

    air_replay: tuple = None
    air_eccentricity: float = None
    air_light_au: float = None
    surface_w: float = None
    retention_cap_pa: float = None
    co2_sink: object = None
    birth_spin: tuple = None
    age_drawn: float = None
    orbit_birth: float = None
    venus_spin: object = None
    resonance_e: float = None
    frame_axis: object = None
    born_moons: list = None
    bare_core: dict = None
    captured_envelope: float = 0.0
    late_core: bool = False
    star_distance_au: float = None
    deferred_climate: tuple = None
    kozai: object = None
    late_accreted_percent: float = None
    plan_eccentricity: float = None
    air_input: tuple = None
    moon_history_unfinished: bool = False

    def __getstate__(self):
        state = dict(self.__dict__)
        state["co2_sink"] = None
        return state

    def __copy__(self):
        clone = object.__new__(type(self))
        clone.__dict__.update(self.__dict__)        # the sink too: a pickle leaves it out, a copy keeps it
        return clone

    def __deepcopy__(self, memo):
        clone = object.__new__(type(self))
        memo[id(self)] = clone
        for name, value in self.__dict__.items():
            clone.__dict__[name] = value if name == "co2_sink" else copy.deepcopy(value, memo)
        return clone


_NOTHING = BuildRecord()
"""BuildRecord: What an object that builds nothing reads (never written: a pass that sets a field reads `body.build` of a real body)."""


def record_of(body):
    """`body.build`, or an empty record for an object that has none (a row restored for display, a stand-in, a belt): a reader of any object asks here."""
    record = getattr(body, "build", None)
    return record if record is not None else _NOTHING


class HasBuildRecord:
    """A class whose objects carry a `BuildRecord` as `build`, made the first time it is read (a body rebuilt from the database never ran `__init__`)."""

    @property
    def build(self):
        record = self.__dict__.get("_build")
        if record is None:
            record = self.__dict__["_build"] = BuildRecord()
        return record

    @build.setter
    def build(self, record):
        self.__dict__["_build"] = record

    def _copy_build_record(self, clone):
        """Gives `clone` (a shallow copy of this body) a record of its own, so a pass on the copy leaves this body's as it was."""
        if "_build" in self.__dict__:
            clone.__dict__["_build"] = copy.copy(self.__dict__["_build"])
