# planetgen/galaxy/components.py

"""
Galaxy components (GEN.159): the discrete stellar components of a galaxy beyond the disk, one record each.

A globular cluster (and later a satellite or a stream, GEN.237) is stored as a density, not as stars: a profile of its
mass about a centre. `Component` is one `galaxy_components` row (`db/store.py` writes and reads them); the profile
itself is `physics/king.py` for the `king` profile. Positions are in parsecs in the galaxy frame (the galactic centre
at the origin), velocities in km/s there, both at `epoch_unix`.
"""

import dataclasses

from planetgen import tuning

KINDS = ("cluster", "satellite", "stream")
"""tuple: The kinds of component a row may have (the column's CHECK)."""

SOURCES = ("catalogue", "generated")
"""tuple: Where a row came from: a real catalogue's entry, or a draw from the galaxy seed."""

PROFILES = ("king", "exponential", "plummer")
"""tuple: The density profiles a row may name: a cluster's `king` (`physics/king.py`), a satellite's `exponential` or
`plummer` (`physics/component_profiles.py`)."""


@dataclasses.dataclass(frozen=True)
class Component:
    """
    One component of a galaxy.

    Attributes:
        kind (str): One of `KINDS`.
        name (str): The component's name, unique within its kind.
        source (str): One of `SOURCES`.
        center_pc (tuple): `(x, y, z)` in parsecs in the galaxy frame.
        velocity_kms (tuple): `(vx, vy, vz)` in km/s in the galaxy frame.
        epoch_unix (float): When the position and velocity hold.
        mass_solar (float): Total mass.
        half_mass_radius_pc (float): The radius holding half the mass.
        profile (str): One of `PROFILES`.
        age_gy (float): The stars' age in Gyr.
        metallicity_fe_h (float): The stars' [Fe/H] in dex against the Sun.
        seed (int): The seed its stars are drawn from.
        king_w0 (int or None): The King central potential, 3 to 12, for the `king` profile and no other.
        tidal_radius_pc (float): The radius the profile is cut at: the King model's own, or the cut of an exponential or
            Plummer profile (`component_profiles.cut_radius_pc`).
        core_collapsed (bool): Whether the core has collapsed.
        catalogue_id (str or None): The catalogue's own name for it, for matching a later import.
        parent_component_id (int or None): The row of the component it lies in.
        id (int or None): Its row, once stored.
    """

    kind: str
    name: str
    source: str
    center_pc: tuple
    velocity_kms: tuple
    epoch_unix: float
    mass_solar: float
    half_mass_radius_pc: float
    profile: str
    age_gy: float
    metallicity_fe_h: float
    seed: int
    tidal_radius_pc: float = None
    king_w0: int = None
    core_collapsed: bool = False
    catalogue_id: str = None
    parent_component_id: int = None
    id: int = None

    def __post_init__(self):
        if self.kind not in KINDS:
            raise ValueError(f"unknown component kind {self.kind!r}")
        if self.source not in SOURCES:
            raise ValueError(f"unknown component source {self.source!r}")
        if self.profile not in PROFILES:
            raise ValueError(f"unknown component profile {self.profile!r}")
        if (self.profile == "king") != (self.king_w0 is not None):
            raise ValueError("a king profile needs king_w0, and no other profile has one")
        if self.tidal_radius_pc is None or self.tidal_radius_pc <= 0.0:
            raise ValueError("a component needs a positive tidal_radius_pc: the radius its profile is cut at")

    @property
    def expected_systems(self):
        """How many systems the component holds: its mass over `tuning.CLUSTER_MEAN_SYSTEM_MASS_SOL`."""
        return self.mass_solar / tuning.CLUSTER_MEAN_SYSTEM_MASS_SOL
