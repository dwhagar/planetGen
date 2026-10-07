# planetgen/physics/units.py

"""
Units
=====

Distance conversions between light-years, parsecs, milliparsecs and AU.
Plain conversions: non-finite values pass through unchanged.
"""

from planetgen.physics import constants as physical_constants


def ly_to_milliparsecs(ly):
    """
    Converts a distance in light-years to milliparsecs (mpc).

    For the database persistence layer only -- sector-scale position/
    geometry columns (star_systems.position_x/y/z_mpc, sectors.edge_mpc;
    see planetgen/db/schema.sql) are stored in milliparsecs specifically,
    distinct from every other distance-shaped column in the schema (which
    is kilometers). Generation/physics code keeps its own native light-year
    units for sector geometry throughout and never calls this.

    Args:
        ly (float): The distance in light-years.

    Returns:
        float: The distance in milliparsecs.
    """
    au = ly * physical_constants.LY_TO_AU
    return au / physical_constants.AU_PER_MILLIPARSEC


def milliparsecs_to_ly(mpc):
    """
    Converts a distance in milliparsecs (mpc) back to light-years -- the
    inverse of `ly_to_milliparsecs`, for reconstructing live objects
    (native light-year sector geometry) from database rows.

    Args:
        mpc (float): The distance in milliparsecs.

    Returns:
        float: The distance in light-years.
    """
    au = mpc * physical_constants.AU_PER_MILLIPARSEC
    return au * physical_constants.AU_TO_LY


def mpc_to_pc(mpc):
    """
    Converts a distance in milliparsecs (mpc) to parsecs (pc). Exact --
    milliparsecs and parsecs are the same unit family, `1/1000` apart, so
    this is a single power-of-ten scaling with no AU round-trip (contrast
    `milliparsecs_to_ly`, which does need one).

    Args:
        mpc (float): The distance in milliparsecs.

    Returns:
        float: The distance in parsecs.
    """
    return mpc / 1000


def pc_to_mpc(pc):
    """
    Converts a distance in parsecs (pc) to milliparsecs (mpc) -- the
    inverse of `mpc_to_pc`. Exact, same reasoning.

    Args:
        pc (float): The distance in parsecs.

    Returns:
        float: The distance in milliparsecs.
    """
    return pc * 1000


def pc_to_ly(pc):
    """
    Converts a distance in parsecs (pc) to light-years (ly), for
    human-readable display (e.g. "~48,923 ly from the galactic core") --
    see `docs/design/galaxy-coordinate-system.md` section 2. Not used by
    the database persistence layer itself (which stores galaxy-scale
    distances in parsecs directly); this exists purely for prose alongside
    the same "display string next to the raw stored value" treatment
    `table_*` columns get elsewhere in this schema.

    Args:
        pc (float): The distance in parsecs.

    Returns:
        float: The distance in light-years.
    """
    return pc * (physical_constants.AU_PER_PARSEC / physical_constants.LY_TO_AU)


def ly_to_pc(ly):
    """
    Converts a distance in light-years (ly) to parsecs (pc) -- the inverse
    of `pc_to_ly`. See that function's docstring.

    Args:
        ly (float): The distance in light-years.

    Returns:
        float: The distance in parsecs.
    """
    return ly * (physical_constants.LY_TO_AU / physical_constants.AU_PER_PARSEC)


def ly_to_au(ly):
    """
    Converts a distance in light-years (ly) to astronomical units (AU),
    for human-readable display at a stellar-phenomenon's own AU scale
    (`planetgen/web/maps/phenomenonmap.py`'s diagram) -- the AU counterpart to
    `pc_to_ly`/`ly_to_pc` above, using the same
    `physical_constants.LY_TO_AU` this package's other ly<->AU
    conversions already share.

    Args:
        ly (float): The distance in light-years.

    Returns:
        float: The distance in astronomical units.
    """
    return ly * physical_constants.LY_TO_AU


def au_to_ly(au):
    """
    Converts a distance in astronomical units (AU) to light-years (ly) --
    the inverse of `ly_to_au`. See that function's docstring.

    Args:
        au (float): The distance in astronomical units.

    Returns:
        float: The distance in light-years.
    """
    return au * physical_constants.AU_TO_LY
