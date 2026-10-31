"""
A `bright_stars` row for `store.insert_bright_stars`, built by column name (TEST.478, TEST.489), so a column added to the table (GEN.846 added
`rotation_period_hours`) cannot leave a hand-written tuple short.
"""

from planetgen.db import store

DEFAULTS = {"population": "young", "star_type": "B", "yerkes_class": "III", "mass_kg": 1e31, "radius_km": 1e7,
            "temperature_k": 15000.0, "luminosity_w": 1e28, "age_gy": 0.05, "lifespan_gy": 0.1, "initial_mass_sol": 8.0,
            "phase_end_age_gy": None, "seed": 1, "metallicity_fe_h": None, "rotation_period_hours": None}
"""dict: A scattered star's values by `store.BRIGHT_STAR_COLUMNS` name, less its address and position; `star_row` fills the rest from it."""

B2V = {"star_type": "B2V", "yerkes_class": "V", "mass_kg": 1.4e31, "radius_km": 3.0e6, "temperature_k": 22000.0, "age_gy": 0.02,
       "lifespan_gy": 0.03, "initial_mass_sol": 7.0, "phase_end_age_gy": 0.03, "seed": 42}
"""dict: The B2V star the bright-star tests place (`luminosity_w` is theirs to give)."""


def star_row(address, position_milli, **values):
    """A `bright_stars` row in `store.BRIGHT_STAR_COLUMNS` order: `address` is `(shell, band, wedge)`, `position_milli` the three position
    columns, `values` the others by name."""
    columns = store.BRIGHT_STAR_COLUMNS
    unknown = set(values) - set(columns)
    assert not unknown, f"not bright_stars columns: {sorted(unknown)}"
    row = {**dict(zip(columns[:3], address)), **dict(zip(columns[3:6], position_milli)), **DEFAULTS, **values}
    return tuple(row[name] for name in columns)
