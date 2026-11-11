# planetgen/describe/influx_row.py

"""
The derived row of a planet's sheet for the belt dust and comets that reach it (GEN.1440, GEN.1453 slice A)
===========================================================================================================

`influx_kg_per_yr(conn, planet_id)` reads the planet, the belts, the comets and the giants stored for its star and gives `physics/influx.planet_influx_kg_per_yr`. It is worked out when the sheet is
made and is not a column: it is a pure function of stored data. `None` for a giant, a planet with no star or orbit, and a system with nothing to send in.
"""

from planetgen.physics import constants, influx
from planetgen.physics import unit_registry as units

INFLUX_KEY = "influx_kg_per_yr"
"""str: The key of the planet's derived row for the mass of belt dust and comets that reaches it each year."""

_kg_to_msun = units.converter("kg", "Msun")
_kg_to_mearth = units.converter("kg", "Mearth")
_km_to_rearth = units.converter("km", "Rearth")


def _au(km):
    return float(km) / constants.AU_TO_KM


def influx_kg_per_yr(conn, planet_id):
    """The belt dust and comet mass, kg a year, that reaches the stored planet `planet_id` (`physics/influx`), or `None` where there is nothing to say."""
    planet = conn.execute("SELECT star_system_id, star_id, body_type, distance_km, radius_km FROM planets WHERE id = ?", (planet_id,)).fetchone()
    if planet is None or planet["body_type"] == "g" or not planet["distance_km"] or not planet["radius_km"]:
        return None
    star = (conn.execute("SELECT mass_kg, age_gy, luminosity_w FROM stars WHERE id = ?", (planet["star_id"],)).fetchone() if planet["star_id"]
            else conn.execute("SELECT mass_kg, age_gy, luminosity_w FROM stars WHERE star_system_id = ? ORDER BY id LIMIT 1", (planet["star_system_id"],)).fetchone())
    if star is None or not star["mass_kg"] or not star["age_gy"]:
        return None
    belts = conn.execute("SELECT star_id, mass_kg, lower_limit_km, upper_limit_km FROM asteroid_belts WHERE star_system_id = ?", (planet["star_system_id"],)).fetchall()
    comets = conn.execute("SELECT star_id, mass_kg, perihelion_distance_km FROM comets WHERE star_system_id = ?", (planet["star_system_id"],)).fetchall()
    giants = conn.execute("SELECT star_id, mass_kg, radius_km, distance_km FROM planets WHERE star_system_id = ? AND body_type = 'g' AND id <> ?",
                          (planet["star_system_id"], planet_id)).fetchall()
    same = lambda row: row["star_id"] == planet["star_id"]
    reservoirs = [influx.Reservoir("belt", float(belt["mass_kg"]), _au(belt["lower_limit_km"]), _au(belt["upper_limit_km"]), float(star["age_gy"]) * 1.0e9)
                  for belt in belts if same(belt) and belt["mass_kg"] and belt["lower_limit_km"] and belt["upper_limit_km"]]
    reservoirs += [influx.Reservoir("comet", float(comet["mass_kg"]), _au(comet["perihelion_distance_km"]), _au(comet["perihelion_distance_km"]))
                   for comet in comets if same(comet) and comet["mass_kg"] and comet["perihelion_distance_km"]]
    if not reservoirs:
        return None
    big = [influx.Giant(_kg_to_mearth(giant["mass_kg"]), _km_to_rearth(giant["radius_km"]), _au(giant["distance_km"]))
           for giant in giants if same(giant) and giant["mass_kg"] and giant["radius_km"] and giant["distance_km"]]
    luminosity = float(star["luminosity_w"]) / constants.SOLAR_LUMINOSITY if star["luminosity_w"] else None
    total, _ = influx.planet_influx_kg_per_yr(float(planet["radius_km"]), _au(planet["distance_km"]), _kg_to_msun(star["mass_kg"]), reservoirs, big, luminosity)
    return total
