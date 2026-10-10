# planetgen/api/idpaths.py

"""
Where each API answer carries an object ID (API.23).

`ID_PATHS` maps an endpoint name to the rules `ids.translate` applies to its JSON answer: the row ids
found at those places leave the API as printed IDs. A path is dotted keys with `[]` for each list level;
a trailing `{}` renames the keys of that object. Keys in `ids.KEY_KINDS` (`sector_id`, `star_system_id`,
`system_id`, `capital_system_id`, `ref`) are translated wherever they appear in an answer that has rules,
so they are not repeated here. An endpoint absent from the table answers no object IDs (an entry with no rules still gets those keys).

Only the kinds in `ids.ACTIVE_KINDS` are translated (every kind has been since API.23 stage 2).
"""

_NEAREST = {"nearest[].id": "system"}

_BODIES = {
    "stars[].id": "star",
    "planets[].id": "planet",
    "planets[].moons[].id": "moon",
    "belts[].id": "belt",
    "comets[].id": "comet",
}

ID_PATHS = {
    "api.sectors": {"items[].id": "sector"},
    "api.sector_detail": {
        "id": "sector",
        "systems[].id": "system",
        "systems[].nearest[].id": "system",
        "phenomena[].id": "phenomenon",
        "phenomena[].nearest[].id": "system",
    },
    "api.systems": {"items[].id": "system"},
    "api.system_detail": {
        "id": "system",
        "nearest_neighbors[].id": "system",
        "sector_siblings[].id": "system",
        **_BODIES,
    },
    "api.system_sections": {"stars{}": "star", "planets{}": "planet", "moons{}": "moon", "belts{}": "belt",
                            "comets{}": "comet"},
    "api.system_scene": {
        "system.id": "system",
        **_BODIES,
        "stars[].orbit.around": "ref",
        "planets[].orbit.around": "ref",
        "planets[].moons[].orbit.around": "ref",
        "planets[].moons[].parent": "ref",
        "belts[].around": "ref",
        "comets[].orbit.around": "ref",
    },
    "api.phenomena": {"items[].id": "phenomenon"},
    "api.phenomenon": {"id": "phenomenon", **_NEAREST},
    "api.object_ref_route": {"id": "by-kind", "siblings[]": "ref"},
    "api.nav": {
        "route.path[]": "navkey",
        "route.hops[].from": "navkey",
        "route.hops[].to": "navkey",
        "route.positions{}": "navkey",
        "legs[].from": "ref",
        "legs[].to": "ref",
    },
    "api.near": {"rows[].id": "by-kind"},
    "api.galaxy_sectors": {"items[].id": "sector"},
    "api.galaxy_phenomena": {"items[].id": "phenomenon"},
    "api.galaxy_made": {"items[].id": "sector"},
    "api.galaxy_locate_route": {"matches[].id": "by-kind"},
    "api.search": {
        "results.sectors.rows[].id": "sector",
        "results.systems.rows[].id": "system",
    },
    "api.nebula_shape_route": {"id": "nebula"},
    "population.polity_detail": {"systems[].id": "system"},
    "population.territories": {"points[].id": "system"},
    "api.system_text": {"id": "system"},
    "api.create_sector": {"id": "sector"},
    "api.create_system": {"id": "system"},
    "api.update_system": {"id": "system"},
    # Answers whose only IDs are the keys every endpoint with rules translates (`star_system_id`, `sector_id`).
    "api.job_status_route": {},
    "api.sector_facilities": {"items[].id": "facility", "facilities[].id": "facility"},
    "api.system_facilities": {"items[].id": "facility", "facilities[].id": "facility"},
    "api.facility": {"id": "facility"},
    "api.create_facility": {"id": "facility"},
    "api.galaxy_cell": {},
    "api.galaxy_uncharted_sector": {},
    "api.rename_star": {"id": "star"},
    "api.rename_planet": {"id": "planet"},
    "api.rename_moon": {"id": "moon"},
    "population.polity_list": {},
    "population.species_list": {},
    "population.species_detail": {},
    "population.system_owner": {},
    "edits.regenerate_planet": {},
    "edits.regenerate_moon": {},
    "edits.regenerate_belt": {},
    "edits.change_planet_class": {},
    "edits.change_moon_class": {},
    "edits.change_system_star": {},
    "edits.system_class_options": {"recommended{}": "ref"},
    "population.planet_species": {},
}
"""dict: Endpoint name -> `{path: kind}` rules for `ids.translate`."""
