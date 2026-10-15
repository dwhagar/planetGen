# planetgen/web/system_species.py

"""
The system page's "Species that could live here" panel (POP.12): for each rocky planet and moon that a species
could live on, the species that could live there unaided and those that could with the equipment their tech level
can make, each with its rating in the six habitability domains. The API (`GET /api/systems/<id>/habitable-species`)
does the scoring; this builds the rows.
"""

from planetgen.physics.habitability_world import DOMAINS
from planetgen.web.lib import apiclient
from planetgen.web.lib.fmt import esc

from .helpers import page_url, population_status, trusted_html

TIER_INDEX = {"Blue": 0, "Green": 1, "Yellow": 2, "Red": 3}
"""dict: A domain's colour word and the `flag-tier-N` class that draws it."""


def _chips(domains):
    """The six domains as chips, each with its name, its colour as a word and its class, so no colour stands alone."""
    return trusted_html("".join(
        f'<span class="flag flag-tier flag-tier-{TIER_INDEX[domains[name]["tier"]]}">'
        f'{esc(name.capitalize())}: {esc(domains[name]["tier"])}</span>' for name in DOMAINS))


def _species_rows(entries):
    return [{"name": e["name"], "url": page_url("species_page", species_id=e["id"]), "chips": _chips(e["domains"]),
             "equipment": e["equipment"], "tech": e.get("tech_level")} for e in entries]


def panel_rows(db, system_id):
    """`[{name, kind, unaided, equipped}]` for the panel: only bodies some species could live on."""
    if not population_status()["species"]:
        return []
    try:
        bodies = apiclient.get_system_habitable_species(db, system_id)
    except apiclient.ApiError:
        return []  # the panel is extra: a system page that cannot read it still opens
    return [{"name": body["name"], "kind": body["kind"], "unaided": _species_rows(body["unaided"]),
             "equipped": _species_rows(body["equipped"])} for body in bodies]
