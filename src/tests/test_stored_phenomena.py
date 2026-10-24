# tests/test_stored_phenomena.py

"""
UX.185: a stored phenomenon says the same paragraph the generator's description says, and its page has the previous and
next phenomenon of its type. Each test saves a generated object, reads its row back the way the API serves it, and
compares the words.
"""

import pytest

from planetgen.db import query, store
from planetgen.describe import describe_phenomenon
from planetgen.describe.html import blocks_html
from planetgen.describe.model import sections
from planetgen.describe.stored_phenomena import TYPES, stored_phenomenon_document
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.asteroid_field import AsteroidField
from planetgen.generation.phenomena.compact_remnant import BlackHole, NeutronStar
from planetgen.generation.phenomena.nebula import Nebula
from planetgen.generation.phenomena.quasar import Quasar
from planetgen.generation.phenomena.rogue import InterstellarComet, RoguePlanet
from planetgen.generation.phenomena.supernova_remnant import SupernovaRemnant

MAKERS = {
    "nebula": (Nebula, store.insert_nebula), "quasar": (Quasar, store.insert_quasar),
    "rogue_planet": (RoguePlanet, store.insert_rogue_planet),
    "interstellar_comet": (InterstellarComet, store.insert_interstellar_comet),
    "asteroid_field": (AsteroidField, store.insert_asteroid_field),
    "supernova_remnant": (SupernovaRemnant, store.insert_supernova_remnant),
    "black_hole": (BlackHole, store.insert_black_hole), "neutron_star": (NeutronStar, store.insert_neutron_star),
}


def _prose(document):
    return [blocks_html([block for block in piece["blocks"] if block["type"] == "prose"]) for piece in sections(document)][:1]


def test_every_type_has_a_maker():
    assert set(MAKERS) == set(TYPES)


@pytest.mark.parametrize("kind", TYPES)
def test_a_stored_phenomenon_reads_the_paragraph_the_generated_one_does(mysql_config, kind):
    cls, insert = MAKERS[kind]
    made = cls(SystemConfig(), name="Testa")
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            row_id = insert(conn, made)
        detail = query.phenomenon_detail(conn, kind, row_id)
    finally:
        conn.close()
    stored = stored_phenomenon_document(kind, detail)
    assert stored["title"] == made.name
    assert _prose(stored)[0], kind
    generated = _prose(describe_phenomenon(made))
    assert _prose(stored) == generated, (kind, _prose(stored), generated)


def test_the_neighbours_are_the_free_phenomena_of_the_type_before_and_after(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            ids = [store.insert_nebula(conn, Nebula(SystemConfig(), name=f"Cloud {n}")) for n in range(3)]
        first, middle, last = (query.phenomenon_detail(conn, "nebula", row_id) for row_id in ids)
    finally:
        conn.close()
    assert first["previous"] is None and first["next"]["name"] == "Cloud 1"
    assert middle["previous"]["name"] == "Cloud 0" and middle["next"]["name"] == "Cloud 2"
    assert last["next"] is None and last["previous"]["name"] == "Cloud 1"
