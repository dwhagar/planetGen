# tests/test_web_facility_types.py

"""
The admin Facility types page and the facility form's type, fields,
affiliation and ratings (POP.10), against a real throwaway database through
the in-process transport.
"""

import re

from planetgen.db import query, store
from planetgen.web import csrf
from tests.publicids import pid
from tests.test_web_facilities import (  # noqa: F401 (fixtures and helpers)
    AU_KM, _csrf, _login_fresh_admin, _panel, _system_with_terrestrial, db_app,
)

PAGE = "/admin/facility-types"


def _post(app, client, **form):
    return client.post(PAGE, data={csrf.FIELD_NAME: _csrf(app, client), **form})


def _followed(client, response):
    assert response.status_code == 303
    return client.get(response.headers["Location"]).get_data(as_text=True)


def test_the_page_needs_an_admin(db_app, mysql_config):
    store.get_connection(mysql_config).close()
    response = db_app.test_client().get(PAGE)
    assert response.status_code == 302 and "/login" in response.headers["Location"]


def test_an_admin_manages_types_and_fields_from_the_page(db_app, mysql_config):
    store.get_connection(mysql_config).close()
    client = db_app.test_client()
    _login_fresh_admin(client, mysql_config)

    empty = client.get(PAGE).get_data(as_text=True)
    assert "Add the example type" in empty and 'aria-current="page">Facility types' in empty

    page = _followed(client, _post(db_app, client, type_action="examples"))
    assert "Added Research outpost." in page and "Research outpost" in page and "<code>staff</code>" in page
    assert "Add the example type" not in page
    assert "The example type is there already." in _followed(client, _post(db_app, client, type_action="examples"))

    page = _followed(client, _post(db_app, client, type_action="create_type", name="Farm world", kind="colony",
                                   icon="tree", description="Grows food.", sort_order="2", active="1"))
    assert "Added the type Farm world." in page and 'name="sort_order"' in page
    type_id = int(re.search(r'id="type-(\d+)"[^>]*>\s*<h2[^>]*><sl-icon name="tree"', page).group(1))

    clash = _followed(client, _post(db_app, client, type_action="create_type", name="Farm world", kind="colony"))
    assert "already exists" in clash and 'class="error"' in clash

    page = _followed(client, _post(db_app, client, type_action="add_field", type_id=str(type_id), key="hectares",
                                   label="Farmland", datatype="decimal", unit="ha", min_value="0", required="1"))
    assert "Added the field hectares." in page and "Farmland" in page and "decimal (ha)" in page
    page = _followed(client, _post(db_app, client, type_action="add_field", type_id=str(type_id), key="crop",
                                   label="Crop", datatype="choice", choices="Grain\nFruit\n\nFungus\n",
                                   default_text="Grain"))
    assert "Grain, Fruit, Fungus" in page
    bad = _followed(client, _post(db_app, client, type_action="add_field", type_id=str(type_id), key="Bad Key",
                                  label="Bad", datatype="text"))
    assert 'class="error"' in bad and "key" in bad
    assert "must be a number" in _followed(client, _post(db_app, client, type_action="add_field", type_id=str(type_id),
                                                         key="ok", label="Ok", datatype="int", min_value="low"))

    conn = store.get_connection(mysql_config)
    try:
        fields = {r["key"]: r["id"] for r in conn.execute(
            "SELECT id, `key` FROM facility_type_fields WHERE type_id = ?", (type_id,)).fetchall()}
    finally:
        conn.close()
    page = _followed(client, _post(db_app, client, type_action="update_field", type_id=str(type_id),
                                   field_id=str(fields["hectares"]), label="Land", datatype="decimal", unit="",
                                   min_value="1", max_value="900", required="1", sort_order="0"))
    assert "Saved the field." in page and "Land" in page and "1.0 to 900.0" in page
    page = _followed(client, _post(db_app, client, type_action="delete_field", type_id=str(type_id),
                                   field_id=str(fields["crop"])))
    assert "Deleted the field." in page and "Grain, Fruit" not in page

    page = _followed(client, _post(db_app, client, type_action="update_type", type_id=str(type_id), name="Farm planet",
                                   kind="colony", icon="", sort_order="2"))      # `active` unchecked: out of use
    assert "Saved Farm planet." in page and "out of use" in page
    page = _followed(client, _post(db_app, client, type_action="delete_type", type_id=str(type_id)))
    assert "Deleted the type." in page and "Farm planet" not in page

    assert _post(db_app, client, type_action="nonsense").status_code == 400
    assert db_app.test_client().post(PAGE, data={"type_action": "examples"}).status_code in (302, 400, 403)


def test_a_typed_facility_with_ratings_is_placed_edited_and_suggested_from_the_system_page(db_app, mysql_config):
    system_id, ids = _system_with_terrestrial(mysql_config)
    client = db_app.test_client()
    _login_fresh_admin(client, mysql_config)
    _followed(client, _post(db_app, client, type_action="examples"))
    conn = store.get_connection(mysql_config)
    try:
        type_id = conn.execute("SELECT id FROM facility_types").fetchone()["id"]
        star = conn.execute("SELECT id FROM stars WHERE star_system_id = ?", (system_id,)).fetchone()["id"]
    finally:
        conn.close()
    url = f"/system/{pid('system', system_id)}"

    page = client.get(url).get_data(as_text=True)
    assert f'<option value="{type_id}" data-kind="outpost">Research outpost</option>' in page

    response = client.post(url, data={
        csrf.FIELD_NAME: _csrf(db_app, client), "facility_action": "save", "name": "Quiet Ear",
        "host": f"star:{pid('star', star)}", "placement": "orbital", "kind": "outpost", "orbit_step": "50",
        "type_id": str(type_id), "field_staff": "30", "field_focus": "Relay", "affiliation_name": "Free Traders",
        "rating_crime": "yellow", "rating_health": "red",
    })
    assert response.status_code == 303, response.get_data(as_text=True)[:2000]
    conn = store.get_connection(mysql_config)
    try:
        (saved,) = query.facilities_for_system(conn, system_id)
    finally:
        conn.close()
    assert (saved["type_name"], saved["affiliation_name"], saved["overall"]) == ("Research outpost", "Free Traders", "red")
    assert saved["fields"] == {"staff": 30, "focus": "Relay"}

    panel = _panel(client.get(url).get_data(as_text=True))
    assert "Research outpost" in panel and "Free Traders" in panel and "focus: Relay, staff: 30" in panel
    assert 'aria-label="Health: Red, crisis"' in panel and 'aria-label="Crime: Yellow, elevated"' in panel

    facility = pid("facility", saved["id"])
    response = client.post(url, data={csrf.FIELD_NAME: _csrf(db_app, client), "facility_action": "suggest",
                                      "facility_id": facility})
    assert "Ratings suggested" in _followed(client, response)
    conn = store.get_connection(mysql_config)
    try:
        (suggested,) = query.facilities_for_system(conn, system_id)
    finally:
        conn.close()
    assert suggested["ratings"]["crime"] == "yellow" and suggested["ratings"]["health"] == "red"   # kept
    assert all(color is not None for color in suggested["ratings"].values())                      # the rest filled

    response = client.post(url, data={csrf.FIELD_NAME: _csrf(db_app, client), "facility_action": "edit",
                                      "facility_id": facility, "affiliation": "", "affiliation_name": "",
                                      "rating_crime": "green", "ratings_note": "Reviewed"})
    assert "Facility changed." in _followed(client, response)
    conn = store.get_connection(mysql_config)
    try:
        (edited,) = query.facilities_for_system(conn, system_id)
    finally:
        conn.close()
    assert edited["ratings"]["crime"] == "green" and edited["ratings"]["health"] is None      # replaced as a set
    assert edited["affiliation_name"] is None and edited["ratings_note"] == "Reviewed"
