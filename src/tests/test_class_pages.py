# tests/test_class_pages.py

"""
The class reference pages (`web/class_pages.py`, `planetgen/web/lib/classref.py`):
every type and class page renders, the catalog's values come from the
generator's own tables, an unknown type or code is a 404, star types
parse into their spectral and luminosity classes, and the system and
phenomenon pages link their class labels here. No database: the system
and phenomenon pages' data layer is faked like `test_web_system_phen.py`.
"""

import re

import pytest

from planetgen.web.app import create_app
from planetgen.api.config import Config

from planetgen.describe import model  # noqa: E402
from planetgen.web.lib import apiclient  # noqa: E402
from planetgen.web.lib import classref  # noqa: E402
from markupsafe import escape  # noqa: E402
from planetgen.physics import constants, planet_classes  # noqa: E402
from planetgen import tuning
from planetgen.util import chem
from planetgen.generation.comet import PERIOD_CLASS_LABELS  # noqa: E402
from planetgen.generation.phenomena.nebula import NEBULA_CLASS_LETTERS, REMNANT_CLASS_LETTERS  # noqa: E402
from planetgen.physics.stellar_evolution import YERKES_CLASS_NAMES  # noqa: E402
from planetgen.web.lib.systempage import star_notes_html, stars_html, system_list_html  # noqa: E402

DB = "planetgen_web_test"


class _Config(Config):
    WEB_DATABASE = DB
    SESSION_COOKIE_SECURE = False
    SECRET_KEY = "test-secret"


@pytest.fixture(scope="module")
def app():
    application = create_app(_Config)
    application.testing = True
    return application


@pytest.fixture
def client(app):
    return app.test_client()


def _html(client, path):
    resp = client.get(path)
    assert resp.status_code == 200, path
    return resp.get_data(as_text=True)


# --- Every page renders ---------------------------------------------------------------

def test_index_lists_every_type_with_its_count(client):
    html = _html(client, "/classes")
    for slug, entry in classref.catalog().items():
        assert f'<a href="/classes/{slug}">{escape(entry["name"])}</a>' in html
        assert f"<td>{len(entry['classes'])}</td>" in html
    # UX.180: Classes is under the Help menu, which is current on it.
    assert '<a href="/classes" aria-current="page">Classes</a>' not in html
    assert 'class="site-help-trigger" aria-haspopup="menu" aria-current="page"><span class="site-help-text">Help</span>' in html
    assert 'data-href="/classes" aria-current="page">Classes</sl-menu-item>' in html


def test_index_lists_the_types_alphabetically(client):
    """UX.77: the class categories are in alphabetical order."""
    html = _html(client, "/classes")
    names = [escape(entry["name"]) for entry in classref.catalog().values()]
    positions = [html.index(f'">{name}</a>') for name in sorted(names, key=str.casefold)]
    assert positions == sorted(positions), "the index is not in alphabetical order"


def test_alphabetical_ignores_case():
    rows = [{"name": "star"}, {"name": "Nebula"}, {"name": "asteroid"}]
    assert [r["name"] for r in classref.alphabetical(rows)] == ["asteroid", "Nebula", "star"]


def test_every_type_and_class_page_renders(client):
    for slug, entry in classref.catalog().items():
        html = _html(client, f"/classes/{slug}")
        for code, item in entry["classes"].items():
            assert f'href="/classes/{slug}/{code}"' in html
            page = _html(client, f"/classes/{slug}/{code}")
            # Breadcrumb back to the type, and every fact shown.
            assert f'<li><a href="/classes/{slug}">{escape(entry["name"])}</a></li>' in page
            for label, text in item["facts"]:
                assert f'<th scope="row">{escape(label)}</th><td>{escape(text)}</td>' in page
            assert "<script>" not in page and 'style="' not in page


def test_unknown_type_or_code_is_404(client):
    assert client.get("/classes/galaxy").status_code == 404
    assert client.get("/classes/nebula/Z").status_code == 404
    assert client.get("/classes/nebula/R").status_code == 404  # R is a remnant, not a nebula
    assert client.get("/classes/star-luminosity/D").status_code == 404  # an alias, not its own page
    assert client.get("/classes/asteroid-field/C3").status_code == 404  # the page is the letter's
    assert classref.class_url_parts("nebula", "Z") is None
    assert classref.class_url_parts("galaxy", "A") is None
    assert classref.class_url_parts("planet", None) is None


# --- The catalog is the generator's own tables -------------------------------------------

def test_catalog_types_cover_the_tables():
    cat = classref.catalog()
    assert list(cat) == ["star-spectral", "star-luminosity", "planet", "nebula", "supernova-remnant",
                         "asteroid-field", "black-hole", "rogue-planet", "comet"]
    assert list(cat["star-spectral"]["classes"]) == list(constants.SPECTRAL_CLASS_COLORS)
    assert set(cat["star-luminosity"]["classes"]) == set(YERKES_CLASS_NAMES) - {"D"}
    assert list(cat["planet"]["classes"]) == sorted(planet_classes.LETTERS)
    assert tuple(cat["nebula"]["classes"]) == NEBULA_CLASS_LETTERS
    assert tuple(cat["supernova-remnant"]["classes"]) == REMNANT_CLASS_LETTERS
    letters = {letter for family in tuning.ASTEROID_FIELD_COMPOSITIONS.values()
               for letter in family["letters"].values()}
    assert set(cat["asteroid-field"]["classes"]) == letters
    assert tuple(cat["black-hole"]["classes"]) == tuning.BLACK_HOLE_MASS_CLASSES
    assert tuple(cat["rogue-planet"]["classes"]) == tuning.ROGUE_PLANET_MASS_BIN_CHOICES
    assert list(cat["comet"]["classes"]) == [*tuning.COMET_PERIOD_CLASSES, classref.PARABOLIC_COMET_CLASS]
    # Built once and cached.
    assert classref.catalog() is cat


def test_class_pages_show_the_constants(client):
    nebula = tuning.NEBULA_CLASSES["D"]
    html = _html(client, "/classes/nebula/D")
    assert str(escape(nebula["name"])) in html
    assert chem.mark_html(nebula["species"]) in html     # UX.149: its formulas written out

    remnant = tuning.NEBULA_CLASSES["T"]
    assert str(escape(remnant["name"])) in _html(client, "/classes/supernova-remnant/T")

    html = _html(client, "/classes/planet/M")  # GEN.645: a class page gives the rule
    assert planet_classes.describe("M").sentence in html
    assert "<th scope=\"row\">Core class</th><td>Yes</td>" in html
    html = _html(client, "/classes/planet/S")
    assert "surface pressure below 10 Pa; mass at least 1.5 Earth masses" in html

    low_k, high_k = constants.TEMP_RANGES["G"]
    html = _html(client, "/classes/star-spectral/G")
    assert f"{classref.number(low_k)} to {classref.number(high_k)} K" in html
    assert constants.SPECTRAL_CLASS_COLORS["G"] in html

    html = _html(client, "/classes/star-luminosity/IA+")
    assert YERKES_CLASS_NAMES["IA+"] in html

    html = _html(client, "/classes/comet/halley_type")
    assert PERIOD_CLASS_LABELS["halley_type"] in html

    family = tuning.ASTEROID_FIELD_COMPOSITIONS["metallic"]
    html = _html(client, f"/classes/asteroid-field/{family['letters']['dense']}")
    assert family["description"][1:] in html
    # One letter can stand for several densities.
    basaltic = tuning.ASTEROID_FIELD_COMPOSITIONS["basaltic"]["letters"]
    assert basaltic["sparse"] == basaltic["typical"]
    facts = dict(classref.class_entry("asteroid-field", basaltic["sparse"])["facts"])
    assert facts["Density"] == "Sparse, typical"
    html = _html(client, "/classes/asteroid-field")
    assert "10^d to 10^(d+1) AU" in html

    low, high = tuning.BLACK_HOLE_INTERMEDIATE_MASS_RANGE_SOLAR
    assert classref.number_range(low, high, "solar masses") in _html(client, "/classes/black-hole/intermediate")

    low, high = tuning.ROGUE_BROWN_DWARF_MASS_RANGE_JUPITER
    assert f"{low:.0f} to {high:.0f} Jupiter masses" in _html(client, "/classes/rogue-planet/brown-dwarf")


def test_number_formatting():
    assert classref.number(150) == "150"
    assert classref.number(4131.4) == "4,131.4"
    assert classref.number(0.08) == "0.08"
    assert classref.number(2_000_000) == "2 × 10⁶"
    assert classref.number(1e-5) == "1 × 10⁻⁵"
    assert classref.number(0) == "0"


# --- Star types --------------------------------------------------------------------------

@pytest.mark.parametrize("star_type, expected", [
    ("G2V", ("G", "V")),
    ("G2V Yellow Main Sequence Star", ("G", "V")),
    ("B0IA", ("B", "IA")),
    ("B0IA+ Blue-White Luminous Supergiant Star", ("B", "IA+")),
    ("K5IAB", ("K", "IAB")),
    ("O5VII", ("O", "VII")),
    ("M3D", ("M", "VII")),  # D is a white dwarf, like VII
    ("A1IV", ("A", "IV")),
    ("M9VI", ("M", "VI")),
    ("F0 0", None),
    ("B00", ("B", "0")),
    ("g2v", ("G", "V")),
    ("G10V", None),
    ("X2V", None),
    ("", None),
    (None, None),
])
def test_star_type_classes(star_type, expected):
    assert classref.star_type_classes(star_type) == expected


def test_class_url_parts_resolve_aliases():
    assert classref.class_url_parts("star-luminosity", "D") == ("star-luminosity", "VII")
    assert classref.class_url_parts("asteroid-field", "C3") == ("asteroid-field", "C")
    assert classref.class_url_parts("nebula", "D") == ("nebula", "D")


# --- Links from the system and phenomenon pages ----------------------------------------------

_SECTIONS = {"overview": [], "stars": {}, "planets": {}, "moons": {}, "belts": {}, "comets": {}}


def _planet(planet_class="M"):
    return {"id": 2, "name": "Terra", "planet_class": planet_class, "body_type": "t", "habitable": True,
            "inhabited": False, "zone": "Habitable Zone", "distance_km": 1.5e8, "period_years": 1.0,
            "gravity_g": 1.0, "orbital_index": 0, "moons": [], "star_id": 1, "radius_km": 6371.0}


def _system():
    star = {"id": 1, "name": "Sol", "role": "single", "star_type": "G2V Yellow Main Sequence Star",
            "mass_kg": 1.989e30, "radius_km": 696_000.0, "temperature_k": 5778.0, "luminosity_w": 3.828e26}
    comet = {"id": 3, "name": "Kohoutek", "orbit_type": "elliptical", "period_class": "halley_type",
             "is_active": True, "perihelion_distance_km": 1e8, "eccentricity": 0.9, "orbital_period_years": 70.0,
             "star_id": 1}
    return {"stars": [star], "planets": [_planet()], "belts": [], "comets": [comet], "binary_configuration": None}


def test_system_list_links_classes(app):
    from planetgen.web.class_pages import class_url
    with app.test_request_context("/system/1"):
        html = system_list_html(_system(), _SECTIONS, class_url)
        table = stars_html(_system()["stars"], class_url)
        notes = star_notes_html(_system(), _SECTIONS, class_url)
    # UX.53: a single star is in the Stars table, not the list; its notes follow the table.
    assert "Spectral class G" not in html + notes
    # UX.182: the class letter or name is the link, with no separate line of class links.
    assert "Planet class M" not in html and "Spectral class" not in notes + table and "Luminosity class" not in notes
    assert 'Class <a href="/classes/planet/M">M</a>' in html
    assert '<a href="/classes/comet/halley_type">Halley-type comet</a>' in html
    assert ('<a href="/classes/star-spectral/G">G2</a><a href="/classes/star-luminosity/V">V</a>'
            " Yellow Main Sequence Star") in table
    # The links are in the row's <summary>, where the class is read.
    summaries = [part.split("</summary>")[0] for part in html.split("<summary>")[1:]]
    assert any('<a href="/classes/planet/M">M</a>' in summary for summary in summaries)
    assert 'class="class-links"' not in html


@pytest.mark.parametrize("orbit_type, period_class", [
    *[("elliptical", code) for code in tuning.COMET_PERIOD_CLASSES], ("parabolic", None),
])
def test_every_comet_links_its_class(app, client, orbit_type, period_class):
    """UX.29: a parabolic comet (no period class) links one too, and
    every link opens a class page."""
    from planetgen.web.class_pages import class_url
    system = _system()
    system["comets"][0].update(orbit_type=orbit_type, period_class=period_class)
    with app.test_request_context("/system/1"):
        html = system_list_html(system, _SECTIONS, class_url)
    code = period_class or classref.PARABOLIC_COMET_CLASS
    name = classref.class_entry("comet", code)["name"]
    assert f'<a href="/classes/comet/{code}">{escape(name)}</a>' in html
    assert name in _html(client, f"/classes/comet/{code}")


def test_system_list_without_hook_or_known_class_stays_text(app):
    html = system_list_html(_system(), _SECTIONS)
    assert "/classes/" not in html and "Class M" in html
    from planetgen.web.class_pages import class_url
    system = _system()
    system["planets"] = [_planet("Z")]
    system["stars"][0]["star_type"] = "Unknown"
    with app.test_request_context("/system/1"):
        html = system_list_html(system, _SECTIONS, class_url)
        table = stars_html(system["stars"], class_url)
    assert "/classes/planet" not in html and "/classes/star-" not in html + table


class _Fake:
    def __init__(self):
        self.phenomenon = None
        self.system = {
            "id": 5, "name": "Sol", "sector_id": None, "quadrant": None, "location": None, "is_binary": 0,
            "binary_type": None, "binary_configuration": None, "wikijs_url": None, "mediawiki_url": None,
            **{key: value for key, value in _system().items() if key != "binary_configuration"},
        }

    def get_system(self, db, system_id):
        return self.system

    def get_system_description(self, db, system_id):
        return model.document("system", "system:1", "Sol", [model.section("system", "system:1", "Sol")])

    def get_system_facilities(self, db, system_id):
        return []

    def get_phenomenon(self, db, phenomenon_type, phenomenon_id):
        return self.phenomenon

    def get_phenomenon_contents(self, db, phenomenon_type, phenomenon_id, limit=None):
        return {"total": 0, "counts": {}, "items": []}

    def auth_me(self, cookie_header):
        return None


@pytest.fixture
def fake(monkeypatch):
    data = _Fake()
    for name in ("get_system", "get_system_description", "get_system_facilities", "get_phenomenon", "get_phenomenon_contents", "auth_me"):
        monkeypatch.setattr(apiclient, name, getattr(data, name))
    return data


def test_system_page_links_classes(client, fake):
    html = _html(client, "/system/5")
    assert '<a href="/classes/star-spectral/G">G2</a><a href="/classes/star-luminosity/V">V</a>' in html
    assert 'Class <a href="/classes/planet/M">M</a>' in html


@pytest.mark.parametrize("kind, detail, label, url", [
    ("nebula", {"nebula_class": "D"}, "D: " + tuning.NEBULA_CLASSES["D"]["name"], "/classes/nebula/D"),
    ("supernova_remnant", {"remnant_class": "T"}, "T: " + tuning.NEBULA_CLASSES["T"]["name"],
     "/classes/supernova-remnant/T"),
    ("asteroid_field", {"field_class": "C3"}, "C3", "/classes/asteroid-field/C"),
    ("black_hole", {"mass_class": "supermassive"}, "Supermassive", "/classes/black-hole/supermassive"),
    ("rogue_planet", {"mass_bin": "brown-dwarf"}, "Brown dwarf", "/classes/rogue-planet/brown-dwarf"),
])
def test_phenomenon_page_links_its_class(client, fake, kind, detail, label, url):
    fake.phenomenon = {"id": 4, "name": "Thing", "sector_id": None, **detail}
    html = _html(client, f"/phenomenon/{kind}/4")
    assert f'<td><a href="{url}">{escape(label)}</a></td>' in html


def test_phenomenon_unknown_class_stays_text(client, fake):
    fake.phenomenon = {"id": 4, "name": "Thing", "sector_id": None, "nebula_class": "Z"}
    html = _html(client, "/phenomenon/nebula/4")
    assert "<td>Z</td>" in html and "/classes/" not in html.split("</header>", 1)[1]


def test_stars_table_hides_role_for_a_single_star(app):
    """UX.53: the Role column is for binaries and multiples."""
    from planetgen.web.class_pages import class_url
    single = _system()["stars"]
    with app.test_request_context("/system/1"):
        assert "<th>Role</th>" not in stars_html(single, class_url)
        pair = [dict(single[0], role="primary"), dict(single[0], role="secondary", name="B")]
        assert "<th>Role</th>" in stars_html(pair, class_url)


# --- UX.90: the habitability levels page --------------------------------------------

def test_habitability_levels_page_explains_the_five_levels_and_the_six_domains(client):
    from planetgen.web.lib import habitability_ref

    html = _html(client, "/help/habitability")
    chem = client.application.jinja_env.filters["chem"]
    for item, _when in habitability_ref.EQUIPMENT_ITEMS:
        assert str(chem(item)) in html
    assert "Shirtsleeve" not in html
    for domain in ("Pressure", "Temperature", "Air", "Water", "Radiation", "Suitability"):
        assert f'<h3 id="domain-{domain.lower()}">{domain}</h3>' in html
    for colour in ("Blue", "Green", "Yellow", "Red"):
        assert f"<dt>{colour}</dt>" in html
    for chip in ("Habitable class", "Habitable-class moon", "Inhabited"):
        assert f'<dt id="chip-{chip.lower().replace(" ", "-")}">{chip}</dt>' in html
    assert "pulsar" in html and "1,000 Sv a year" in html


def test_habitability_levels_page_draws_each_symbol_in_every_rating_colour(client):
    from planetgen.web.lib import habitability_symbols

    html = _html(client, "/help/habitability")
    for name, label in habitability_symbols.LABELS.items():
        for tier, colour in enumerate(("Blue", "Green", "Yellow", "Red")):
            assert f'aria-label="{label}: {colour}"' in html, (name, tier)
            assert f'#hab-{name}"' in html


def test_habitability_levels_page_draws_the_four_life_scores(client):
    html = _html(client, "/help/habitability")
    for label in ("Microbial life", "Flora", "Fauna", "Sentient life"):
        assert f'aria-label="{label}: none"' in html and f'aria-label="{label}: 9 of 9"' in html
    for name in ("microbial", "flora", "fauna", "sentient"):
        assert f'#life-{name}"' in html


def test_habitability_levels_page_says_where_the_symbols_appear(client):
    from planetgen.web.lib import habitability_ref

    html = _html(client, "/help/habitability")
    assert 'id="where-heading"' in html
    for place, _what in habitability_ref.WHERE_ROWS:
        assert place in html
    assert "/lists/planet" in html and "/lists/moon" in html


def test_habitability_levels_page_states_the_habitable_zone_from_the_generators_edges(client):
    from planetgen.physics import orbits
    from planetgen.web.lib import habitability_ref

    inner, outer = orbits.calculate_habitable_zone(constants.SOLAR_LUMINOSITY, orbits.KOPPARAPU_REFERENCE_TEFF_K)
    text = habitability_ref.zone_text()
    assert f"{inner:.2f} and {outer:.2f} AU" in text
    assert (inner, outer) == pytest.approx((0.95, 1.68), abs=0.01)
    low, high = orbits.KOPPARAPU_TEFF_RANGE_K
    assert f"{low:,.0f} to {high:,.0f} K" in text
    html = _html(client, "/help/habitability")
    assert 'id="zone-heading"' in html and "runaway greenhouse" in html


def test_the_habitability_levels_numbers_are_the_constants_the_score_uses(client):
    """GEN.219: no number on the page is typed in; it is the limit the score reads."""
    from planetgen.physics import habitability
    from planetgen.physics.protection_limits import limit
    from planetgen.web.lib import habitability_ref

    rows = {variable: cells for _domain, table in habitability_ref.domains() for variable, *cells in table}
    assert rows["Air pressure"] == ["50 to 250 kPa", "14.3 to 400 kPa", "6.3 to 10,000 kPa"]
    assert rows["Air temperature"] == ["0 to 31 C", "-20 to 45 C", "-50 to 122 C"]
    assert rows["Carbon dioxide"][0] == f"under {habitability.GAS_LIMITS_KPA['CO2'][0]:g} kPa"
    assert rows["Chaotropicity (how strongly dissolved salts pull the water's structure apart)"][0] == (
        f"under {habitability.CHAOTROPICITY_CRIT_KJ_KG:g} kJ/kg")
    assert rows["Water activity (how freely water is available)"][2] == f"{habitability.WATER_ACTIVITY_CRIT:g} or more"
    assert f"{limit('heat_shielded_c'):g} C" in habitability_ref.EQUIPMENT_RULE
    # A changed limit moves the page.
    wide = habitability.bands({"pressure_no_gear_kpa": 60.0, "dose_blue_sv_yr": 0.03})
    moved = {variable: cells for _domain, table in habitability_ref.domains(wide) for variable, *cells in table}
    assert moved["Air pressure"][0] == "60 to 250 kPa" and moved["Dose per year"][0] == "under 0.03 Sv"
    html = _html(client, "/help/habitability")
    for text in ("Chaotropicity", "too much", "too little", "Hydrogen sulfide"):
        assert text in html


def test_the_help_menu_holds_classes_and_habitability_and_the_old_links_are_gone(client):
    """UX.180: Classes leaves the section row for the Help menu beside the search button, with Habitability."""
    html = _html(client, "/classes")
    bar = html.split('<header class="site-header">')[1].split("</header>")[0]
    assert 'href="/classes"' not in bar.split("site-help-menu")[0]
    help_menu = re.search(r'<sl-dropdown class="site-help-menu".*?</sl-dropdown>', bar, re.S).group(0)
    assert re.findall(r'data-href="([^"]*)"[^>]*>([^<]*)<', help_menu) == [
        ("/classes", "Classes"), ("/help/habitability", "Habitability")]
    assert bar.index("site-help-menu") < bar.index("site-search-link") < bar.index("site-gear")
    assert "Habitability levels" not in html
    page = _html(client, "/help/habitability")
    assert 'class="site-help-trigger" aria-haspopup="menu" aria-current="page"><span class="site-help-text">Help</span>' in page
    assert 'data-href="/help/habitability" aria-current="page">Habitability</sl-menu-item>' in page
    assert "<title>Habitability" in page


def _scored_system():
    system = _system()
    system["planets"][0].update({
        "habitable": True, "equipment": {"required": [], "recommended": [], "line": "", "code": "", "items": [
            {"slot": "supply", "item": "gas-mix mask", "severity": "required", "domain": "Air"}]},
        "phi4": 1.0, "tier_pressure": 0, "tier_temperature": 0, "tier_air": 0, "tier_water": 0,
        "tier_radiation": 0, "tier_suitability": 0, "life_microbial": 3, "life_flora": 0, "life_fauna": 0,
        "life_sentient": 0})
    return system


def test_the_system_list_symbols_and_chips_link_to_their_sections_of_the_habitability_page(app, client):
    """UX.181: no separate "Habitability levels" line; each symbol, life score and chip is the link."""
    url = "/help/habitability"
    html = system_list_html(_scored_system(), _SECTIONS, habitability_url=url, sprite_url="/static/icons.svg")
    assert "Habitability levels" not in html
    for anchor in ("overall", "domain-pressure", "domain-temperature", "domain-air", "domain-water",
                   "domain-radiation", "domain-suitability", "life-microbial", "life-flora", "life-fauna",
                   "life-sentient", "chip-habitable-class", "equipment-gas-mix-mask"):
        assert f'<a class="help-link" href="{url}#{anchor}">' in html, anchor
    page = _html(client, url)
    for anchor in re.findall(r'href="/help/habitability#([^"]+)"', html):
        assert f'id="{anchor}"' in page, anchor
    # A list built with no page to link to (and the page's own legend) shows the shapes alone.
    assert "help-link" not in system_list_html(_scored_system(), _SECTIONS, sprite_url="/static/icons.svg")
    assert "help-link" not in page


def test_every_equipment_item_name_has_a_row_on_the_habitability_page(client):
    """UX.181: a chip's item name reaches a row of the equipment table, never the section's fallback."""
    from planetgen.physics import equipment
    from planetgen.web.lib import habitability_ref, habitability_symbols

    names = {*equipment.BODY_NAMES.values(), *equipment.SUPPLY_NAMES.values(), *equipment.GAS_FILTER_NAMES.values(),
             *equipment.SHELTER_NAMES.values(), *equipment.RAD_SUIT_NAMES.values(), "powered frame", "hardened suit",
             "dust filter", "dust-tight seals", "CO2 scrubber", "hardened sealed suit", "hardened heat-shielded suit"}
    page = _html(client, "/help/habitability")
    for name in names:
        anchor = habitability_symbols.equipment_anchor(name)
        assert anchor != habitability_symbols.EQUIPMENT_FALLBACK, name
        assert f'<tr id="{anchor}">' in page, name
    assert sorted(habitability_ref.EQUIPMENT_ROW_IDS.values()) == sorted(anchor for anchor, _ in habitability_symbols.EQUIPMENT_ANCHORS)
    assert set(habitability_ref.EQUIPMENT_ROW_IDS) == {label for label, _when in habitability_ref.EQUIPMENT_ITEMS}


def test_the_equipment_search_tags_no_longer_carry_a_levels_link(app):
    from planetgen.web import searchpage
    from planetgen.web.searchpage import SearchState

    options = [{"value": "0", "label": "Ideal", "count": 2}]
    with app.test_request_context("/search"):
        groups = searchpage.tag_groups(SearchState(), {"equipment": options})
    assert all("help" not in group for group in groups)


def test_a_range_gives_both_ends_the_one_style():
    """UX.107: a class page range never mixes plain and scientific ("5,000 to 1 × 10⁴ km")."""
    assert classref.number_range(5000, 1e4, "km") == "5 × 10³ to 1 × 10⁴ km"
    assert classref.number_range(100, 900) == "100 to 900"
    assert classref.number_range(0.5, 20000) == "0.5 to 2 × 10⁴"
