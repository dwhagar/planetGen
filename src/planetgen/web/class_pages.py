# planetgen/web/class_pages.py

"""
The class reference pages, built from `planetgen/web/lib/classref.py`'s catalog (made
once when the site starts, from the generator's own tables):

- `/classes`: every class type, with how many classes each has.
- `/classes/<type_slug>`: one type's classes, each linking to its page.
- `/classes/<type_slug>/<code>`: one class's facts.

The habitability page (`/help/habitability`) is under the Help menu with the class pages.

An unknown type or code is a 404. `class_url` is what the system,
phenomenon and other pages use to link a class label here.
"""

from flask import abort
from markupsafe import Markup

from planetgen.web.lib import habitability_ref, life_help, property_help
from planetgen.web.lib.classref import alphabetical, catalog, class_entry, class_type, class_url_parts

from . import bp
from .helpers import crumb, page_url, population_status, render_page


def class_url(type_slug, code):
    """
    The URL of a class's page, or `None` when `code` is not a class of
    that type (an old or unknown value stays plain text). Needs a request
    or app context (`page_url`).
    """
    parts = class_url_parts(type_slug, code)
    if parts is None:
        return None
    return page_url("class_page", type_slug=parts[0], code=parts[1])


def _plain(summary):
    """A class summary as the page's plain-text description: one with formulas written out (UX.149) loses its tags."""
    return summary.striptags() if isinstance(summary, Markup) else summary


@bp.route("/classes")
def classes():
    """Every class type."""
    rows = [{
        "name": entry["name"],
        "url": page_url("class_type_page", type_slug=slug),
        "summary": entry["summary"],
        "count": len(entry["classes"]),
    } for slug, entry in catalog().items()]
    rows = alphabetical(rows)
    return render_page(
        "classes.html",
        title="Classes",
        section="classes",
        breadcrumbs=[crumb("Classes")],
        description="Every class of star, planet, nebula and other object this galaxy survey assigns.",
        rows=rows,
        population=_population_links(),
    )


def _population_links():
    """The Species and Polities pages, listed beside the classes once
    population data exists."""
    status = population_status()
    links = []
    if status["species"]:
        links.append({"label": "Species", "url": page_url("species")})
    if status["polities"]:
        links.append({"label": "Polities", "url": page_url("polities")})
    return links


@bp.route("/help/habitability")
def habitability_levels():
    """UX.90, UX.180: what the equipment levels and PHI-6 colours on the planet pages mean. The symbols, life scores
    and chips on every page link to its sections (UX.181)."""
    return render_page(
        "habitability.html",
        title="Habitability",
        section="help",
        breadcrumbs=[crumb("Habitability")],
        description="What the equipment levels, the Blue to Red domain colours and the Habitable chips mean.",
        **habitability_ref.page(),
    )


@bp.route("/help/properties")
def property_help_page():
    """UX.187: where the physics values on the star, planet and moon pages come from: each value's unit, who has it and the
    one-line source its page shows under it."""
    return render_page(
        "property_help.html",
        title="Properties",
        section="help",
        breadcrumbs=[crumb("Properties")],
        description="Where the orbit, heat, spin, field, ring and life values on a star, planet or moon page come from.",
        **property_help.page(),
    )


def _life_help_href(topic, number=None):
    """The address of a life help page (`life_help.page`'s `href`)."""
    if number is not None:
        return page_url("life_help_class", number=number)
    return page_url("life_help_page") if topic == "kinds" else page_url("life_help_topic", topic=topic)


def _render_life_help(topic, number=None):
    built = life_help.page(topic, number, _life_help_href)
    if built is None:
        abort(404)
    return render_page(
        "life_help.html",
        title=built["title"],
        section="help",
        breadcrumbs=[crumb("Life", "life_help_page")] + ([crumb(built["title"])] if built["topic"] != "kinds" else []),
        description=built["description"],
        page=built,
    )


@bp.route("/help/life")
def life_help_page():
    """DOC.84: the four kinds of life (microbes, flora, fauna, sentient life), with links to the pages on the scales, the stages and the classes."""
    return _render_life_help("kinds")


@bp.route("/help/life/<topic>")
def life_help_topic(topic):
    """DOC.84: one page of the life help: `abundance`, `trend`, `security`, `suitability`, `stages` or `classes`."""
    if topic == "kinds":
        abort(404)
    return _render_life_help(topic)


@bp.route("/help/life/classes/<int:number>")
def life_help_class(number):
    """DOC.84: one life class (1 to 9): what it drinks, eats, breathes and gives off, where it can live and how sure we are."""
    return _render_life_help("classes", number)


@bp.route("/classes/<type_slug>")
def class_type_page(type_slug):
    """One type's table of classes."""
    entry = class_type(type_slug)
    if entry is None:
        abort(404)
    rows = [{
        "code": code,
        "url": page_url("class_page", type_slug=type_slug, code=code),
        "name": item["name"],
        "summary": item["summary"],
    } for code, item in entry["classes"].items()]
    return render_page(
        "class_type.html",
        title=entry["name"],
        section="classes",
        breadcrumbs=[crumb("Classes", "classes"), crumb(entry["name"])],
        description=_plain(entry["summary"]),
        class_type=entry,
        rows=rows,
    )


@bp.route("/classes/<type_slug>/<code>")
def class_page(type_slug, code):
    """One class's facts."""
    entry = class_type(type_slug)
    item = class_entry(type_slug, code)
    if item is None:
        abort(404)
    # A letter reads well after its label ("Planet class M"); a word code
    # ("jupiter_family") is shown by its name instead.
    title = f"{entry['label']} {code}" if len(code) <= 3 else item["name"]
    return render_page(
        "class.html",
        title=title,
        section="classes",
        breadcrumbs=[crumb("Classes", "classes"), crumb(entry["name"], "class_type_page", type_slug=type_slug),
                     crumb(title)],
        description=_plain(item["summary"]),
        class_type=entry,
        type_url=page_url("class_type_page", type_slug=type_slug),
        item=item,
    )
