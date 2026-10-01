# html/web/class_pages.py

"""
The class reference pages, built from `lib/classref.py`'s catalog (made
once when the site starts, from the generator's own tables):

- `/classes`: every class type, with how many classes each has.
- `/classes/<type_slug>`: one type's classes, each linking to its page.
- `/classes/<type_slug>/<code>`: one class's facts.

An unknown type or code is a 404. `class_url` is what the system,
phenomenon and other pages use to link a class label here.
"""

from flask import abort

from classref import catalog, class_entry, class_type, class_url_parts

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


@bp.route("/classes")
def classes():
    """Every class type."""
    rows = [{
        "name": entry["name"],
        "url": page_url("class_type_page", type_slug=slug),
        "summary": entry["summary"],
        "count": len(entry["classes"]),
    } for slug, entry in catalog().items()]
    return render_page(
        "classes.html",
        title="Classes",
        section="classes",
        breadcrumbs=[crumb("Classes")],
        description="Every class of star, planet, nebula and other object this galaxy generator assigns.",
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
        description=entry["summary"],
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
        description=item["summary"],
        class_type=entry,
        type_url=page_url("class_type_page", type_slug=type_slug),
        item=item,
    )
