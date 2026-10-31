# planetgen/web/lib/reference_urls.py

"""
Where the References pages are and how a page points into them (DOC.87, References R1).

The site speaks in universe. Where a number came from, what is a stand-in, what was fitted to a result and how sure a part is are said
in one place: the References pages (`/references`, Help > References), which carry an Out of universe banner. Every other page keeps a
quiet "Where this comes from" link to the entry that answers it. This module is the one place that knows the layout, so a page never
writes a References address by hand: it calls `web.helpers.reference_link(kind, id)` (or `reference_path` here) and the build, the
crawl test (`tests/test_references_crawl.py`) and the pages themselves cannot drift apart.

The address of an entry is `/references/<page>` and, for one entry on a page, `#<anchor>`. The kinds:

* `overview`: `/references`, the start page (no id);
* `papers`: `/references/papers`; the id of a source is its citation ID (`R0123`, DOC.100), which gives `#paper-r0123`;
* `data`: `/references/data`; the id of a dataset, catalog or archive is its citation ID too, and gives `#paper-<id>`;
* `method`: `/references/method/<area>`; the id is the area (`atmospheres`) or the area and a value (`atmospheres/surface_pressure`),
  which gives `#<anchor of the value>`;
* `stand-in`: `/references/stand-ins` for the whole registry and `/references/stand-ins/<anchor of the id>` for one (`GR1` is
  `/references/stand-ins/gr1`);
* `fitted`: `/references/fitted`, with `#<anchor>` for one value;
* `speculative`: `/references/speculative`, with `#class-<id>` for one life class (`3a` gives `#class-3a`);
* `gaps`: `/references/gaps`, with `#<area>` for one area;
* `scorecard`: `/references/scorecard`, with `#<area>` for one area;
* `about`: `/references/about`, credits and version (no id).

`anchor(text)` is the one way an id becomes part of an address: lower case, every run of characters that is not a letter or a digit
becomes one hyphen (`10.1086/151682` is `10-1086-151682`). Standard library only.
"""

import re

AREAS = ("stars", "galaxy", "orbits", "bodies", "moons-belts", "atmospheres", "climate", "radiation", "rings-magnetism", "life", "appearance")
"""tuple: The areas of the method pages (`/references/method/<area>`): stars; the galaxy; orbits; rocky and giant bodies; moons and belts;
atmospheres; climate; radiation; rings and magnetism; life; appearance (color and maps). The one list the method pages, the gaps page and
the scorecard share."""

ROOT = "/references"
"""str: The address of the References start page; every other address begins with it and a slash."""

_RUN = re.compile(r"[^a-z0-9]+")


def anchor(text):
    """`text` as the part of an address that names one entry: lower case, each run of other characters one hyphen, none at the ends."""
    slug = _RUN.sub("-", str(text).strip().lower()).strip("-")
    if not slug:
        raise ValueError(f"no letters or digits to name an entry by: {text!r}")
    return slug


def _page(name):
    return lambda entry: (f"{ROOT}/{name}", f"#{anchor(entry)}" if entry is not None else "")


PAPER_ID = re.compile(r"R\d{4,}")
"""re.Pattern: A citation ID (DOC.100): `R` and at least four digits, given once to every paper, book, dataset and archive of the reference
list (`docs/reference/references-data.json`) and never reused. The one key a text, a page and a link use for a source."""


def is_paper_id(text):
    """True when `text` is the shape of a citation ID (it may still name no source; the build checks that)."""
    return bool(PAPER_ID.fullmatch(str(text)))


def paper_anchor(paper_id):
    """The fragment that names one source on the Papers page and in `docs/REFERENCES.md`: `paper-r0123` for `R0123`."""
    if not is_paper_id(paper_id):
        raise ValueError(f"{paper_id!r} is not a citation ID: an ID is R and at least four digits (R0123)")
    return f"paper-{paper_id.lower()}"


def _paper(entry):
    return (f"{ROOT}/papers", f"#{paper_anchor(entry)}") if entry is not None else (f"{ROOT}/papers", "")


def _data(entry):
    return (f"{ROOT}/data", f"#{paper_anchor(entry)}") if entry is not None else (f"{ROOT}/data", "")


def _method(entry):
    if entry is None:
        raise ValueError("a method entry names its area: 'atmospheres' or 'atmospheres/surface_pressure'")
    area, _slash, value = str(entry).partition("/")
    if area not in AREAS:
        raise ValueError(f"{area!r} is not an area of the method pages: the areas are {', '.join(AREAS)}")
    return f"{ROOT}/method/{area}", f"#{anchor(value)}" if value else ""


def _stand_in(entry):
    return (f"{ROOT}/stand-ins/{anchor(entry)}", "") if entry is not None else (f"{ROOT}/stand-ins", "")


def _speculative(entry):
    return f"{ROOT}/speculative", f"#class-{anchor(entry)}" if entry is not None else ""


def _whole(name):
    return lambda entry: (f"{ROOT}/{name}", "") if entry is None else _refuse(name)


def _refuse(name):
    raise ValueError(f"the {name} page has no entries to point at")


KINDS = {
    "overview": lambda entry: (ROOT, "") if entry is None else _refuse("overview"),
    "papers": _paper,
    "data": _data,
    "method": _method,
    "stand-in": _stand_in,
    "fitted": _page("fitted"),
    "speculative": _speculative,
    "gaps": _page("gaps"),
    "scorecard": _page("scorecard"),
    "about": _whole("about"),
}
"""dict: Each kind of entry and the function that gives its `(path, "#anchor" or "")` from an id (or `None` for the whole page)."""


def reference_path(kind, entry=None):
    """
    The address of a References entry: `reference_path("stand-in", "GR1")` is `/references/stand-ins/gr1`, `reference_path("papers",
    "R0123")` is `/references/papers#paper-r0123` and `reference_path("speculative", "3a")` is
    `/references/speculative#class-3a`. `entry` is `None` for the page as a whole. A kind that does not exist, or an id the kind cannot
    use, is a `ValueError`, so a typo in a template fails the page and the crawl test instead of leaving a dead link.
    """
    if kind not in KINDS:
        raise ValueError(f"{kind!r} is not a kind of References entry: the kinds are {', '.join(KINDS)}")
    path, fragment = KINDS[kind](entry)
    return path + fragment


def is_reference_path(path):
    """True when `path` is a References page or an address inside one."""
    return path == ROOT or path.startswith(ROOT + "/")
