# planetgen/web/references_pages.py

"""
The References pages (DOC.88, References R2; design: `docs/design/interface/references.md`).

Help > References is the one section of the site that speaks out of universe: where each number came from, what is a stand-in, what was
fitted to a result and how sure each part is. Every other page keeps a quiet "More information" link to the entry here
(`helpers.reference_link`; the addresses are `lib/reference_urls.py`).

- `/references`: the start page: what the section is, the headline numbers and the tier and tweak legends.
- `/references/stand-ins`: the registry of every stand-in (DOC.90), filtered by tier, tweak level, section and words; `/references/stand-ins/<slug>`:
  the page of one stand-in; `/references/fitted`: the stand-ins at tweak levels L2 to L4 and the tuning targets.
- `/references/papers` (DOC.89): every paper, book, and data set the project cites, with its DOI, arXiv, and ADS links, where we use it (area, code,
  data, stand-ins, design notes), views by A to Z, area, field, year, and type, totals, and a search over title, authors, DOI, and ID. A source has
  the anchor `#paper-r0123`;
- `/references/data` (DOC.89): the data sets fetched at setup, the archives and catalogs, the standards and web pages, and the terms of use that
  `DATA-LICENSE.md` records.
- `/references/speculative` (DOC.92): how sure the air chemistry of each life class is, with one section for each class.

The pages are static: no database, no script. Their numbers are read from `docs/references/*.json`, built at setup from the registries by
`scripts/build_reference_pages.py`. Every page opens with the Out of universe banner (`id="out-of-universe-banner"`) and the site header
takes another tint (`body.references-mode`), so a reader always knows they have left the story.
"""

from flask import abort, request

from planetgen.web.lib import help_registry, property_help, references_site, references_usage

from . import bp
from .helpers import crumb, render_page
from .lib.reference_urls import AREAS, reference_path


_BANNER = {"banner": references_site.BANNER, "banner_text": references_site.BANNER_TEXT}
"""dict: The Out of universe banner every page of the section opens with, as `render_page` takes it."""


def _crumbs(*names):
    """The breadcrumbs of a page of the section: References, then each of `names` (the last one plain, the others linking to their page)."""
    crumbs = [crumb("References", "references_page")]
    for name, route in names:
        crumbs.append(crumb(name, route) if route else crumb(name))
    return crumbs


def _built(*names):
    """The built data files in `names` (`docs/references/*.json`), or `None` where one is missing, so the page can say the build has not run."""
    found = [references_site.load(name) for name in names]
    return found if all(item is not None for item in found) else None


@bp.route("/references")
def references_page():
    """DOC.88: the start page of the References section."""
    return render_page(
        "references.html",
        title="References",
        section="references",
        breadcrumbs=[crumb("References")],
        description="Where the numbers come from: the papers and data sets, the stand-ins and how sure each part is.",
        page=references_site.overview(references_site.load("overview.json")),
        **_BANNER,
    )


@bp.route("/references/papers")
def references_papers_page():
    """DOC.89: the sources, with links, where we use each, views, totals, and a search (`?view=area&group=stars`, `?q=lodders`)."""
    page = references_site.papers_page(
        references_site.load("papers.json"), references_usage.uses(), view=request.args.get("view", "az"),
        group=request.args.get("group", ""), q=request.args.get("q", "")[:200])
    return render_page(
        "references_list.html",
        title="Papers",
        section="references",
        breadcrumbs=_crumbs(("Papers", None)),
        description="Every paper, book, and data set the project cites, with its links and where it is used.",
        page=page,
        search_action=reference_path("papers"),
        **_BANNER,
    )


@bp.route("/references/speculative")
def references_speculative_page():
    """DOC.92: how sure the air chemistry of each life class is, and which classes are speculative."""
    return render_page(
        "references.html",
        title="Speculative chemistry",
        section="references",
        breadcrumbs=_crumbs(("Speculative chemistry", None)),
        description="How sure the air chemistry of each life class is, and which classes are speculative.",
        page=references_site.speculative(references_site.load("speculative.json")),
        **_BANNER,
    )


@bp.route("/references/data")
def references_data_page():
    """DOC.89: the data sets, archives, catalogs, standards, and their terms of use."""
    return render_page(
        "references_list.html",
        title="Data sets and archives",
        section="references",
        breadcrumbs=_crumbs(("Data sets and archives", None)),
        description="The data sets downloaded at setup, the archives and catalogs, the standards, and the terms of use.",
        page=references_site.data_page(references_site.load("data.json"), references_site.licence_text()),
        search_action="",
        **_BANNER,
    )


def _not_built(title, crumbs):
    return render_page(
        "references.html", title=title, section="references", breadcrumbs=crumbs, **_BANNER,
        page={"title": title, "intro": [], "sections": [{
            "id": "not-built", "heading": "The pages are not built", "table": None,
            "paragraphs": ["The numbers on this page are built at setup by scripts/build_reference_pages.py. Run it, or run update.sh, and reload."]}]})


@bp.route("/references/stand-ins")
def references_stand_ins():
    """DOC.90: the registry of every stand-in, with filters."""
    crumbs = _crumbs(("Stand-ins", None))
    built = _built("stand-ins.json", "overview.json")
    if built is None:
        return _not_built("Stand-ins", crumbs)
    args = request.args
    page = references_site.standins(built[0], built[1], tier=args.get("tier") or None, tweak=args.get("tweak") or None,
                                    section=args.get("section") or None, text=(args.get("text") or "").strip()[:100] or None)
    return render_page(
        "references_stand_ins.html", title="Stand-ins", section="references", breadcrumbs=crumbs,
        description="Every value that stands in for a measurement, with its source tier and tweak level.",
        page=page, **_BANNER)


@bp.route("/references/stand-ins/<slug>")
def references_stand_in(slug):
    """DOC.90: the page of one stand-in."""
    built = _built("stand-ins.json", "overview.json")
    crumbs = _crumbs(("Stand-ins", "references_stand_ins"))
    if built is None:
        return _not_built("Stand-in", crumbs + [crumb("Stand-in")])
    page = references_site.standin(built[0], built[1], slug)
    if page is None:
        abort(404)
    return render_page(
        "references_stand_in.html", title=page["title"], section="references", breadcrumbs=crumbs + [crumb(page["title"])],
        description=f"Where stand-in {page['row']['id']} comes from, how sure it is, and what would replace it.",
        page=page, **_BANNER)


@bp.route("/references/fitted")
def references_fitted():
    """DOC.90: the stand-ins fitted to a result, and the tuning targets."""
    crumbs = _crumbs(("Fitted to a result", None))
    built = _built("stand-ins.json", "overview.json", "fitted.json")
    if built is None:
        return _not_built("Fitted to a result", crumbs)
    return render_page(
        "references_fitted.html", title="Fitted to a result", section="references", breadcrumbs=crumbs,
        description="The values that were bent toward a result, and the published statistics they were tuned to.",
        page=references_site.fitted(built[0], built[2], built[1]), **_BANNER)


@bp.route("/references/method/<area>")
def references_method(area):
    """DOC.91: how one area is worked out: the values of the body pages, the stand-ins and the sources of the area."""
    if area not in AREAS:
        abort(404)
    crumbs = _crumbs((references_usage.AREA_TITLES[area], None))
    built = _built("stand-ins.json")
    if built is None:
        return _not_built(references_usage.AREA_TITLES[area], crumbs)
    page = references_site.method_page(area, built[0], property_help.by_area()[area], references_usage.uses())
    return render_page(
        "references.html", title=page["title"], section="references", breadcrumbs=crumbs,
        description=f"How the {references_usage.AREA_TITLES[area].lower()} part of the galaxy is worked out, and what it rests on.",
        page=page, **_BANNER)


@bp.route("/references/gaps")
def references_gaps():
    """DOC.91: what we have and what we lack, by area."""
    crumbs = _crumbs(("What we have and lack", None))
    built = _built("stand-ins.json")
    if built is None:
        return _not_built("What we have and lack", crumbs)
    return render_page(
        "references.html", title="What we have and lack", section="references", breadcrumbs=crumbs,
        description="For each part of the galaxy, how much rests on a source and which values are weakest.",
        page=references_site.gaps_page(built[0], references_usage.uses()), **_BANNER)


@bp.route("/references/scorecard")
def references_scorecard():
    """DOC.91: how well known each area is, computed from the registries."""
    crumbs = _crumbs(("Scorecard", None))
    built = _built("stand-ins.json")
    if built is None:
        return _not_built("Scorecard", crumbs)
    return render_page(
        "references.html", title="Scorecard", section="references", breadcrumbs=crumbs,
        description="How well known each part of the galaxy is, counted from the stand-in registry and the sources.",
        page=references_site.scorecard_page(built[0], references_usage.uses()), **_BANNER)


help_registry.register_category("references", "Sources", "Where the numbers come from, and how sure each part is.", order=1000)
help_registry.register_topic("references_page", "References", "references",
                             "The papers and data sets behind the numbers, and how sure each part is.", order=1000)
