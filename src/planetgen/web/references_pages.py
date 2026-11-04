# planetgen/web/references_pages.py

"""
The References pages (DOC.88, References R2; design: `docs/design/interface/references.md`).

Help > References is the one section of the site that speaks out of universe: where each number came from, what is a stand-in, what was
fitted to a result and how sure each part is. Every other page keeps a quiet "Where this comes from" link to the entry here
(`helpers.reference_link`; the addresses are `lib/reference_urls.py`).

- `/references`: the start page: what the section is, the headline numbers and the tier and tweak legends.
- `/references/stand-ins`: the registry of every stand-in (DOC.90), filtered by tier, tweak level, section and words; `/references/stand-ins/<slug>`:
  the page of one stand-in; `/references/fitted`: the stand-ins at tweak levels L2 to L4 and the tuning targets.
- `/references/papers` (DOC.89): every paper, book, and data set the project cites, with its DOI, arXiv, and ADS links, where we use it (area, code,
  data, stand-ins, design notes), views by A to Z, area, field, year, and type, totals, and a search over title, authors, DOI, and ID. A source has
  the anchor `#paper-r0123`;
- `/references/data` (DOC.89): the data sets fetched at setup, the archives and catalogs, the standards and web pages, and the terms of use that
  `DATA-LICENSE.md` records.

The pages are static: no database, no script. Their numbers are read from `docs/references/*.json`, built at setup from the registries by
`scripts/build_reference_pages.py`. Every page opens with the Out of universe banner (`id="out-of-universe-banner"`) and the site header
takes another tint (`body.references-mode`), so a reader always knows they have left the story.
"""

from flask import abort, request

from planetgen.web.lib import help_registry, references_site, references_usage

from . import bp
from .helpers import crumb, render_page
from .lib.reference_urls import reference_path


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
        banner=references_site.BANNER,
        banner_text=references_site.BANNER_TEXT,
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
        breadcrumbs=[crumb("References", "references_page"), crumb("Papers")],
        description="Every paper, book, and data set the project cites, with its links and where it is used.",
        page=page,
        search_action=reference_path("papers"),
        banner=references_site.BANNER,
        banner_text=references_site.BANNER_TEXT,
    )


@bp.route("/references/data")
def references_data_page():
    """DOC.89: the data sets, archives, catalogs, standards, and their terms of use."""
    return render_page(
        "references_list.html",
        title="Data sets and archives",
        section="references",
        breadcrumbs=[crumb("References", "references_page"), crumb("Data sets and archives")],
        description="The data sets downloaded at setup, the archives and catalogs, the standards, and the terms of use.",
        page=references_site.data_page(references_site.load("data.json"), references_site.licence_text()),
        search_action="",
        banner=references_site.BANNER,
        banner_text=references_site.BANNER_TEXT,
    )


def _not_built(title, crumbs):
    return render_page(
        "references.html", title=title, section="references", breadcrumbs=crumbs, banner=references_site.BANNER,
        banner_text=references_site.BANNER_TEXT,
        page={"title": title, "intro": [], "sections": [{
            "id": "not-built", "heading": "The pages are not built", "table": None,
            "paragraphs": ["The numbers on this page are built at setup by scripts/build_reference_pages.py. Run it, or run update.sh, and reload."]}]})


@bp.route("/references/stand-ins")
def references_stand_ins():
    """DOC.90: the registry of every stand-in, with filters."""
    crumbs = [crumb("References", "references_page"), crumb("Stand-ins")]
    built = _built("stand-ins.json", "overview.json")
    if built is None:
        return _not_built("Stand-ins", crumbs)
    args = request.args
    page = references_site.standins(built[0], built[1], tier=args.get("tier") or None, tweak=args.get("tweak") or None,
                                    section=args.get("section") or None, text=(args.get("text") or "").strip()[:100] or None)
    return render_page(
        "references_stand_ins.html", title="Stand-ins", section="references", breadcrumbs=crumbs,
        description="Every value that stands in for a measurement, with its source tier and tweak level.",
        page=page, banner=references_site.BANNER, banner_text=references_site.BANNER_TEXT)


@bp.route("/references/stand-ins/<slug>")
def references_stand_in(slug):
    """DOC.90: the page of one stand-in."""
    built = _built("stand-ins.json", "overview.json")
    crumbs = [crumb("References", "references_page"), crumb("Stand-ins", "references_stand_ins")]
    if built is None:
        return _not_built("Stand-in", crumbs + [crumb("Stand-in")])
    page = references_site.standin(built[0], built[1], slug)
    if page is None:
        abort(404)
    return render_page(
        "references_stand_in.html", title=page["title"], section="references", breadcrumbs=crumbs + [crumb(page["title"])],
        description=f"Where stand-in {page['row']['id']} comes from, how sure it is, and what would replace it.",
        page=page, banner=references_site.BANNER, banner_text=references_site.BANNER_TEXT)


@bp.route("/references/fitted")
def references_fitted():
    """DOC.90: the stand-ins fitted to a result, and the tuning targets."""
    crumbs = [crumb("References", "references_page"), crumb("Fitted to a result")]
    built = _built("stand-ins.json", "overview.json", "fitted.json")
    if built is None:
        return _not_built("Fitted to a result", crumbs)
    return render_page(
        "references_fitted.html", title="Fitted to a result", section="references", breadcrumbs=crumbs,
        description="The values that were bent toward a result, and the published statistics they were tuned to.",
        page=references_site.fitted(built[0], built[2], built[1]), banner=references_site.BANNER, banner_text=references_site.BANNER_TEXT)


help_registry.register_category("references", "Sources", "Where the numbers come from, and how sure each part is.", order=1000)
help_registry.register_topic("references_page", "References", "references",
                             "The papers and data sets behind the numbers, and how sure each part is.", order=1000)
