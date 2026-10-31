# planetgen/web/references_pages.py

"""
The References pages (DOC.88, References R2; design: `docs/design/interface/references.md`).

Help > References is the one section of the site that speaks out of universe: where each number came from, what is a stand-in, what was
fitted to a result and how sure each part is. Every other page keeps a quiet "Where this comes from" link to the entry here
(`helpers.reference_link`; the addresses are `lib/reference_urls.py`).

- `/references`: the start page: what the section is, the headline numbers and the tier and tweak legends.

The pages are static: no database, no script. Their numbers are read from `docs/references/*.json`, built at setup from the registries by
`scripts/build_reference_pages.py`. Every page opens with the Out of universe banner (`id="out-of-universe-banner"`) and the site header
takes another tint (`body.references-mode`), so a reader always knows they have left the story.
"""

from planetgen.web.lib import references_site

from . import bp
from .helpers import crumb, render_page


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
