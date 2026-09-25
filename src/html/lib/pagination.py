# html/lib/pagination.py

"""
The site's one pagination control: every paged table in `html/` renders
its pager through `render_pagination` below, so they all look and behave
the same -- a "Showing X-Y of Z" summary, then First/Prev, a window of
numbered pages (with `...` gaps out to page 1 and the last page),
Next/Last. At most `PAGE_SIZE` rows per page, site-wide.

Pages carry a 1-based page *number* (e.g. `sectors_page=3`), not a raw
row offset, as a `page.post_link` hidden field like every other
navigational parameter here (see `lib/page.py`'s module docstring); each
table on a page has its own page parameter, and every pager link carries
the page's other parameters through unchanged, so paging one table never
resets another.

Two ways a page feeds this, depending on where its rows come from:

  - Paged through the API (`GET /api/sectors`/`/api/systems`/
    `/api/phenomena`/`/api/search`): `fetch_page` asks for exactly one
    page's rows by `limit`/`offset`, clamping a stale/out-of-range page
    number back to the last real page.
  - Rows the page already needs in full for something else (a sector's
    systems, which its 3D map also plots; the Galaxy Map's placed
    sectors, which its per-Quadrant summary also counts): `page_slice`
    cuts the same page out of that list in memory, so only one page of
    table rows is ever rendered.
"""

from urllib.parse import urlencode

from fmt import esc, post_link

PAGE_SIZE = 50
"""int: Most rows any paged table shows at once, site-wide."""

_WINDOW = 2
"""int: Numbered page links shown either side of the current page."""


def parse_page(raw):
    """
    Parses a page-number nav param into a 1-based int -- defensively: a
    missing, non-numeric, zero or negative value (a stale or hand-edited
    request) falls back to page 1 rather than raising, since this only
    ever controls which rows are shown.
    """
    try:
        value = int(raw)
    except (TypeError, ValueError):
        return 1
    return max(1, value)


def page_count(total, page_size=PAGE_SIZE):
    """Number of pages `total` rows fill (at least 1, so an empty table
    is still "page 1 of 1")."""
    return max(1, -(-int(total) // page_size))


def clamp_page(page, total, page_size=PAGE_SIZE):
    """`page` pulled back into `1..page_count(total)`."""
    return min(max(1, page), page_count(total, page_size))


def page_offset(page, page_size=PAGE_SIZE):
    """Row offset of `page`'s first row."""
    return (page - 1) * page_size


def page_slice(items, page, page_size=PAGE_SIZE):
    """
    One page of an in-memory list.

    Returns:
        tuple[list, int]: `(rows, page)` -- `page` clamped to the list's
            real page range first, so a stale page number past the end
            shows the last page rather than an empty table.
    """
    page = clamp_page(page, len(items), page_size)
    start = page_offset(page, page_size)
    return items[start:start + page_size], page


def fetch_page(fetch, page, page_size=PAGE_SIZE):
    """
    One page of a paged API listing.

    Args:
        fetch (callable): `fetch(limit, offset)` returning a paginated
            envelope (`items`/`total`, see `docs/api.md`'s "Pagination").
        page (int): The requested 1-based page.

    Returns:
        tuple[dict, int]: `(envelope, page)`. A page past the end (the
            table shrank since the link was rendered) is fetched again
            as the last real page instead of coming back empty.
    """
    envelope = fetch(page_size, page_offset(page, page_size))
    clamped = clamp_page(page, envelope["total"], page_size)
    if clamped != page:
        envelope = fetch(page_size, page_offset(clamped, page_size))
    return envelope, clamped


def _page_numbers(page, last):
    """The page numbers to show as links, with `None` marking a `...`
    gap: always 1 and `last`, plus `_WINDOW` pages either side of
    `page`. A gap of exactly one page shows that page instead of `...`."""
    shown = {1, last} | set(range(max(1, page - _WINDOW), min(last, page + _WINDOW) + 1))
    numbers = []
    previous = 0
    for number in sorted(shown):
        if number - previous == 2:
            numbers.append(number - 1)
        elif number - previous > 2:
            numbers.append(None)
        numbers.append(number)
        previous = number
    return numbers


def render_pagination(action, params, page_param, page, total, page_size=PAGE_SIZE, anchor=None, label="Pages",
                      method="post"):
    """
    Builds the pager for one table.

    Args:
        action (str): The page's own script, e.g. `"browse.py"`.
        params (dict or list[tuple]): Every parameter the page needs to
            re-render as it is now (db, id, filters, other tables' page
            numbers) -- `post_link`'s own `params` shape. `page_param`
            itself is added per link, so leave it out.
        page_param (str): This table's page parameter, e.g.
            `"sectors_page"`.
        page (int): The current (already clamped) 1-based page.
        total (int): Total row count across every page.
        anchor (str, optional): An element id to land on after following
            a link -- the table's own panel, so paging a table lower on
            the page doesn't jump back to the top.
        label (str): Accessible name for the `<nav>` landmark, e.g.
            `"Sector pages"` (several pagers can share one page).
        method (str): `"post"` (the default, every CGI page) renders each
            link as a `post_link` form button; `"get"` renders plain
            `<a href="{action}?{params}&{page_param}=N#{anchor}">` links
            instead, for the Flask pages (`html/web/`), whose URLs are
            ordinary bookmarkable GET URLs. In GET mode `action` is the
            page's own URL path (e.g. `url_for("web.index")`).

    Returns:
        str: A `<nav class="pagination">` block, or `""` when every row
            already fits on one page.
    """
    if total <= page_size:
        return ""

    last = page_count(total, page_size)
    base = list(params.items()) if isinstance(params, dict) else list(params)
    base = [(key, value) for key, value in base if key != page_param]
    target = f"{action}#{anchor}" if anchor else action

    def _link(number, text, css_class="page-link", attrs=""):
        if method == "get":
            query = urlencode([(key, value) for key, value in base + [(page_param, number)]
                               if value is not None and value != ""])
            href = f"{action}?{query}" + (f"#{anchor}" if anchor else "")
            return f'<a class="{esc(css_class)}" href="{esc(href)}" {attrs}>{text}</a>'
        return post_link(target, base + [(page_param, number)], text, css_class=css_class, attrs=attrs)

    def _disabled(text):
        return f'<span class="page-link disabled" aria-disabled="true">{text}</span>'

    controls = []
    if page > 1:
        controls.append(_link(1, "&laquo; First", attrs='aria-label="First page"'))
        controls.append(_link(page - 1, "&lsaquo; Prev", attrs='aria-label="Previous page"'))
    else:
        controls.append(_disabled("&laquo; First"))
        controls.append(_disabled("&lsaquo; Prev"))
    for number in _page_numbers(page, last):
        if number is None:
            controls.append('<span class="page-gap" aria-hidden="true">&hellip;</span>')
        elif number == page:
            controls.append(f'<span class="page-link current" aria-current="page">{number:,}</span>')
        else:
            controls.append(_link(number, f"{number:,}", attrs=f'aria-label="Page {number}"'))
    if page < last:
        controls.append(_link(page + 1, "Next &rsaquo;", attrs='aria-label="Next page"'))
        controls.append(_link(last, "Last &raquo;", attrs='aria-label="Last page"'))
    else:
        controls.append(_disabled("Next &rsaquo;"))
        controls.append(_disabled("Last &raquo;"))

    first_row = page_offset(page, page_size) + 1
    last_row = min(page * page_size, total)
    summary = f"Showing {first_row:,}&ndash;{last_row:,} of {total:,}"
    return (
        f'<nav class="pagination" aria-label="{esc(label)}">'
        f'<span class="pagination-summary">{summary}</span>'
        f'<span class="pagination-controls">{"".join(controls)}</span>'
        "</nav>"
    )
