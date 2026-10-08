# planetgen/web/tables.py

"""
The data tables' JSON route (UX.41) and the helper list pages call to draw
one. A page module describes its table with `lib/datatable.Table`, calls
`register` on it, and renders it with `render`; `static/datatable.js` then
fetches further pages from `/table/<name>` as the visitor scrolls, sorts or
filters. See `lib/datatable.py` for the whole picture.
"""

from flask import jsonify, request, url_for

from planetgen.util import log
from planetgen.web.lib import apiclient
from planetgen.web.lib.datatable import parse_state, view
from planetgen.web.lib.pagination import PAGE_SIZE, clamp_page, page_offset

from . import bp

TABLES = {}
"""dict[str, Table]: Every registered table by name."""


def register(table):
    """Add `table` to the registry (its JSON route is then live) and return it."""
    TABLES[table.name] = table
    return table


def render(table, path, anchor=None):
    """
    Load the page of `table` the current request asks for and return its
    view model (`lib/datatable.view`) for `partials/datatable.html`. A page
    past the end shows the last real one.
    """
    state = parse_state(table, request.args)
    # Carry along only the other tables' parameters, never whatever else is in the address.
    others = set().union(*(other.owned_params() for other in TABLES.values())) - table.owned_params()
    extra = [(key, value) for key, value in request.args.items(multi=True) if key in others]
    result = table.load(state, PAGE_SIZE, page_offset(state.page), True)
    page = clamp_page(state.page, result.total)
    if page != state.page:
        state = state._replace(page=page)
        result = table.load(state, PAGE_SIZE, page_offset(page), True)
    return view(table, state, result, path, url_for("web.table_data", name=table.name), anchor=anchor, extra=extra)


def _whole(raw, default, lowest):
    """`raw` as an int of at least `lowest`, or `default` when it is missing or not a number."""
    try:
        return max(lowest, int(raw))
    except (TypeError, ValueError):
        return default


def _json_error(message, status):
    response = jsonify({"error": message})
    response.status_code = status
    return response


@bp.route("/table/<name>")
def table_data(name):
    """
    One page of a data table as JSON, for `static/datatable.js`:
    `?sort=&order=&<filters>&offset=&limit=` (limit at most one page), and
    `facets=1` to add the filter menus' option counts. Answers
    `{"rows", "total", "facets"}`; rows are lists of cells (see
    `lib/datatable.py`). An unknown table is a 404 and an API failure a 502,
    both as `{"error": ...}`.
    """
    table = TABLES.get(name)
    if table is None:
        return _json_error("No such table.", 404)
    state = parse_state(table, request.args, prefix="")
    limit = min(_whole(request.args.get("limit"), PAGE_SIZE, 1), PAGE_SIZE)
    offset = _whole(request.args.get("offset"), 0, 0)
    try:
        result = table.load(state, limit, offset, request.args.get("facets") == "1")
    except apiclient.ApiError as exc:
        log.exception(f"API error while loading the {name} table: {exc}")
        return _json_error("The table could not be loaded. Please try again shortly.", 502)
    response = jsonify({"rows": result.rows, "total": result.total, "facets": result.facets})
    response.headers["Cache-Control"] = "no-store"
    return response


table_data.json_only = True  # not a page: tests/test_web_a11y.py skips it
