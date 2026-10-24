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
from planetgen.web.errors import RETRY_AFTER_SECONDS
from planetgen.web.lib import apiclient, slowpages
from planetgen.web.lib.slowroute import RETRY_HEADER, percent
from planetgen.web.lib.datatable import parse_state, view
from planetgen.web.lib.pagination import PAGE_SIZE, clamp_page, page_offset

from . import bp

TABLES = {}
"""dict[str, Table]: Every registered table by name."""


def register(table):
    """Add `table` to the registry (its JSON route is then live) and return it."""
    TABLES[table.name] = table
    return table


def render(table, path, anchor=None, keep=(), **source_args):
    """
    Load the page of `table` the current request asks for and return its
    view model (`lib/datatable.view`) for `partials/datatable.html`. A page
    past the end shows the last real one. `source_args` ride on the table's
    JSON address (a table about one thing, such as a polity's systems, names
    it there). `keep` names further page parameters to carry in every link.
    """
    state = parse_state(table, request.args)
    # Carry along only the other tables' parameters, never whatever else is in the address.
    others = set().union(*(other.owned_params() for other in TABLES.values())) - table.owned_params()
    extra = [(key, value) for key, value in request.args.items(multi=True) if key in others or key in keep]
    result = table.load(state, PAGE_SIZE, page_offset(state.page), True)
    page = clamp_page(state.page, result.total)
    if page != state.page:
        state = state._replace(page=page)
        result = table.load(state, PAGE_SIZE, page_offset(page), True)
    return view(table, state, result, path, url_for("web.table_data", name=table.name, **source_args), anchor=anchor, extra=extra)


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
    `?sort=&order=&<filters>&offset=&limit=` (limit at most one page; a keyset table also takes `after=`, the `next`
    of the page before, PERF.75), and
    `facets=1` to add the filter menus' option counts. Answers
    `{"rows", "total", "capped", "next", "facets", "note"}`; rows are lists of cells (see
    `lib/datatable.py`). An unknown table is a 404, a read that hit the database's statement time limit a 504
    with `Retry-After` (and `"busy": true`, so the table keeps its rows and asks again), and any other API
    failure a 502, all as `{"error": ...}`. A slow read (`slowpages.cost`) is not held in the request: the answer is
    `202` with `{"pending": true, "job": {...}, "status_url", "percent"}` (the job's `pagejobs.status`), and the table
    polls `status_url` and asks again when the job has succeeded (UX.136). `X-Slow-Retry: 1` replaces a failed job.
    """
    table = TABLES.get(name)
    if table is None:
        return _json_error("No such table.", 404)
    state = parse_state(table, request.args, prefix="")
    limit = min(_whole(request.args.get("limit"), PAGE_SIZE, 1), PAGE_SIZE)
    offset = _whole(request.args.get("offset"), 0, 0)
    slowpages.take_notice()      # one left from another request on this thread isn't this table's
    try:
        with slowpages.allow_pending(retry=request.headers.get(RETRY_HEADER) == "1"):
            result = table.load(state, limit, offset, request.args.get("facets") == "1",
                                **({"after": request.args.get("after") or None} if table.keyset else {}))
    except slowpages.Pending as pending:
        job = pending.job
        response = jsonify({"pending": True, "job": job, "percent": percent(job),
                            "status_url": url_for("api.page_job_status_route", job_id=job["id"])})
        response.status_code = 202
        response.headers["Cache-Control"] = "no-store"
        response.headers["Retry-After"] = "2"
        return response
    except apiclient.ApiError as exc:
        if exc.status_code == 504:
            log.debug(f"The {name} table hit the statement time limit: {exc}")
            response = jsonify({"error": "The database is busy; trying again shortly.", "busy": True})
            response.status_code = 504
            response.headers["Retry-After"] = str(RETRY_AFTER_SECONDS)
            response.headers["Cache-Control"] = "no-store"
            return response
        log.exception(f"API error while loading the {name} table: {exc}")
        return _json_error("The table could not be loaded. Please try again shortly.", 502)
    response = jsonify({"rows": result.rows, "total": result.total, "capped": bool(result.capped), "next": result.next,
                       "facets": result.facets, "note": result.note})
    response.headers["Cache-Control"] = "no-store"
    return response


table_data.json_only = True  # not a page: tests/test_web_a11y.py skips it
