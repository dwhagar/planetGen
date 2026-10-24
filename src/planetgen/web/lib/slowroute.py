# planetgen/web/lib/slowroute.py

"""
The loading state of a slow page (UX.118; design in docs/design/operations/slow-page-answers.md).

`@slow_route(...)` on a page's view lets the view's slow requests (`slowpages.cost`) wait for a job instead of
holding the request (PERF.105). When one does, the view is abandoned and the visitor gets the same page frame
(`title`, breadcrumbs, header) with a progress bar in place of the content, status `202`, never cached or indexed.
`static/slowload.js` polls `GET /api/page-jobs/<id>`; when the job has finished it fetches the page again with
the same address and swaps the finished content in without reloading, and if that page is waiting for another job
(a search asks for several answers) it carries on with that one. A failed or cancelled job shows its reason and a
Try again, which asks for the page with the `X-Slow-Retry` header: a new job replaces the failed one. Without
JavaScript the page refreshes itself every few seconds (`<noscript>`).

A view that built a slow answer in the request because the queue could not take it (`slowpages.take_notice`) says
so at the top of its page (`slow_notice`, a template global: `slowpages.take_notice`).

The cache warm-up (PERF.104, `config.WARMUP`) is not a visitor: it visits to fill the page cache, so its requests
never wait for a job (PERF.156). A slow page is built in the warm-up's own request, as it always was, and the unit
that visited it has the answer cached when it ends; waiting on `202` would count the loading page as a failed visit.
"""

import functools

from flask import current_app, request

from planetgen.util.format import format_duration_seconds, format_number
from planetgen.web.helpers import render_page
from planetgen.web.lib import slowpages

RETRY_HEADER = "X-Slow-Retry"
"""str: The header the loading page's Try again sends: replace a failed or cancelled job at once."""

REFRESH_SECONDS = 5
"""int: How often the page refreshes itself without JavaScript."""

STATE_TEXT = {
    "queued": "Waiting for a worker to start building this page.",
    "running": "Building this page.",
    "failed": "This page could not be built.",
    "cancelled": "An administrator stopped building this page.",
    "succeeded": "Built. Showing it now.",
}
"""dict: What each job state says, in the page and (the same words) in `slowload.js`."""


def progress_text(job):
    """The line under the bar for `job` (`pagejobs.status`): an estimate, said so. `slowload.js` writes the same."""
    if job["state"] != "running" or job["total"] is None or job["done"] is None:
        return STATE_TEXT.get(job["state"], "")
    percent = 0 if not job["total"] else min(99, int(100 * job["done"] / job["total"]))
    left = ("still working" if job["eta_seconds"] is None
            else f"about {format_duration_seconds(job['eta_seconds'])} left")
    return (f"About {percent}% (an estimate: {format_number(job['done'])} of {format_number(job['total'])} rows), "
            f"{left}.")


def percent(job):
    """The bar's value 0 to 100 for `job`, or `None` when it can't say (queued, or its size is not known)."""
    if job["state"] == "succeeded":
        return 100
    if job["state"] != "running" or not job["total"] or job["done"] is None:
        return None
    return min(99, int(100 * job["done"] / job["total"]))


def loading_page(pending, **page):
    """The loading state for `pending` (`slowpages.Pending`) as a `render_page` result, in the page frame."""
    job = pending.job
    html, _status = render_page("slow_loading.html", status=202, job=job, state_text=STATE_TEXT,
                                progress_text=progress_text(job), percent=percent(job),
                                refresh_seconds=REFRESH_SECONDS, **page)
    return html, 202, {"Cache-Control": "no-store", "X-Robots-Tag": "noindex", "Retry-After": "2"}


def slow_route(title, section=None, breadcrumbs=(), description=None):
    """
    Decorator for a view whose answers may be slow to build.

    Args:
        title (str): The page's title, as `render_page` takes it.
        section (str, optional): The active section.
        breadcrumbs (callable or list): `crumb(...)` items, or a function that returns them (called when needed).
        description (str, optional): The page's description.
    """
    def decorate(view):
        @functools.wraps(view)
        def wrapper(*args, **kwargs):
            if current_app.config.get("WARMUP"):
                return view(*args, **kwargs)
            slowpages.take_notice()      # one left from another request on this thread isn't this page's
            try:
                with slowpages.allow_pending(retry=request.headers.get(RETRY_HEADER) == "1"):
                    return view(*args, **kwargs)
            except slowpages.Pending as pending:
                crumbs = breadcrumbs() if callable(breadcrumbs) else breadcrumbs
                return loading_page(pending, title=title, section=section, breadcrumbs=crumbs,
                                    description=description)
        return wrapper
    return decorate
