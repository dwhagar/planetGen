# planetgen/web/lib/datatable.py

"""
The site's shared data table (UX.41): sortable columns, faceted filters and
a page of rows fed by the API, 50 at a time. Every list page that adopts it
describes its table once, as a `Table`, and the pieces fit together like so:

  - The page renders the first page of the table on the server
    (`partials/datatable.html`): sort links in the headers, a filter form of
    checkboxes, and the shared pager. That works with scripts off.
  - `static/datatable.js` then takes over the same markup: TanStack Table
    keeps the sorting and filter state, TanStack Virtual scrolls the rows,
    and the next 50 rows come from the table's JSON route (`/table/<name>`,
    `web/tables.py`) as they scroll into view. Sorting and filtering fetch
    again from the first row.

Both read the same query parameters: `sort`, `order` (`asc`/`desc`), one
repeatable parameter per facet, and `page` on the server side. A table with
a `prefix` (several tables on one page) puts it in front of `sort`, `order`
and `page`.

A table's `load(state, limit, offset, want_facets)` returns a `Result` from
the API; its rows are lists of cells, one per column, each a dict with
`text`, and optionally `href` (a link) and `muted` (shown as a dash or "None"
in italics). A cell may instead carry `parts`, a list of plain strings and
`{"text", "href"}` links (a Location cell's nearest-neighbor links), and
`map_target` (a key the Sector Map knows, which adds a "Show on map" button
after the text). No cell ever holds markup: the browser fills it with text
nodes and the page template escapes it.
"""

import html
from collections import namedtuple
from html.parser import HTMLParser
from urllib.parse import urlencode

from markupsafe import Markup

from planetgen.web.lib.pagination import PAGE_SIZE, parse_page, render_pagination
from planetgen.web.lib.tabledisplay import to_plain_text

MAX_FILTER_VALUES = 60
"""int: Most values one filter parameter may carry."""

MAX_VALUE_LENGTH = 80
"""int: Longest filter value kept."""

Column = namedtuple("Column", "key label sortable", defaults=(True,))
"""A column: its key (the API's `sort` name), header text and whether it sorts."""

Facet = namedtuple("Facet", "param label")
"""A filter menu: the query parameter its checked values travel in, and its name."""

Result = namedtuple("Result", "rows total facets", defaults=(None,))
"""One loaded page: `rows` (lists of cells), the `total` rows that pass the
filters and `facets` (`{param: [{"value", "label", "count"}]}`, or None when
not asked for)."""

State = namedtuple("State", "sort descending filters page")
"""What a visitor chose: the sort key, its direction, `{param: [values]}` and the page."""


class Table:
    """
    One data table.

    Args:
        name (str): Its id in `/table/<name>` and in the page.
        label (str): What the rows are, for the accessible names ("Phenomena").
        columns (list[Column]): Left to right; the first sortable column is the
            default sort.
        load (callable): `load(state, limit, offset, want_facets)` -> `Result`.
        facets (list[Facet]): The filter menus.
        prefix (str): Put before `sort`, `order` and `page` in the page's URL.
        default_sort (str, optional): The sort key when none is asked for
            (the first sortable column if omitted).
        noun (tuple[str, str]): Singular and plural for the count line.
    """

    def __init__(self, name, label, columns, load, facets=(), prefix="", noun=("row", "rows"), default_sort=None):
        self.name = name
        self.label = label
        self.columns = list(columns)
        self.load = load
        self.facets = list(facets)
        self.prefix = prefix
        self.noun = noun
        self.default_sort = default_sort or next(column.key for column in columns if column.sortable)

    def owned_params(self):
        """The page query parameters this table reads; every other one belongs to the page or
        another table and is carried along unchanged."""
        return {f"{self.prefix}sort", f"{self.prefix}order", f"{self.prefix}page",
                *(facet.param for facet in self.facets)}

    def sort_keys(self):
        return [column.key for column in self.columns if column.sortable]


def parse_state(table, args, prefix=None):
    """
    The `State` a request's query parameters ask for, defensively: an unknown
    sort key, order or page falls back to the default; filter values are
    trimmed and capped.

    Args:
        table (Table): The table.
        args (werkzeug.datastructures.MultiDict): `request.args`.
        prefix (str, optional): Overrides `table.prefix` (the JSON route takes none).
    """
    prefix = table.prefix if prefix is None else prefix
    sort = args.get(f"{prefix}sort")
    if sort not in table.sort_keys():
        sort = table.default_sort
    filters = {}
    for facet in table.facets:
        values = []
        for raw in args.getlist(facet.param):
            value = raw.strip()[:MAX_VALUE_LENGTH]
            if value and value not in values:
                values.append(value)
        filters[facet.param] = values[:MAX_FILTER_VALUES]
    return State(sort, args.get(f"{prefix}order") == "desc", filters, parse_page(args.get(f"{prefix}page")))


def state_params(table, state, page=None, sort=None, descending=None, extra=()):
    """
    The query pairs that put `state` (with optional overrides) back in a URL.
    A default sort order is left out; `page` is added only above 1; `extra`
    (the page's other parameters) is added last.
    """
    sort = state.sort if sort is None else sort
    descending = state.descending if descending is None else descending
    pairs = []
    if sort != table.default_sort or descending:
        pairs.append((f"{table.prefix}sort", sort))
    if descending:
        pairs.append((f"{table.prefix}order", "desc"))
    for facet in table.facets:
        pairs.extend((facet.param, value) for value in state.filters.get(facet.param, ()))
    if page is not None and page > 1:
        pairs.append((f"{table.prefix}page", page))
    return pairs + list(extra)


def _href(path, pairs, anchor):
    query = urlencode(pairs)
    return path + (f"?{query}" if query else "") + (f"#{anchor}" if anchor else "")


def view(table, state, result, path, source, anchor=None, extra=()):
    """
    Everything `partials/datatable.html` needs to draw the table.

    Args:
        table (Table): The table.
        state (State): What was asked for.
        result (Result): The loaded page, with facets.
        path (str): The page's own path (links and the filter form go there).
        source (str): The table's JSON route.
        anchor (str, optional): Element id the links land on.
        extra (list[tuple]): The page's other query pairs, kept in every link.
    """
    columns = []
    for column in table.columns:
        header = {"key": column.key, "label": column.label, "sortable": column.sortable, "aria_sort": None}
        if column.sortable:
            current = column.key == state.sort
            header["aria_sort"] = ("descending" if state.descending else "ascending") if current else "none"
            header["arrow"] = ("▼" if state.descending else "▲") if current else ""
            header["href"] = _href(path, state_params(table, state, sort=column.key,
                                                      descending=(not state.descending) if current else False,
                                                      extra=extra), anchor)
        columns.append(header)
    facets = []
    for facet in table.facets:
        checked = set(state.filters.get(facet.param, ()))
        options = [{"value": option["value"], "label": option["label"], "count": option["count"],
                    "checked": option["value"] in checked}
                   for option in (result.facets or {}).get(facet.param, ())]
        known = {option["value"] for option in options}
        options.extend({"value": value, "label": value, "count": 0, "checked": True}
                       for value in state.filters.get(facet.param, ()) if value not in known)
        facets.append({"param": facet.param, "label": facet.label, "options": options, "selected": len(checked)})
    return {
        "id": table.name,
        "label": table.label,
        "path": path,
        "source": source,
        "anchor": anchor,
        "columns": columns,
        "rows": result.rows,
        "total": result.total,
        "noun_one": table.noun[0],
        "noun_many": table.noun[1],
        "facets": facets,
        "filtered": any(state.filters.get(facet.param) for facet in table.facets),
        "sort": state.sort,
        "descending": state.descending,
        "default_sort": table.default_sort,
        "prefix": table.prefix,
        "sort_fields": state_params(table, state._replace(filters={}), extra=extra),
        "clear_href": _href(path, state_params(table, state._replace(filters={}), extra=extra), anchor),
        "page_size": PAGE_SIZE,
        "pager": Markup(render_pagination(
            path, state_params(table, state, extra=extra), f"{table.prefix}page", state.page, result.total,
            anchor=anchor, label=f"{table.label} pages")),
    }


def plain(markup):
    """A formatter's HTML (`<sup>`, `&sup3;`) as the plain text a table cell shows: the cell is
    filled with `textContent` in the browser, so it can't hold markup."""
    return html.unescape(to_plain_text(str(markup)))


class _PartsParser(HTMLParser):
    """Collects text and `<a href>` links; any other tag is dropped, its text kept."""

    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.parts = []
        self._href = None

    def handle_starttag(self, tag, attrs):
        if tag == "a":
            self._href = dict(attrs).get("href")

    def handle_endtag(self, tag):
        if tag == "a":
            self._href = None

    def handle_data(self, data):
        if self._href:
            self.parts.append({"text": data, "href": self._href})
        else:
            self.parts.append(data)


def parts_of(markup):
    """The `parts` of a cell from trusted link markup (`<a href>` and text, as the site's
    location formatters write it)."""
    parser = _PartsParser()
    parser.feed(str(markup))
    parser.close()
    return parser.parts


def in_memory(items, state, limit, offset, want_facets, sorts, facet_values, to_cells):
    """
    A `Result` for a table whose rows are a list the page already holds (the sector's
    Contents, a Quadrant's sectors) rather than something the API pages: filter by
    `state.filters`, order by `state.sort`, then slice.

    Args:
        items (list): The table's rows in any form.
        state (State): What was asked for.
        limit (int), offset (int): The slice.
        want_facets (bool): Also count each filter menu's options.
        sorts (dict): Sort key -> `fn(item, descending)`, a value to order by or `None` for
            "no value", which always lists last.
        facet_values (dict): Facet parameter -> `fn(item)`, the item's value there or `None`.
        to_cells (callable): `to_cells(item)` -> the row's cells.
    """
    def passes(item, skip=None):
        for param, fn in facet_values.items():
            chosen = state.filters.get(param)
            if param != skip and chosen and fn(item) not in chosen:
                return False
        return True

    kept = [item for item in items if passes(item)]
    key = sorts[state.sort]
    keyed = [(key(item, state.descending), item) for item in kept]
    present = sorted((pair for pair in keyed if pair[0] is not None), key=lambda pair: pair[0],
                     reverse=state.descending)
    ordered = [item for _, item in present] + [item for value, item in keyed if value is None]
    facets = None
    if want_facets:
        facets = {}
        for param, fn in facet_values.items():
            counts = {}
            for item in items:
                value = fn(item)
                if value is not None and passes(item, skip=param):
                    counts[value] = counts.get(value, 0) + 1
            facets[param] = [{"value": value, "label": value, "count": counts[value]}
                             for value in sorted(counts, key=str.casefold)]
    return Result([to_cells(item) for item in ordered[offset:offset + limit]], len(ordered), facets)
