// html/static/datatable.js
//
// The site's data table (UX.41), taking over the table a list page rendered
// (partials/datatable.html, web/lib/datatable.py). The page already shows
// its first rows, sort links and filter form without scripts; this file
// makes them live:
//
//   - TanStack Table (vendor/tanstack/table-core.min.js) holds the sort and
//     the filters, and builds the rows and cells that are on screen.
//   - TanStack Virtual (vendor/tanstack/virtual-core.min.js) works out which
//     rows are on screen, so the table scrolls through every row while only
//     a few dozen are in the page.
//   - The rows come from the table's JSON route, one server page (50 rows)
//     at a time, as they scroll into view. Sorting or filtering starts over
//     from the first row. datatablestate.js keeps the fetched pages.
//
// Cells are filled with textContent and href only, never innerHTML. The
// address bar follows the sort and the filters (history.replaceState), so a
// reload or a shared link gives the same table.

// Imported with this module's own `?v=` query, as galaxymap3d.js explains.
const VERSION_QUERY = new URL(import.meta.url).search;
const { createTable, getCoreRowModel } = await import(`./vendor/tanstack/table-core.min.js${VERSION_QUERY}`);
const { Virtualizer, elementScroll, observeElementOffset, observeElementRect } =
  await import(`./vendor/tanstack/virtual-core.min.js${VERSION_QUERY}`);
const { PAGE_SIZE, createPageStore, dataQuery, pageAddress, pagesFor } =
  await import(`./datatablestate.js${VERSION_QUERY}`);
const { formatNumber } = await import(`./numberformat.js${VERSION_QUERY}`);

const ROW_HEIGHT_ESTIMATE = 37;
const OVERSCAN = 8;
const SVG_NS = "http://www.w3.org/2000/svg";

function el(tag, className, text) {
  const node = document.createElement(tag);
  if (className) {
    node.className = className;
  }
  if (text !== undefined) {
    node.textContent = text;
  }
  return node;
}

// The "Show on map" button of a cell with a `map_target` (the sector page's
// Contents). Hidden until the Sector Map, which listens for the
// "datatable:rows" event, knows the target.
function mapButton(cell, iconUrl) {
  const button = el("button", "icon-btn");
  button.type = "button";
  button.hidden = true;
  button.title = "Show on map";
  button.setAttribute("aria-label", `Show ${cell.text} on the map`);
  button.dataset.mapTarget = cell.map_target;
  const svg = document.createElementNS(SVG_NS, "svg");
  svg.setAttribute("class", "icon");
  svg.setAttribute("aria-hidden", "true");
  svg.setAttribute("focusable", "false");
  const use = document.createElementNS(SVG_NS, "use");
  use.setAttribute("href", `${iconUrl}#show-on-map`);
  svg.append(use);
  button.append(svg);
  return button;
}

// One table cell from a served cell: `{text, href?, muted?, swatch?, parts?, map_target?, form?}`.
function cellNode(cell, iconUrl) {
  const td = el("td");
  if (!cell) {
    td.className = "datatable-pending";
    td.textContent = "…";
    return td;
  }
  if (cell.swatch) {
    // A polity's map color: a small square, as the page draws it with no scripts.
    const svg = document.createElementNS(SVG_NS, "svg");
    svg.setAttribute("class", "polity-swatch");
    svg.setAttribute("viewBox", "0 0 10 10");
    svg.setAttribute("aria-hidden", "true");
    svg.setAttribute("focusable", "false");
    const square = document.createElementNS(SVG_NS, "rect");
    square.setAttribute("width", "10");
    square.setAttribute("height", "10");
    square.setAttribute("fill", cell.swatch);
    svg.append(square);
    td.append(svg, " ");
  }
  if (cell.form) {
    // A POST button, as the admin's API keys Revoke; its fields carry the CSRF token.
    const form = el("form", "table-form");
    form.method = "post";
    form.action = cell.form.action;
    cell.form.fields.forEach(([name, value]) => {
      const input = el("input");
      input.type = "hidden";
      input.name = name;
      input.value = value;
      form.append(input);
    });
    const button = el("button", "btn", cell.form.button);
    button.type = "submit";
    button.setAttribute("aria-label", cell.form.label);
    form.append(button);
    td.append(form);
  } else if (cell.parts) {
    // Text and links, as a Location cell lists its nearest systems.
    cell.parts.forEach((part) => {
      if (typeof part === "string") {
        td.append(part);
      } else {
        const link = el("a", "", part.text);
        link.href = part.href;
        td.append(link);
      }
    });
    td.classList.add("inline-links");
  } else if (cell.href) {
    const link = el("a", "", cell.text);
    link.href = cell.href;
    if (cell.label) {
      link.setAttribute("aria-label", cell.label);
    }
    td.append(link);
  } else if (cell.muted) {
    td.append(el("em", "", cell.text));
  } else {
    td.append(cell.text);
  }
  if (cell.map_target) {
    td.append(" ", mapButton(cell, iconUrl));
  }
  return td;
}

function enhance(root) {
  const scroller = root.querySelector(".datatable-scroll");
  const tableEl = scroller.querySelector("table");
  const tbody = tableEl.querySelector("tbody");
  const headers = Array.from(tableEl.querySelectorAll("th[data-col]"));
  const countLine = root.querySelector(".datatable-count-line");
  const filterForm = root.querySelector(".datatable-filters");
  const facetNodes = Array.from(root.querySelectorAll(".datatable-facet"));
  const noun = (count) => (count === 1 ? root.dataset.nounOne : root.dataset.nounMany);
  const clearLink = root.querySelector(".datatable-clear");
  const source = root.dataset.source;
  const defaultSort = root.dataset.defaultSort;
  const prefix = root.dataset.prefix || "";
  const iconUrl = root.dataset.icons || "";
  const owned = [`${prefix}sort`, `${prefix}order`, `${prefix}page`]
    .concat(facetNodes.map((node) => node.dataset.param));
  let rowHeight = ROW_HEIGHT_ESTIMATE;

  const columns = headers.map((th, index) => ({
    id: th.dataset.col,
    header: th.textContent.trim(),
    accessorFn: (row) => (row.cells ? row.cells[index] : null),
    enableSorting: Boolean(th.querySelector(".datatable-sort")),
  }));

  // --- TanStack Table: the sort and the filters ------------------------
  const table = createTable({
    data: [],
    columns,
    state: {},
    onStateChange: () => {},
    renderFallbackValue: null,
    manualSorting: true,
    manualFiltering: true,
    manualPagination: true,
    autoResetAll: false,
    enableMultiSort: false,
    enableSortingRemoval: false,
    sortDescFirst: false,
    getRowId: (row) => String(row.index),
    getCoreRowModel: getCoreRowModel(),
  });
  let tableState = {
    ...table.initialState,
    sorting: [{ id: root.dataset.sort, desc: root.dataset.descending === "true" }],
    columnFilters: facetNodes.map((node) => ({ id: node.dataset.param, value: checkedValues(node) }))
      .filter((filter) => filter.value.length > 0),
  };
  table.setOptions((previous) => ({
    ...previous,
    state: tableState,
    onStateChange: (updater) => {
      tableState = typeof updater === "function" ? updater(tableState) : updater;
      table.setOptions((options) => ({ ...options, state: tableState }));
      stateChanged();
    },
  }));

  function checkedValues(node) {
    return Array.from(node.querySelectorAll("input[type=checkbox]:checked")).map((input) => input.value);
  }

  function currentState() {
    const sorting = tableState.sorting[0] || { id: defaultSort, desc: false };
    const filters = {};
    tableState.columnFilters.forEach((filter) => {
      filters[filter.id] = filter.value;
    });
    return { sort: sorting.id, descending: Boolean(sorting.desc), filters };
  }

  // --- Fetched pages -----------------------------------------------------
  let wantFacets = true;
  const store = createPageStore(PAGE_SIZE, (page, generation, after) => {
    const withFacets = page === 0 && wantFacets;
    return fetch(`${source}${source.includes("?") ? "&" : "?"}${dataQuery(currentState(), page * PAGE_SIZE, PAGE_SIZE, withFacets, after)}`,
      { headers: { Accept: "application/json" } }).then((response) => {
      if (!response.ok) {
        throw new Error(`table ${response.status}`);
      }
      return response.json();
    });
  });

  // --- TanStack Virtual: which rows are on screen -------------------------
  const virtualizer = new Virtualizer({
    count: Number(root.dataset.total) || 0,
    getScrollElement: () => scroller,
    estimateSize: () => rowHeight,
    overscan: OVERSCAN,
    scrollToFn: elementScroll,
    observeElementRect,
    observeElementOffset,
    onChange: () => scheduleRender(),
  });
  virtualizer._didMount();

  // The virtualizer reports a change for every row it measures, and rows of uneven
  // height (a Location cell that wraps) make those changes feed back into the render;
  // one render per frame lets them settle.
  let frame = 0;
  function scheduleRender() {
    if (!frame) {
      frame = requestAnimationFrame(() => {
        frame = 0;
        render();
      });
    }
  }

  let drawn = "";
  function render() {
    const items = virtualizer.getVirtualItems();
    // Drawing replaces the rows, and the new rows are measured again, which reports another
    // change; skipping a draw that would look the same ends that loop.
    const signature = `${virtualizer.getTotalSize()}|${table.options.state.sorting.map((s) => s.id + s.desc)}|` +
      items.map((item) => `${item.index}:${item.start}:${store.row(item.index) ? 1 : 0}`).join(",");
    if (signature === drawn) {
      return;
    }
    drawn = signature;
    const rows = items.map((item) => ({ index: item.index, cells: store.row(item.index) }));
    table.setOptions((options) => ({ ...options, data: rows }));
    if (items.length) {
      pagesFor(items[0].index, items[items.length - 1].index, PAGE_SIZE)
        .filter((page) => !store.has(page))
        .forEach((page) => store.ensure(page).then((loaded) => loaded && render()));
    }
    const nodes = [];
    const colSpan = columns.length;
    const spacer = (height) => {
      const tr = el("tr", "datatable-spacer");
      tr.setAttribute("aria-hidden", "true");
      const td = el("td");
      td.colSpan = colSpan;
      td.style.height = `${height}px`;
      tr.append(td);
      return tr;
    };
    if (items.length && items[0].start > 0) {
      nodes.push(spacer(items[0].start));
    }
    const built = table.getRowModel().rows.map((row, position) => {
      const tr = el("tr");
      const index = items[position].index;
      tr.dataset.index = String(index);
      tr.setAttribute("aria-rowindex", String(index + 2));
      if (!row.original.cells) {
        tr.className = "datatable-pending";
        tr.setAttribute("aria-busy", "true");
      }
      row.getVisibleCells().forEach((cell) => tr.append(cellNode(cell.getValue(), iconUrl)));
      return tr;
    });
    nodes.push(...built);
    if (items.length) {
      const tail = virtualizer.getTotalSize() - items[items.length - 1].end;
      if (tail > 0) {
        nodes.push(spacer(tail));
      }
    }
    if (!items.length) {
      const tr = el("tr");
      const td = el("td", "", "None");
      td.colSpan = colSpan;
      td.className = "datatable-none";
      tr.append(td);
      nodes.push(tr);
    }
    tbody.replaceChildren(...nodes);
    built.forEach((tr) => virtualizer.measureElement(tr));
    document.dispatchEvent(new CustomEvent("datatable:rows"));
  }

  // --- Keeping the page around the table in step -------------------------
  function showHeaders() {
    const state = currentState();
    headers.forEach((th) => {
      const link = th.querySelector(".datatable-sort");
      if (!link) {
        return;
      }
      const active = th.dataset.col === state.sort;
      th.setAttribute("aria-sort", active ? (state.descending ? "descending" : "ascending") : "none");
      link.querySelector(".datatable-arrow").textContent = active ? (state.descending ? "▼" : "▲") : "";
    });
  }

  function showCount() {
    const state = currentState();
    const filtered = Object.keys(state.filters).some((param) => state.filters[param].length > 0);
    const more = root.dataset.capped === "true" || store.capped ? "+" : "";
    countLine.textContent = `${formatNumber(store.total)}${more} ${noun(store.total + (more ? 1 : 0))}${filtered ? " match" : ""}`;
    if (clearLink) {
      clearLink.hidden = !filtered;
    }
    tableEl.setAttribute("aria-rowcount", String(store.total + 1));
  }

  function updateFacets(facets) {
    facetNodes.forEach((node) => {
      const options = (facets && facets[node.dataset.param]) || [];
      const list = node.querySelector(".datatable-options");
      const items = new Map(Array.from(list.children).map((li) => [li.querySelector("input").value, li]));
      options.forEach((option) => {
        let li = items.get(option.value);
        if (!li) {
          li = el("li");
          const label = el("label");
          const input = el("input");
          input.type = "checkbox";
          input.name = node.dataset.param;
          input.value = option.value;
          label.append(input, el("span", "datatable-option-label", option.label), document.createTextNode(" "),
            el("span", "datatable-count", ""));
          li.append(label);
          const after = Array.from(list.children).find((other) =>
            other.querySelector(".datatable-option-label").textContent.localeCompare(option.label) > 0);
          list.insertBefore(li, after || null);
          items.set(option.value, li);
        }
        li.querySelector(".datatable-count").textContent = `(${option.count})`;
        li.hidden = false;
      });
      const present = new Set(options.map((option) => option.value));
      items.forEach((li, value) => {
        if (!present.has(value)) {
          li.querySelector(".datatable-count").textContent = "(0)";
          li.hidden = !li.querySelector("input").checked;
        }
      });
      const chosen = checkedValues(node).length;
      const summary = node.querySelector("summary");
      let note = summary.querySelector(".datatable-selected");
      if (chosen && !note) {
        note = el("span", "datatable-selected");
        summary.append(document.createTextNode(" "), note);
      }
      if (note) {
        note.textContent = chosen ? `(${chosen} chosen)` : "";
      }
    });
  }

  function showAddress() {
    try {
      history.replaceState(history.state, "", pageAddress(location.pathname, location.search, currentState(),
        defaultSort, prefix, owned) + location.hash);
    } catch (error) {
      // The address bar is only a convenience.
    }
  }

  // Everything has to be fetched again from the first row.
  let loadToken = 0;
  function reload() {
    loadToken += 1;
    const mine = loadToken;
    store.reset();
    drawn = "";
    wantFacets = true;
    root.setAttribute("aria-busy", "true");
    return store.ensure(0).then((loaded) => {
      if (mine !== loadToken || !loaded) {
        return;
      }
      root.removeAttribute("aria-busy");
      virtualizer.setOptions({ ...virtualizer.options, count: store.total });
      virtualizer._willUpdate();
      scroller.scrollTop = 0;
      virtualizer.scrollToOffset(0);
      updateFacets(store.facets);
      showHeaders();
      showCount();
      render();
    });
  }

  // TanStack Table also reports state this table doesn't use; only a new
  // sort or filter fetches again.
  let shown = JSON.stringify(currentState());
  function stateChanged() {
    const now = JSON.stringify(currentState());
    if (now === shown) {
      return;
    }
    shown = now;
    showAddress();
    reload();
  }

  // --- Events ---------------------------------------------------------------
  tableEl.querySelector("thead").addEventListener("click", (event) => {
    const link = event.target.closest && event.target.closest(".datatable-sort");
    if (!link) {
      return;
    }
    event.preventDefault();
    table.getColumn(link.closest("th").dataset.col).toggleSorting(undefined, false);
  });

  function filtersFromForm() {
    table.setColumnFilters(facetNodes.map((node) => ({ id: node.dataset.param, value: checkedValues(node) }))
      .filter((filter) => filter.value.length > 0));
  }

  if (filterForm) {
    filterForm.addEventListener("change", filtersFromForm);
    filterForm.addEventListener("submit", (event) => {
      event.preventDefault();
      filtersFromForm();
    });
    if (clearLink) {
      clearLink.addEventListener("click", (event) => {
        event.preventDefault();
        filterForm.querySelectorAll("input[type=checkbox]").forEach((input) => {
          input.checked = false;
        });
        filtersFromForm();
      });
    }
  }

  // The first fetch decides it: the rows the server rendered stay until the
  // table has its own, so a failed fetch leaves a working static table.
  store.ensure(0).then((loaded) => {
    if (!loaded) {
      return;
    }
    root.classList.add("datatable-enhanced");
    root.dataset.enhanced = "true";
    virtualizer.setOptions({ ...virtualizer.options, count: store.total });
    virtualizer._willUpdate();
    updateFacets(store.facets);
    showCount();
    render();
    const first = tbody.querySelector("tr[data-index]");
    if (first) {
      rowHeight = first.getBoundingClientRect().height || ROW_HEIGHT_ESTIMATE;
    }
  });
}

document.querySelectorAll("[data-datatable]").forEach(enhance);
