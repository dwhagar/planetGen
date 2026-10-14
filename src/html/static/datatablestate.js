// html/static/datatablestate.js
//
// The data table's logic without the page (UX.41): the query a state turns into, which pages the visible rows
// need, and the store that keeps the pages already fetched. datatable.js
// wires these to the DOM, TanStack Table and TanStack Virtual; plain node
// can test this file alone.

// Rows the server sends per page (lib/pagination.PAGE_SIZE).
export var PAGE_SIZE = 50;

// The query pairs for a state: `sort`/`order` (after the table's `prefix`,
// for a page with several tables), one pair per chosen filter value.
// `defaultSort` leaves the sort out when it is the page's default (so the
// address bar stays short).
export function stateParams(state, defaultSort, prefix) {
  var pairs = [];
  prefix = prefix || "";
  if (state.sort && (state.sort !== defaultSort || state.descending)) {
    pairs.push([prefix + "sort", state.sort]);
  }
  if (state.descending) {
    pairs.push([prefix + "order", "desc"]);
  }
  Object.keys(state.filters || {}).forEach(function (param) {
    (state.filters[param] || []).forEach(function (value) {
      pairs.push([param, value]);
    });
  });
  return pairs;
}

function encode(pairs) {
  return pairs.map(function (pair) {
    return encodeURIComponent(pair[0]) + "=" + encodeURIComponent(pair[1]);
  }).join("&");
}

// The table route's query for one page of rows: the state, the `offset` and
// `limit`, `facets=1` when the filter menus' counts are wanted, and `after`
// (the previous page's `next` key) when it is known, so a page far down the
// list is read from there instead of by skipping rows (PERF.75).
export function dataQuery(state, offset, limit, withFacets, after) {
  var pairs = stateParams(state, null);
  pairs.push(["offset", offset], ["limit", limit]);
  if (after) {
    pairs.push(["after", after]);
  }
  if (withFacets) {
    pairs.push(["facets", "1"]);
  }
  return encode(pairs);
}

// The address the page shows for a state: the page's own path, the query
// parameters that aren't this table's (`owned` names the ones that are:
// its sort, order, page and filters) and the same pairs a visitor would get
// from the no-script links. `search` is the address's current query string.
export function pageAddress(path, search, state, defaultSort, prefix, owned) {
  var pairs = [];
  new URLSearchParams(search).forEach(function (value, key) {
    if (owned.indexOf(key) < 0) {
      pairs.push([key, value]);
    }
  });
  var query = encode(pairs.concat(stateParams(state, defaultSort, prefix)));
  return path + (query ? "?" + query : "");
}

// The 0-based pages that hold rows `first` to `last` (inclusive).
export function pagesFor(first, last, pageSize) {
  var pages = [];
  for (var page = Math.floor(first / pageSize); page <= Math.floor(last / pageSize); page += 1) {
    pages.push(page);
  }
  return pages;
}

// How long to wait before asking again after the database was too busy to answer (a 504) `attempt` times in a
// row: 2 s, 4 s, 8 s ... 30 s, the Galaxy Map's pause (UX.96).
export function busyPause(attempt) {
  return Math.min(30000, 2000 * Math.pow(2, Math.max(1, Math.min(attempt, 5)) - 1));
}

// Fetches with `request()` (a promise of a Response), asking again after `busyPause` while the answer is a 504,
// for as long as `wanted()` says the page is still wanted. `onBusy(true)` is called when a retry begins waiting
// and `onBusy(false)` once an answer other than a 504 arrives (or the page is no longer wanted), so the page can
// show "database busy, trying again". Any other failing status throws as before.
export function fetchWhileBusy(request, wanted, onBusy, sleep) {
  var wait = sleep || function (ms) {
    return new Promise(function (resolve) { setTimeout(resolve, ms); });
  };
  function attempt(count) {
    return request().then(function (response) {
      if (response.status === 504 && wanted()) {
        onBusy(true);
        return wait(busyPause(count + 1)).then(function () {
          if (wanted()) {
            return attempt(count + 1);
          }
          onBusy(false);
          return Promise.reject(new Error("table no longer wanted"));
        });
      }
      onBusy(false);
      if (!response.ok) {
        throw new Error("table " + response.status);
      }
      return response.json();
    });
  }
  return attempt(0);
}

// The rows fetched so far, by page. `fetchPage(page, generation)` returns a
// promise of `{rows, total, facets, next}`; its third argument is the `next`
// of the page before, when that page was fetched. `reset()` forgets everything (the
// sort or a filter changed) and makes every answer still on its way stale,
// so a slow old answer can't land in the new table.
export function createPageStore(pageSize, fetchPage) {
  var pages = new Map();
  var keys = new Map();
  var pending = new Map();
  var generation = 0;
  var store = {
    total: 0,
    capped: false,
    facets: null,
    reset: function () {
      generation += 1;
      pages = new Map();
      keys = new Map();
      pending = new Map();
      store.total = 0;
      store.capped = false;
      store.facets = null;
    },
    // True when an answer fetched for `generation` is no longer wanted (the table was reset since).
    stale: function (mine) {
      return mine !== generation;
    },
    row: function (index) {
      var page = pages.get(Math.floor(index / pageSize));
      return page ? page[index % pageSize] : undefined;
    },
    has: function (page) {
      return pages.has(page);
    },
    // Resolves true when the page is now loaded (or already was), false when
    // the answer was stale or failed.
    ensure: function (page) {
      if (pages.has(page)) {
        return Promise.resolve(true);
      }
      if (pending.has(page)) {
        return pending.get(page);
      }
      var mine = generation;
      var request = fetchPage(page, mine, keys.get(page - 1)).then(function (body) {
        if (mine !== generation) {
          return false;
        }
        pages.set(page, body.rows);
        if (body.next) {
          keys.set(page, body.next);
        }
        store.total = body.total;
        store.capped = Boolean(body.capped);
        if (body.facets) {
          store.facets = body.facets;
        }
        return true;
      }, function () {
        return false;
      }).then(function (ok) {
        if (mine === generation) {
          pending.delete(page);
        }
        return ok;
      });
      pending.set(page, request);
      return request;
    },
  };
  return store;
}
