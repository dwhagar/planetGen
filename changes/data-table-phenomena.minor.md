### Added

- The Phenomena list is the site's first data table on TanStack Table and
  TanStack Virtual (UX.41). Click a column header to sort by it (click again
  to reverse), filter by type and by descriptor (the nebula class, remnant
  shape, rogue planet kind and so on) from the menus above the table, and
  scroll through every row: the next 50 arrive as they come into view, so
  the page never holds more than a few dozen rows. The address bar follows
  the sort and filters, so a reload or a shared link shows the same table.
  Without scripts the table still sorts and filters with plain links and a
  form, and pages with the usual pager.
- `GET /api/phenomena` takes `sort`, `order`, `type`, `descriptor`, `placed`
  and `facets=1` (the option counts for the filter menus).
