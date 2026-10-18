**Headline:** The Systems table has a "Starting at" box that begins the list at a typed letter, name or sector without paging down to it.

### Changed
- **The Systems table jumps by value (PERF.75).** Under the name sort and the sector sort a "Starting at" box begins the list at the first name, or the first sector, from the typed text on: `M` starts at the first name from M, and under a descending sort the names that begin with M come first and the list runs down from them. The count line says how many systems are left ("30 systems from “Drifting 100”"), the page address carries the text (`systems_from`), the scroll bar and the pager count from it, and another sort drops it and hides the box. It works with scripts off as a plain form. The jump is one more condition on the same read, so a jump deep into millions of systems costs what the first page does.
- **`GET /api/systems` takes `from=<text>`** with the same meaning (`total` counts the systems from there; `offset` and `after` count from there). The other sorts ignore it. `docs/api.md` and `docs/design/slow-reads-and-timeouts.md` describe it.

TODO items resolved: PERF.75. Added: none.
