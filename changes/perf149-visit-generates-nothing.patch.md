**Headline:** A visit to an uncharted sector generates nothing and queues nothing, and the old generate-on-visit entry point is gone.

### Changed
- **No page or API route generates a sector while it answers.** `ensure_sector_generated`, the "recalculate on visit" entry point nothing in the site or API called any more, is removed so it cannot come back; sectors are made by runs (`planetgen galaxy`, the Generate page, `regenerate_sector_cell` for an API job). Tests check that a visit to an uncharted cell changes no table and starts no run, and that no site or API module calls the functions that save a sector.
- An admin's visit queuing a one-sector run stays open as PERF.204.
