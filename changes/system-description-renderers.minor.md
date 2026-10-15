### Added

- `GET /api/systems/<id>/description` returns a system's description as structured data (sections, fact and rating tables, compositions, prose, and every number as a typed quantity), built from the database rows on each request.
- A system can be written as plain text as well as Markdown or wikitext (`/api/systems/<id>/text?format=text|markdown|wikitext`). Markdown and wikitext keep the headings and the wiki's data templates; the habitability rating is a table, and amounts are written by the same unit ladders as the pages.

### Changed

- The text, Markdown and wikitext of a system, its wiki upload and the `--markdown` output of `planetgen system` come from that description instead of each class writing its own text. Some wording is shorter or reordered (a star's age now reads "13.8 Gyr").
