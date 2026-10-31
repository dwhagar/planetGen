**Headline:** Every paper, book, dataset, and archive in the reference list has one fixed citation ID, shown on the References pages.

### Changed
- **Citation IDs (DOC.100).** Each source in `docs/reference/references-data.json` gets an ID such as `R0123`, given once and never reused. A preprint that later gets a DOI keeps its ID. `docs/REFERENCES.md` and `.html` print the ID before every entry with the anchor `#paper-r0123`, and `web/lib/reference_urls.py` addresses the Papers page by it. The in-line citation style (a bracketed ID, `[R0123]`, after the words it supports) is written into the References design and the plain-words guide. Display only: no schema change.
