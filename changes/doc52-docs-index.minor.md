**Headline:** The docs have one index page that lists every HTML page with what it contains and links back from each.

### Changed
- **`docs/index.html` (DOC.52).** One entry per HTML page under `docs/` (the execution plan, the version plan, the TODO tree and reference, the release roadmap, the architecture page and the architecture map) with its title, a one- or two-sentence description of what it holds and who it is for, and a relative link, so it opens offline. It also lists the reference notes (`docs/*.md`, each with its first sentence), the design notes' index, the goal issues' note and the GitHub issue list. Each page now carries a `<meta name="description">` and a "Docs index" link in its header, in the TODO page templates too.
- **`scripts/build_docs_index.py` writes it from the pages themselves**, so a new page cannot be forgotten: `src/tests/test_docs_index.py` fails with the command to run when a page is missing from the index, has no description or has no link back, or when the index is out of date. Run `python scripts/build_docs_index.py` after adding a page.

TODO items resolved: DOC.52. Added: none.
