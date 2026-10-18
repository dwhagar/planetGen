**Headline:** The system page and the one-off system page now show every body's detail from the system's description, not from Markdown.

### Changed
- The system page reads `GET /api/systems/<id>/description` and splits the Document per body (`describe.model.split_by_body`); each row, the star notes and the overview are written as HTML by `describe.html` (`blocks_html`), with a body's name linking to its row. `GET /api/systems/<id>/sections` and `db.render.render_system_sections` are removed, and `/description` now answers printed object IDs in every `ref`.
- The one-off system page (`/admin/generate/system`) shows the HTML page, with a Download button for Markdown, wikitext and plain text instead of a Markdown code box and preview. The command line gains `--format document` (the Document as JSON), which the page runs.
- Neither page module calls `markdown_to_html` any more.

TODO items resolved: UX.120, UX.105 / Added: none
