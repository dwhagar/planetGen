**Headline:** The old per-class description text, the Markdown renderer and the `markdown` dependency are gone; every page and download now comes from the one description Document.

### Changed
- Removed the old description path: each body's `to_paragraph_list` and `__str__`, the `MARKDOWN` text, `web/lib/mdrender.py`, and the `markdown` package from `setup.py` and the lock files.
- The `markdown` recipe key of a system run is gone (`docs/api.md`, `docs/system-file-format.md`); Markdown, wikitext and text are written from the Document.
- Tests that read the old text now read the Document through `tests/render_support.py`; `test_mdrender.py` is removed.
- The docs no longer describe `mdrender`; the Document HTML renderer escapes every string it writes.

TODO items resolved: UX.109
