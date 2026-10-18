**Headline:** A system's description can now be written as HTML tables, ready for the system page to use in the next change.

### Changed
- New `planetgen/describe/html.py`, the fourth renderer beside text, Markdown and wikitext: the Document's facts, ratings, composition and prose as HTML, every ladder quantity a `qty` span so the reader's unit choice applies, every piece of text escaped, and references linked through a function the page supplies. Tier names read in words and colours (`tier-blue` to `tier-red`); the six symbols arrive with UX.101. No page uses it yet (UX.120).
- `style.css` styles it from the site's colour tokens, so both themes and the `data-theme` override follow, and on a phone its tables scroll inside their section. `tests/test_describe_html.py` snapshots every block type and checks escaping; `tests/test_web_browser_document.py` checks contrast in the four theme combinations and that nothing runs past a 390 px screen.

TODO items resolved: UX.119 / Added: none
