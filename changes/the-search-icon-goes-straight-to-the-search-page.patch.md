**Headline:** The search icon in the top bar goes straight to the Search page instead of dropping a box down.

### Changed
- **The search icon opens the Search page (UX.160, Boss's live observation).** The header's search control was an icon button that dropped a small name box down (UX.146). It is now a plain icon link to `/search`, whose own Name box and filters do the searching. `static/headersearch.js` and the dropdown's styles are gone. Tests in `test_web_header_menus.py`, `test_web_pages.py` and `test_web_search.py`.
TODO items resolved: UX.160. Added: none.
