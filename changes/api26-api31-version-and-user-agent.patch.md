**Headline:** The API version is 7 and records the shell, band and wedge address change, and every outgoing request now sends the standard User-Agent (the page scripts send the API version header).

### Changed
- **`API_VERSION` is 7 (API.26).** The switch of every sector address to `shell`/`band`/`wedge` (2026-10-15) never raised the version, and versions 4, 5 and 6 were taken by later changes, so 7 carries it; `docs/reference/api.md` lists it and the route snapshot moves with it. No route that takes the old `ring`/`layer`/`slot` names is kept.
- **One User-Agent everywhere (API.31).** `planetgen.util.user_agent.user_agent()` (release, project URL, API version) was already used by the dataset downloads; the web pages' API client and `scripts/build_element_tables.py` now send it too.
- **The page scripts send the API version.** A browser script cannot set `User-Agent`, so `static/apiversion.js`, loaded first by every page, adds `X-PlanetGen-API-Version` to each `fetch` and `XMLHttpRequest` to the site's own `/api/`. The server keeps answering with the header (API.22).
- `tests/test_user_agent.py` pins the format (product tokens, release, URL comment, API version), that it follows both versions, that every outgoing path calls the one function and that no other module opens a URL; `tests/js/apiversion.test.mjs` covers the page script.

TODO items resolved: API.26, API.31. Added: none.
