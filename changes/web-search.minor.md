### Changed
- **The search page moved to `/search` (Flask), with bookmarkable GET
  URLs.** Every filter is a query parameter: `q` (the header search box)
  searches sector, system, star, planet and moon names at once; the
  per-object name fields (`sector_q`, ...), size ranges
  (`planet_min_radius_km`, ...), repeated tag facets
  (`spectral=G&spectral=K`) and each result panel's page
  (`stars_page=2`) follow it. Tags and "remove filter" chips are plain
  links, the per-object fields fold into a "Search by object and size"
  section, and results now appear above the tag browser. A submitted
  form's empty fields are dropped by a redirect to the short URL.
  `search.py` is now a shim that 301-redirects to `/search`, keeping every
  search parameter from an old link, bookmark or form post.
