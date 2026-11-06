**Headline:** The sector gravity data route is now classed in the search-engine indexing matrix as a utility route, so it is never offered for indexing.

### Fixed
- **Indexing class of the sector gravity route (TEST.653).** `web.sector_gravity`, the JSON route behind the Sector Map gravity layer (MAP.168), had no entry in the indexing matrix of `web/seo.py`, so `test_every_web_page_is_classified_on_purpose` failed. It joins the other data routes of the maps (`web.sector_scene`, `web.system_scene`) in the public utility class: reachable, not listed in a sitemap, and marked not to be indexed.
