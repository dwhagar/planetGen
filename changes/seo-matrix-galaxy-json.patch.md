**Headline:** The search-engine page matrix lists the two newest galaxy JSON endpoints on purpose.

### Fixed
- **`web.galaxy_system_nearest` and `web.galaxy_regions` were in no search-engine class.** Both are public utility JSON for the galaxy map, so they join the other galaxy utility endpoints in `OTHER_PUBLIC_UTILITY` and `test_every_web_page_is_classified_on_purpose` passes again.

TODO items resolved: TEST.429 (the web_seo matrix part).
