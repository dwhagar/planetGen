**Headline:** The search-engine page matrix lists the two newest galaxy JSON endpoints on purpose.

### Fixed
- **`web.galaxy_system_nearest`, `web.galaxy_regions` and `web.galaxy_near` were in no search-engine class.** All three are public utility JSON for the galaxy map, so they join the other galaxy utility endpoints in `OTHER_PUBLIC_UTILITY` and `test_every_web_page_is_classified_on_purpose` passes again.
- **The settings schema file was out of date, and the batch-fill work-units test gave one fill the other's command.** `config.schema.json` is regenerated, and the test runs each fill as its own command's run (TEST.446).

TODO items resolved: TEST.429 (the web_seo matrix part), TEST.446.
