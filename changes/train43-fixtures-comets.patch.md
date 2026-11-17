**Headline:** The seeded test galaxy and the API link audit no longer depend on a star holding comets by chance.

### Fixed
- **Two test fixtures assumed comets are common (GEN.1775 made them rarer).** `rich_galaxy_support` adds one parabolic comet when the seeded galaxy's few comets are all elliptical (`comets.parabolic_mean_anomaly` was NULL in every row, so `test_every_column_is_written` failed), and `test_api_link_ids` asks for comets (`COMETS = True`) and fails with the count tried instead of running 400 systems for 43 minutes and indexing a planetless one.

generation-output: unchanged
generation-affects: none
