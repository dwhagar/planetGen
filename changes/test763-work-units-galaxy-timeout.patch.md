**Headline:** The work-units galaxy test has the hour it needs, so the suite's 900 seconds no longer ends it.

### Changed

- TEST.763: `test_a_galaxy_run_stores_each_of_its_stages_once_and_follows_its_sector_units` runs a whole 36-sector galaxy through the real queue (21 min 12 s on a 4-core lane box), and the suite's default 900 seconds with the hang guard (TEST.621) ended it with no assertion. It now carries its own `pytest.mark.timeout(3600)`, as `test_sky_hyg_calibration` and `test_prevalence_life_census` do. Test only.

generation-output: unchanged
generation-affects: none
