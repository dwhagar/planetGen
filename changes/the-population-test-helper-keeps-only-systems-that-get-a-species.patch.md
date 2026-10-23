**Headline:** The population tests' civilized-system helper keeps only systems the population pass gives a species.

### Fixed
- **`test_a_generated_galaxys_species_are_all_inside_their_limits` and `test_territories_only_skips_new_worlds` each failed once in a gate (TEST.296).** The `_civilized_system` helper of `test_population.py` kept any system with a life timeline. The population pass gives a species only to a planet at the technological-civilization stage whose limits can hold it (POP.14: none where undrinkable water rules it out), so about 1 system in 75 got none (2 of 150 measured) and a test that expected two new species, or one, saw one, or none. The helper now checks both conditions on the saved system and deletes a system the pass would leave without a species before it tries again. 25 runs of the two tests in a row pass at `-n 4`.
- **`test_the_line_is_the_estimate_only_output_for_the_same_values` failed on a machine without the downloaded data files (found by this gate).** The size line runs `planetgen galaxy --estimate-only` in a process of its own, which refuses when a dataset is missing (GEN.599) because conftest's "nothing is missing" does not reach a subprocess. The test now sends the form's own `disregard_missing_data` switch.

TODO items resolved: TEST.296. Added: none.
