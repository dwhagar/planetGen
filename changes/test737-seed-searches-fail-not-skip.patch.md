**Headline:** No test skips quietly any more when its seed did not draw the body it needs.

### Changed
- **Seed-dependent skips are gone.** About 60 tests ended in `pytest.skip("no ... in these draws")` when an unseeded search found nothing, so a run could pass with the check never made. Each search now draws from fixed streams (`draw.bound(...)`), and a search that still finds nothing fails with a message that names how many fixed streams it tried. Searches for a rare body (a giant, a world that makes oxygen, a world that draws a trace gas) run over a wide range so a generator change does not empty them; the scored sample of `tests/scored_site.py` now picks its four systems by the moons they hold instead of three fixed seeds.
- **The oxygen-breather schedule test is built by hand.** Generation no longer makes a class 4a world short of oxygen, so `test_a_class_that_breathes_oxygen_but_makes_none_stops_at_photosynthesis_in_an_air_without_it` dates a class 4a on a generated world with a thin air (same assertion).
- **The belt page test counts the planets that were bodies of the belt** (their designations are in `<code>` too).

generation-output: unchanged
generation-affects: none
