**Headline:** No test skips quietly any more when its seed did not draw the body it needs.

### Changed
- **Seed-dependent skips are gone.** About 60 tests ended in `pytest.skip("no ... in these draws")` when an unseeded
  search found nothing, so a run could pass with the check never made. Each search now draws from fixed streams
  (`draw.bound(...)`) that are known to produce the body, and a search that still finds nothing fails with a message
  that names how many fixed streams it tried.

generation-output: unchanged
generation-affects: none
