**Headline:** The scatter's share of recycled neutron-star pulsars follows the new age-based pulsar rule, which fixes a broken import left when two merges crossed.

### Fixed
- `NEUTRON_STAR_RECYCLED_SHARE` (the share of scattered neutron stars that are recycled millisecond pulsars, which decides their companions) was built from `NEUTRON_STAR_PULSAR_CHANCE` and `PULSAR_MILLISECOND_CHANCE`, which GEN.153 removed, so importing `generation/phenomenon_scatter.py` failed. It is now `PULSAR_RECYCLED_CHANCE * PULSAR_RATE_SCALE`, about 1e-4, the real fraction of old neutron stars that are recycled pulsars (was 21%).

TODO items resolved: none / Added: none
