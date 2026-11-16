**Headline:** The tables a sector's generation writes are declared in full, so the touches check no longer depends on which bodies a test galaxy happens to generate.

### Changed
- OPS.241: `queue/touches.py` `SECTOR_WRITES` (what `generate-unplaced-sector` and the other commands that save a sector write) lacked `black_holes`, `neutron_stars`, `planet_reflection_spectrum`, `moon_reflection_spectrum` and `moon_interior_layers`. `test_command_touches` failed on Train 37 when its generated galaxy held a neutron star or a body with a reflection spectrum, and passed on main only because the galaxy there did not. No generation output changes.
