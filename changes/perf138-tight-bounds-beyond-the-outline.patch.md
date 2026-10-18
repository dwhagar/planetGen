**Headline:** The shells beyond the galaxy's outline generate in seconds instead of hours, because the polar reach groups no longer take the dense core's nuclear density as their ceiling.

### Changed
- **A reach group's nuclear bound is read where the group really is** (`generation/object_first.py`, `bright_stars.py`, `phenomenon_reach.py`). A group beyond the outline reaches all the way to the pole, so its smallest cylindrical radius is near zero and the nuclear cluster and nuclear disc were bounded by their central density, thousands of times the density anywhere in the group. The scatter then drew millions of candidates to keep a few. The bound now uses the group's inner spherical radius as well as its cylindrical one, and cuts the height range of the nuclear disc into 12 steps, so the polar groups' bounds sit within 3 times the largest density a group holds. A shell-5 call that took over an hour (still unfinished after 60 minutes) on main now takes 43.5 s for 2,382,680 rows.
- A looser or tighter bound changes how many candidates are drawn, so the rows a seed places in these shells change; the generator epoch is 35.
generation-output: changed

### Added
- `tests/test_gen_reach_groups.py`: a polar group's bound stays within 3 times its largest true value.

TODO items resolved: PERF.138. Added: none.
