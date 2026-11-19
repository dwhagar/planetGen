**Headline:** Two of the slowest gate files skip the first fill's draw of every sector they make, which their assertions never read.

### Fixed

- `tests/test_gen1254_scatter_uid.py` (50 minutes in the Train 49 gate) and `tests/test_deferred_sectors.py` (55 minutes) made whole sectors and so ran each sector's first fill down to its deepest floor, thousands of systems, though they check the IDs of the scatter stars' systems and which sectors a plan cuts and saves. Both now use `first_fill.without_the_deepest_floor`, as `test_gen_resume.py` and eight other files already do (TEST.417). The assertions are unchanged and every test still passes: the scatter file runs in 2 minutes 18 seconds (was 50 minutes), the deferred file in 22 minutes under a gate's load (was 55 minutes).

generation-output: unchanged
