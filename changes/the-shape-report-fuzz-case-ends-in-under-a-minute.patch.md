**Headline:** The shape-report route's fuzz case, which ran past 900 s, now ends in under a minute.

### Fixed
- **`test_admin_pages_survive_hostile_query_when_logged_in[GET /admin/generate/shape-report]` ran past the 900 s test limit (TEST.275).** The fuzz sends many queries, and every one of them, junk keys or not, planned the default outline (`/admin/generate/shape-report` answers the line under the dials), which took 9 s to 10 s: the mean-system-mass age cache held 512 ages where one plan asks for about 23,000 (fixed by PERF.163), and the cell count of every shell was added up wedge by wedge again for each of the 1,532 shells, 1.8 million wedge calculations each time. `spherical.cells_in_band_range` is now cached like `shell_layout`, so the repeats of the same plan cost the rest of the plan only (about 2.4 s). The case alone takes 47 s, down from over 1,780 s.

TODO items resolved: TEST.275. Added: none.
