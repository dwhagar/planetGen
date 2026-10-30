**Headline:** The "what is near here" search works out each candidate's distance once instead of twice with the subtraction repeated, so a 50 pc count is about 2.4 times faster.

### Changed

- PERF.216: `planetgen.db.region.region_statement` computes the squared distance as `POW(axis - c, 2)` once per candidate and tests it with `HAVING d2 <= r²`; `capped_count` (the early stop count of NAV.61) no longer computes a distance it does not use (`with_distance=False`). On a 3 million row table the count of a 50 pc sphere (62,000 systems) takes 178 ms instead of 428, the same sphere with the rows returned 425 instead of 600, and 25 pc 32 instead of 65. The systems returned and the counts are the same as before, a system exactly on the radius included; the plan is the same. The design note `docs/design/galaxy/spherical-region-queries.md` section 8.2 has the profile. No stored value changes.
