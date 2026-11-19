**Headline:** The dense-cluster stripping test tries 120 seeds instead of 40 to find a system that loses a planet to its birth cluster.

### Changed
- **`test_cluster_stripping` tries 120 seeds, not 40.** A seed qualifies when its plain system ejects nothing, its dense system strips a planet and holds fewer planets; that is about one seed in 40 (seed 7 on main; the growth limit by star mass moved it to seed 68). The three conditions and the assertions are unchanged.
