**Headline:** Six constants that the constants-completeness test named now carry their tag, and one of them is no longer a second copy of another's value.

### Fixed
- **`test_constants_completeness` was red on main** for six constants with no registry row or tag. `compaction.DRAWN_POROSITY_BELOW_KG` (GEN.1637's compaction) is tagged as derived from `belt_ledger.L1_MASS_KG` (it cannot import it without a cycle), and a new test holds the two equal. `disk_dispersal.YEAR_PER_MYR` (GEN.244) is tagged a unit conversion by definition. `gas_equilibrium.SOLVE_MEMORY` and `thermosphere.EXOBASE_MEMORY` (PERF.276, #6993) are tagged cache sizes that change no result. `interior_heat.FLOOR_STIFFNESS_MARGIN` (PERF.277) is tagged a numerical safety factor derived from the integrator's stability limit. `tuning.SLOW_PAGE_WANTED_SECONDS` (OPS.230) is tagged as derived from the background-tab timer. No physical value changes and no row is exempted.

generation-output: unchanged
generation-affects: none
