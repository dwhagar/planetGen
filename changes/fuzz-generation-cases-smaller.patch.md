**Headline:** The system-generation fuzz tests are split or sampled so none passes the 900 second limit with every flag forced on; no generated output changes.

### Changed
- **Smaller fuzz cases (TEST.417, PERF.241).** A system with every orbit flag forced on costs 15 to 96 s for an O-type star (profile in reports/perf241-all-flags-o0v-profile-2026-11-05.md), so ten of them in one case passed 900 s. `test_every_star_type_generates_a_sane_system` now runs each letter and class in two cases, subclasses 0 to 4 and 5 to 9, with the same systems and seeds (all 70 combinations are still built); the slowest case of the M class took 215 s. `test_extreme_orbit_counts` (30 to 12 examples), `test_every_forced_flag_on_together` (25 to 10) and `test_random_valid_configs_generate_sane_systems` (20 to 10) draw fewer examples in the ordinary profile; the deep profile (`PLANETGEN_FUZZ_PROFILE=deep`) draws thousands as before. Restore the counts when GEN.1282 and PERF.242 cut the cost of a system.

generation-output: unchanged

TODO items resolved: none (PERF.241 continues).
