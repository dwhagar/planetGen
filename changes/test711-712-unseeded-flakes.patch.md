**Headline:** Two tests that failed about half the time on any machine, the deferred-climate moons test and the wide-binary re-roll console test, no longer depend on luck.

### Changed

- **TEST.711 and TEST.712 were coin tosses, not regressions.** Neither is a change of main: the deferred-climate test asked `moon_count=2` and expected exactly two moons, but a planet built the usual way with the same request also has one moon about a third of the time (a moon that cannot stay is not made), and the run stream is seeded from the system's entropy, so the test failed on about 40 percent of runs (11 of 30 seeds at the commit that added it, 19 of 40 on main). It now draws from fixed seeds until the planet is given its two moons and checks the state before and after `finish_climate` for each. The wide-binary re-roll test expected the console line to contain "orbits beyond the farthest stable orbit", but the line names only the first three problems (`REPORT_PROBLEMS`), and a belt that stays at its place is a problem that comes first; it now expects the first three problems the scan itself found, and that the system has a "beyond the farthest stable orbit" problem among them. No code changes.

generation-output: unchanged
