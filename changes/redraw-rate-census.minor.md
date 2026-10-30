**Headline:** A test now counts how often the generator throws a drawn star tree, planet, moon or belt away to make it fit, and fails when a rate doubles.

### Added
- **Redraw counters (TEST.476).** `draw.redrawn(kind)` and `draw.counting()` in `planetgen/util/draw.py` count, without reading a stream, every place a system build redraws, moves or cuts an item to make it fit: the pair orbits of a star tree (Kozai, stability, gravitational waves, interaction gates, overlap, no allowed period), the tree restarts and the dropped member, the whole companion set drawn again (the GEN.1004 clearance), the planet list drawn again for a forced habitable world or belt, planets and belts pushed out by the Hill and belt spacing, the circumbinary floor, the orbit-ceiling and cross-star trims, moons that would be too large, comets, giant heavy-element masses, snow-line scatter and the primary star's model draws. Thirty-five kinds are registered in `planetgen/generation/redraw_kinds.py`, and a test fails when a call names a kind that is not registered or a kind has no call.
- **The redraw-rate census** (`tests/audit/redraw_rates.py`, `tests/test_redraw_rates.py`). Seeded systems of nine groups (ordinary, one to six stars forced, a habitable world and belt requested, habitable worlds barred) give the events per generated system for every kind, by group, by the number of stars the system came out with and by the class of the placed star, with a 95% interval of the share of systems that redraw at least once. The report is `redraw-rate-latest.md` beside the other census reports. `python -m tests.audit.redraw_rates measure --count 1500 --groups stars-3,stars-4,stars-5,stars-6` draws the big sample the baseline was measured on.
- **Two gates.** A ceiling per kind and group in `tests/data/redraw_rate_expectations.json` (registered as stand-ins of tier T9, tweak L3: twice the baseline measured before the overlap-free build plus one event of the reference's size) fails when a rate is above it by more than the run's own sampling error. The overlap kinds (`pair.overlap`, `tree.restart_overlap`, `companions.overlap`) are held at exactly zero by an invariant that is an expected failure until GEN.1015 deletes the GEN.1004 quick fix.

No generator output changes: a counted system is the same system.

generation-output: unchanged

TODO items resolved: TEST.476.
