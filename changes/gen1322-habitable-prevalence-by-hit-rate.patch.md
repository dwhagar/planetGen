**Headline:** Asking for more habitable worlds in a run now forces enough systems that the share with one comes to the wanted figure, up to the 74.7% a forced system can give, and the Generate page says a habitable world is only asked for.

### Changed

- `src/planetgen/generation/prevalence.py` (GEN.1322): `resolve` forced a system on in `p * base / (1 - base)` of cases, which assumes every forced system has the feature. A habitable world is asked for, not placed (GEN.1130: the conditions' odds are multiplied by ten and the request never forces), and the count is of planets and moons of a class that can carry life, so 224 of 300 forced systems had one (0.747). The chance is now `p * base / (hit - base)`, with `hit = tuning.PREVALENCE_FORCED_HIT_RATES` (1 for the belt, large star and planets, so those are as before), and every system is forced from the share a forced system can give (`reachable_share`). +50 had read 0.62 to 0.66 against 0.876 asked; it now forces every system and should read about 0.75, the most that can be had; +20 asks for 0.70 and forces 80% of systems.
- `src/planetgen/tuning.py`: `PREVALENCE_FORCED_HIT_RATES` `{"habitable_world": 0.747}`, measured over 300 systems with `HABITABLE_WORLD` True (run seeds 801 to 803), remeasured with the base shares.
- `src/planetgen/web/generate_page.py`, `partials/generate_fields.html`: the habitable-world prevalence field says "A higher share raises the odds and never forces one; the most it gives is 74.7%."
- `src/tests/test_prevalence.py`: the measured-feature test expects at most the reachable share; tests of the forced chance at +20, +50 and -100 and of the page's limit.

generation-output: unchanged
