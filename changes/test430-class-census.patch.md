**Headline:** A standing class census draws planets and moons per star type and zone and flags a class share, a habitability colour or a life ceiling that leaves its stored range.

### Added
- **The class census (TEST.430, with TEST.315's class shares).** `src/tests/audit/prevalence.py` draws seeded planets and moons around M2V, K5V, G2V, F5V, A0V and giant stars in the hot, ecosphere and cold zones, and `src/tests/data/prevalence_expectations.json` holds a min and max for every class share, every Blue, Green, Yellow and Red share and every life ceiling share, each marked observed, researched norm or stand-in (here the prevalence before 4c, commit `0b4d0cc`, 200 draws per cell). It flags a share out of range, a class missing where expected, a class present where barred, and a draw over the time limit; a known fault carries its item (GEN.855 for the oxygen loss since 4c, GEN.829 for class T never drawn, GEN.866 for the interior-heat hang) and the report lists it apart. `python -m tests.audit.prevalence_baseline` re-baselines the file. Documented in `docs/reference/testing.md`.

TODO items resolved: none (TEST.430 stays open for its later parts). Added: none.
