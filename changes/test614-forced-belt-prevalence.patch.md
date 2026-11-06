**Headline:** The test of the forcing rules asks for a belt in every system with a prevalence large enough to do it, and passes again.

### Changed

- `test_prevalence.py::test_resolve_keeps_forced_flags_and_the_forcing_rules` (TEST.614) asked for a belt at +1000 percent and expected one in all 300 draws. A positive prevalence forces a feature in `p/100 * share / (1 - share)` of the systems, "all of them" only from 100 percent of that; the belt's usual share is 0.039 (`tuning.PREVALENCE_BASE_SHARES`), so +1000 percent forces it in about 41 percent of systems and the test failed on main whenever a draw missed. The rule was right and the test's number was stale. It asks +10,000 percent now, and the second half of the test (a belt forced on that the small star undoes) asks the same, so the undo is what leaves the belt unset in every run, where at +1000 it was exercised in about four runs in ten.

generation-output: unchanged
