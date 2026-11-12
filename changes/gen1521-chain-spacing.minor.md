**Headline:** Planets in a tight chain of four or more are now spaced further apart, so the chains an N-body check showed scattering no longer come out of the generator.

### Changed

- A planet that would be the fourth or later in a chain of planets is held to an effective spacing of 12 mutual Hill radii over the gaps of its window of four (`tuning.CHAIN_K_EFF`, `clearance_ledger.chain_distance_au`, GEN.1521): the mean of the three gaps less half their spread, the combination that sets a system's survival in Pu and Wu 2015, at their minimum for a gigayear of stability. The pairwise floor of 10 stays. Lane 04's N-body check (GEN.1468) found three chains of four to twelve small planets, every pair at 10 to 11.4, that scattered within 20,000 years; all three have an effective spacing of 10.1 to 10.5 and are refused.
- On the 145 systems of the check with four or more planets, 108 (74 percent) have a planet moved outward, by 21 percent at the median; at an effective spacing of 10.5 it would be 61 systems (42 percent) and 3.7 percent. Systems of three planets or fewer do not change. The setting is one number in `tuning.py`; 0 turns the test off.

generation-output: changed
generation-affects: plan
