**Headline:** A rocky body's air history asks for its interior heat at a hundred ages and each is no longer run from the formation: the box takes up the run of an earlier age, 4 times faster, with the same results to the last bit.

### Fixed
- PERF.289: building a system with planets took 12 to 38 s where it took 3 to 12 s before the air history (GEN.1745), and in a profile a quarter of the history's time was the mantle box (`interior_heat.heat`) run again from the formation for each of the history's ages (about 220 runs for each body with air). The steps of the box do not depend on the age asked for except the one that lands on it, so a run now keeps the steps and the states of the box (every 16th step, the core copied) and an older or younger age takes it up at the last step it shares and goes on from there. Same output: 8,100 results (60 bodies with and without a tide, three orders of ages) are equal to the last bit to the runs from the formation, and the new test compares the two in whatever order the ages come. Heat for those 60 bodies went from 15.1 s to 3.9 s.

generation-output: unchanged
