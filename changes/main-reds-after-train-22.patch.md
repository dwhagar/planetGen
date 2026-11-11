**Headline:** Three groups of tests that went red on main after Train 22 now look for the case they need instead of taking the first one.

### Changed

- TEST.714: `test_web_generate.py::test_redo_scatters_clears_then_draws_each_ticked_scatter_with_its_setting` claimed a scatter step's own stage comes first; since OPS.123 slice 3 (#5876 for the bright stars, #5896 for the phenomena) a scatter lists what it empties first. It now says the step's own stage is in the list, only clears come before it, and the rest follows. Nothing in the code changed.
- TEST.714: `test_effect_ledger.py` took the first binary system with planets, but a lone source is implied away and writes no row, and which binaries have a planet that two stars light is the draw's (seed 0 after Train 22 had none). The helper takes the first that has one.
- TEST.714: `test_api_life_air.py` took the first airy rocky planet; the history has a methane step only when the air can hold methane and still show it, which that planet's draw decides (its two tests failed with no CH4). The helper takes the first of the same draws whose history through the API has methane. No claim of either test is weaker.

generation-affects: none
