**Headline:** The atmosphere stand-in registry names the four tests whose names DOC.103 respelled, and its test passes again.

### Changed

- `docs/design/research/gap-filling-plan.md` rows AM10, AM15, the impactor-escape row and the alumina row named tests by their British spelling (`vapour`, `oxidised`); DOC.103 respelled the tests (`vapor`, `oxidized`) and left the rows, so `test_atmosphere_standin_registry.py::test_a_row_names_a_constant_and_a_test_that_exist` failed in four cases on main since train 10. The rows now name the tests as they are. `stand-in-values.md` and `docs/references/` are rebuilt from the plan.

generation-output: unchanged
