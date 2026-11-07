**Headline:** A new setting, `generation.cell_cap` in config.json, sets how many systems a crowded cell may hold before it is only partly filled (5,000 unless changed).

### Changed

- `generation.cell_cap` (GEN.303, Boss ruled 2:21 PM EST on 2026-11-07): the deployment's cell star ceiling, from 500 to 50,000, default 5,000. A plan reads it when the run starts (`run_plan.skeleton_params`) and stores it with the galaxy as that galaxy's own `cell_star_ceiling`; `--cell-star-ceiling` still sets one galaxy's value and wins over it. At the default a galaxy generates exactly as before, with no migration and no epoch. A change applies only to galaxies planned afterwards, never to cells or galaxies already stored. It is on the settings page under "Galaxy generation" and in `docs/reference/config.md`, `config.json.example` and `config.schema.json`; the `--cell-star-ceiling` help names it.
- `src/tests/test_gen_cell_cap.py` (new): the bounds and default, the plan taking the setting and the command line winning, a lower cap leaving fewer systems resolved in a crowded cell, and the option in the example file, the docs table and the help.

generation-output: unchanged
