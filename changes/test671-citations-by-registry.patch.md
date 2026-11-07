**Headline:** Five sources named in code and docs (Kasting 1988, Simpson's limit, Van Eylen 2019, Hadden and Lithwick 2014, Fulton 2017) are written by their registry IDs, so the citation guard passes outside the TODO text.

### Changed

- **Sources written by registry ID (TEST.671).** `radiative_convective.py`, `scripts/validation_sweep_planets.py` and `docs/design/research/gap-filling-plan.md` no longer name these papers by author and year alone: Kasting 1988 is [R2275], Van Eylen and others 2019 is [R1307], Fulton and others 2017 is [R1478], Hadden and Lithwick 2014 is the new [R5137], and Simpson's limit is cited through Goldblatt and Watson 2012, the second-hand source, as the new [R1608]. Only comments, docstrings, a report line and the registry change.

generation-output: unchanged
