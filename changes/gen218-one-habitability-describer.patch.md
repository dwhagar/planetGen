**Headline:** The habitability wording now comes from one module: the second explainer, which nothing used, is gone, and the numbers in the sentences use the unit ladders' own formatter.

### Changed
- `physics/protection.py` is the one describer of the six PHI-6 categories (rows of domain, colour, score and sentence, the paragraph, and the gear). `physics/habitability_explain.py` (its `factors` and helpers were used by one test only, its `paragraphs` was a one-line wrapper) is deleted; `Planet` calls `protection.paragraph` directly and the useful tests moved to `tests/test_protection.py`.
- The sentences print numbers through `util/format.three_figures_or_tiny` (the formatter the ladders use; it was private) and name gases from `atmosphere.SPECIES_NAMES`, so there is one number formatter and one gas-name table. A very small dose now reads "2.00 × 10⁻⁷ Sv a year" instead of "2.0e-07".
- The docstrings of `habitability.py` and `physics/__init__.py`, `docs/html-interface.md` and `docs/design/architecture.md` say PHI-6 and no longer name the deleted module. The stored names (`phi4`, the `"PHI-4"` key) stay.

TODO items resolved: GEN.218 / Added: `docs/design/architecture.html` still lists `habitability_explain.py` (generated page, DOC.67 owns it)
