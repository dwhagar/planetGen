**Headline:** A 3+1 quadruple star system can now be tested for stability with a published classifier that is right about nine times in ten.

### Changed
- New `physics/quadruple_stability.py` (GEN.1029): the 3+1 quadruple stability network of Vynatheya, Mardling and Hamers 2023 as plain numbers and a numpy forward pass (no machine-learning library). On my N-body check of 331 systems it agreed with the integration 95 percent of the time against 49 percent for the nested-triple test. Nothing in generation calls it yet (GEN.1016 wires it), so no generated world changes.
- Registry rows QS1 to QS4 (section 6.65), the design note `docs/design/systems/quadruple-stability-3p1.md`, and `scripts/build_quadruple_stability.py`, which reads the authors' model file without running its code.
