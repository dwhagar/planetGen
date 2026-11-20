**Headline:** The request-dial test checks the request draw alone, not the sparse tier's other dials.

### Fixed
- **TEST.806: the test that a tier with no request share makes the same system as science failed on main.** Sparse's water, dry-part, quiet-interior and carbon-cycle dials move the planets by design (GEN.1720, GEN.1721, GEN.1798 to GEN.1801), so seed 103 gave 18 planets on science and 16 on sparse. The test now holds the sparse row at no dial, which is what it means to check: a request share of 0 draws nothing.
- generation-affects: none
