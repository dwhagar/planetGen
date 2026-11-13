**Headline:** The independent overlap verifier stays outside the generator on purpose; the wiring register now says so in one row instead of two.

### Changed
- **`physics/clearance.py` is by design read only by the tests (GEN.1685).** It must not share a bug with the generator it checks, so the stability gate of GEN.1016 does not call it; the seed sweep in `tests/test_clearance.py` builds, stores and audits systems with it. The register's second row (waiting on the gate) is gone and the first says what reads the verifier.

generation-output: unchanged
generation-affects: none
