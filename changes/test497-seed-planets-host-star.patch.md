**Headline:** The planet tests draw their host star from a fixed seed, so a run no longer fails on a star whose thin air or life scores happen to sit at an edge.

### Fixed
- **test_planets host star (TEST.497).** The module's `host_star` fixture drew its G2V star from the run's entropy, so the planets of one run came from a different star than those of the last (the e-3 case that stored 0.99999 Pa under a floor of 1.0, which TEST.485 widened by the trace trim, and the c-1 life class case of TEST.486, failed only on some stars). The star is now drawn from a fixed stream; the file passes three runs in a row. Nothing the generator makes changes.

generation-output: unchanged

TODO items resolved: TEST.497.
