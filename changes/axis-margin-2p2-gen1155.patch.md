**Headline:** The 2+2 quadruple-star stability margin is now fitted on our own N-body runs (1.5, was 1.3 carried over from the 3+1 case).

### Changed
- **2+2 axis margin (GEN.1155, GM14).** 302 new REBOUND 2+2 systems (masses log-uniform over a decade, axis ratios 0.02 to 0.5, importance-sampled around the network's boundary) and the 972 of GEN.1031 show 1.3 lets too many unstable systems through (precision of the stable calls 0.81 on the paper-like sample against the paper's 0.94). 1.5 reaches 0.94 there and 0.99 on the wider sample. `quadruple_stability_2p2.json` `axis_margin` goes from 1.3 (T6, L4) to 1.5 (T6, L2); the GM14 row, the burn-down line (now `replaced`) and `quadruple-stability-2p2.md` section 5 carry the numbers. The sampler is `nbody2p2_margin.py`. Nothing calls the 2+2 verdict from the generator yet, so no stored number changes.

generation-output: unchanged
