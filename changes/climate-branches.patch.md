**Headline:** The ice-albedo feedback now has a tested module that finds the snowball, warm and unstable climates of one flux.

### Added
- **GEN.1046: the branches of a climate with an ice-albedo feedback.** `physics/climate_branches.py` holds the albedo law of the ice (ice-free 0.30 above 283 K, covered 0.62 below 263 K) and finds every fixed point of `T = response(albedo(T))`, with its stability and the flux window in which two climates exist. Design in `docs/design/composition/climate-branches.md`; stand-in rows CC15 and CC16. Nothing in the generator reads it yet (GEN.972 and GEN.1051 wire it), so no world changes.
- generation-output: unchanged
