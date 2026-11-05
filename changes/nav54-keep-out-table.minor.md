**Headline:** The keep-out table is complete: asteroid fields, the central black hole, quasars, and stars inside their own system now have a keep-out radius.

### Changed
- An asteroid field keeps out its stored radius. A supermassive black hole keeps out its sphere of influence (3.1 parsecs for four million solar masses) instead of its event horizon, and a quasar keeps out the larger of that and the radius where its light is 10 kW per square meter.
- A star carries a second keep-out for a course inside its own system, `keep_out_in_system`: the larger of ten of its radii and the distance where its light is 50 kW per square meter (0.17 AU for the Sun). The Galaxy keep-out of a star is unchanged.
- The note on a nebula or a supernova remnant now says why a course can pass through it (thin gas), and an asteroid belt carries a note too. The knobs are in `tuning.py` (`KEEP_OUT_*`).

generation-output: unchanged
