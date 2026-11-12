**Headline:** A sub-Neptune that lost half of its envelope or more to its star now ends with the helium that hydrogen-first escape leaves behind, and the planets that reach a helium share of 0.4 are named as having a helium-rich envelope (oddity W16).

### Changed

- GEN.1527: `helium_fractionation` is read by the generator. At the capture step the loss history of the envelope (`envelope_escape.envelope_left_earth` with a new `trace`) is run through it for every planet that lost at least half of its starting envelope (`envelope_capture.HELIUM_LOSS_GATE`, registry HF4); the helium share it gives replaces the star's in the stored `envelope_helium_fraction` and in the helium budget of the envelope mixture (`mixture_giant.giant_mixture`, new subject quantity `envelope_helium_y`), so the air of such a planet is helium-richer than its star's gas. Oddity W16 (lever L8, no draw of its own) names the planets at 0.4 or over. A planet that lost less keeps the star's gas as before.
- Planets with a deeply stripped envelope draw a different air, so seeds change for them: the golden and the epoch are Lane 02's to regenerate.

generation-output: changed
generation-affects: fill
