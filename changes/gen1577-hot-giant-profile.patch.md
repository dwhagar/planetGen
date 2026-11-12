**Headline:** A hot giant's atmosphere can now be given a temperature profile from its light, its own heat and its gravity, ready for the cloud deck to be placed on.

### Changed

- GEN.1577 (built for Lane 14): new `physics/hot_giant_profile.py`, the Guillot 2010 semi-grey irradiated profile (`temperature_k`, `profile_function`, `skin_k`) over 1e-6 to 1e3 bar, checked against the research's worked example (1473 K at 66.7 mbar, 1647 K at 1 bar, 1648 K at 10 bar); stand-ins PF1 to PF3 (registry 6.277) and the design note `docs/design/composition/hot-giant-profile.md`. Nothing calls it yet; Lane 14 wires it behind the cloud deck.

generation-output: unchanged
generation-affects: none
