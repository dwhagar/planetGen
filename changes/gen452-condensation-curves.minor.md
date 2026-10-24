**Headline:** Where an atmosphere freezes out is now worked out from the checked saturation curves, so hydrogen sulphide, carbon dioxide and sulphur dioxide land at the right temperatures.

### Changed
- **The freeze-out cap reads the checked saturation curves (GEN.452).** `physics/atmosphere.py` no longer carries the two-point condensation table, whose hydrogen sulphide row was a copy of sulphur dioxide's and which put carbon dioxide, sulphur dioxide and the cold gases 2 to 20 K off. The most a gas can be vapour at the surface temperature is now its saturation pressure from `physics/saturation.py` (the Fray and Schmitt and laboratory solid curves, the Span and Wagner carbon dioxide curve, the Wagner liquid curves, nothing above the critical temperature). Hydrogen sulphide boils at 213 K at 100 kPa instead of 263 K, nitrogen at 37 K holds 1.2 Pa (Pluto's air) instead of 10 Pa, and the cold gases follow the laboratory curves of the cold ends. New `saturation.condensation_temperature_k` gives the inverse. The cloud condensation of `physics/aerosols.py` reads the same curves. Stored partial pressures and total surface pressures change for cold bodies.

generation-output: changed
