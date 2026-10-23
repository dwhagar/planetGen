**Headline:** The saturation-pressure function now covers liquids as well as ices, up to the critical point.

### Changed
- **The liquid branch of the saturation-pressure function (GEN.452).** `physics/saturation.py` now answers above the triple point too, with the Wagner fits to the NIST saturation tables for hydrogen, neon, argon, krypton, xenon, nitrogen, oxygen, water, carbon monoxide, carbon dioxide, methane, ethylene, ethane, methanol, ammonia, hydrogen sulphide, sulphur dioxide and benzene and a correlation for hydrogen cyanide; hydrogen's solid comes from its triple point. Each answer says whether the condensate is solid or liquid, and above the critical temperature it says there is no saturation. At the triple point the two branches meet within a few percent (benzene 7), which is documented, not smoothed. Nothing reads the function yet, so no generated value changes.

generation-output: unchanged
