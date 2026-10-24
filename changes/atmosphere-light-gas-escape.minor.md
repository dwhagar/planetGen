**Headline:** A small warm world no longer keeps hydrogen, helium or other light gases in its drawn air.

### Added
- **Thermal escape filter for the atmosphere mixture (GEN.693, `physics/mixture_escape.py`).** After the chemistry of a rocky, Venus-like, Mars-like or abiotic-oxygen body, each gas is kept or lost by its Jeans parameter at the body's exobase against the escape speed and age, so a 3.5 km/s world at 310 K ends airless instead of holding a hydrogen and carbon monoxide air, and Earth keeps a few ppm of helium. The exobase temperature is a labelled stand-in. Not called by the generator yet: no generator epoch.
- **An oxygen cap (`physics/mixture_oxygen.py`).** A lifeless world that lost its water is drawn with at most 95 percent oxygen and at least 0.5 percent nitrogen plus argon; `BIOTIC_O2_MAX` (35 percent) is the cap the life stage applies to a living world. No air is pure oxygen.

### Changed
- `mixture_draw.required_inputs` for the rocky classes now also names `v_esc_kms`, `teq_k` and `xuv_index`.
