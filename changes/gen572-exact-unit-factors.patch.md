**Headline:** A parsec is exactly 3.2615637772 light-years and a light-year exactly 63241.0770843 AU wherever the generator used rounded values.

### Changed
- **The three rounded unit factors are the exact ones (GEN.572).** `physics/hii.py` typed 3.26156 light-years per parsec, and `generation/phenomena/asteroid_field.py` and `generation/phenomena_plausibility.py` each typed 63241.077 AU per light-year; all three now read the definitions in `physics/constants.py` (`PARSEC_M / LIGHTYEAR_M` and `LIGHTYEAR_M / AU_M`), which differ from the old values by 1.2e-6 and 1.3e-9 in relative size. `hii.region` returns its radius in light-years with the factor, and the nebula generator stores it, so a nebula's radius changes in the last digits for a seed; an asteroid field's size digit can differ only for a radius at its floor boundary. Nothing else changes.
- **The generator epoch is 25** (a pending epoch until a release names it); the full wiped rescatter that carries the epoch's other changes carries this one.
generation-output: changed

TODO items resolved: GEN.572. Added: none.
