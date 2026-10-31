**Headline:** The property sheet now says why a heavy body shows no porosity and why a shallow ice shell shows no high-pressure ice.

### Changed
- **GEN.780.** `porosity` has been stored for every non-giant body lighter than 1e21 kg and for none above it (checked on 60 generated systems: no exceptions), because above that mass gravity has closed the voids. The sheet now carries a rule with its sentence of reason (the compaction taper end), and the open gap that excused the column is removed.
- **GEN.781.** `hp_ice_km` has been stored exactly when a body has an ice shell state and the water's base pressure reaches 620 MPa, where ice VI begins (60 generated systems: 47 with, 120 without, 0 exceptions). The sheet now carries that rule with its reason, and the open gap is removed.

Nothing a generator draws or stores changes.
generation-output: unchanged
