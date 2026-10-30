**Headline:** The property sheet audit no longer excuses the moons' induced magnetic moment as an open gap; the column has been filled since the ocean-moon induction was built.

### Fixed

- **The induced moment of a moon is no longer listed as a known gap (GEN.779).** The fill audit excused the moons' `induced_moment_a_m2` as an open gap, written when the column was believed empty. A moon's induced dipole is set in the final pass of every system (`magnetism.generate_induced_field`, GEN.178: the moon's ocean answering its parent's turning field, 0 for a moon with no liquid ocean, no parent field or a parent that turns with it); a sample of 30 systems has it filled on every moon and non-zero for the ice-covered-ocean moons. The excuse is removed, so the audit now fails if the column is ever empty again. A planet keeps its sheet rule (it orbits its star alone, and its field answer is the part of GEN.779 that rides a train). No stored value or generated output changes.
