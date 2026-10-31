**Headline:** Stars now have a proper integrator for their path through the galaxy, ready for the orbital update to use.

### Changed
- **The star stepper (GEN.340).** `planetgen/galaxy/stepper.py` moves stars through the galactic potential with second-order kick-drift-kick steps, on arrays so a whole sector's stars are one call. Each star's step is a twentieth of its own crossing time (no more than 5 Myr), cut into equal power-of-two blocks so stars of one level step together; a star that would need more than 4096 steps in a span is stepped at that cap and flagged. A callable adds an extra pull to the same kick (the influence set of GEN.109 goes in this way). A star at the Sun's kind of radius keeps its radius to 1 pc and its angular momentum to 1e-10 over 5 Gyr, and a century in one step agrees with a hundred one-year steps. Nothing calls it yet: the orbital update run (GEN.105) will replace the turn about the galactic pole with it, which changes where stars are, so that rides a train; generation output does not change now. Stand-in rows ST1 and ST2.

generation-output: unchanged
