**Headline:** A very hot outgassed world no longer crashes its whole system: atomic chlorine, atomic fluorine and fluorine are now in the gas table (GEN.1337).

### Changed

- At a surface near 2500 K the hydrogen halides dissociate and a few thousandths of the outgassed gas is Cl, F and F2. The gas equilibrium already carried their free energies, but the gas table did not list them, so `mixture_base.finish` refused the draw (a MixtureError, "0.00143 of the gas is species the gas table lacks") and the whole StarSystem failed.
- The three rows (Cp from the NASA Glenn fits, as for the other atoms) are in part 01 section 5 of the atmosphere research, in `scripts/build_gas_table.py` and in `gases.json`, now 115 rows. The error itself stays: a gas outside the table still stops the draw.
- A hot surface's gas now lists Cl, F and F2 where it did not, so those mixtures change (epoch).
