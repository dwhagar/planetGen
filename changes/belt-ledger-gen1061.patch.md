**Headline:** A belt's biggest bodies can now be drawn one by one from the belt's mass and size law, and the mass left over is kept, so nothing is lost or counted twice (not used by the generator yet).

### Changed
- **Belt ledger (GEN.1061).** New `physics/belt_ledger.py`: the size laws of an icy and a rocky belt, the heaviest-body order statistic, the mass cut rule (every body above 1e19 kg, up to 90 percent of the belt's mass, at most 300), mass conserved, ranks stable as the belt ages, and the largest body settled with the belt's collisional mass. The belt's mass is now the primordial one when the largest body is too strong for collisions to break (it was a negative number). The size-law note had Ivezic's main belt indices as cumulative; they are differential (cumulative 1.3 and 3.0). Nothing calls the ledger yet.

generation-output: unchanged
