**Headline:** A neighbourhood fill now generates the centre sector first and fills outward in a fixed, space-filling order.

### Changed
- **Fill order for neighbourhoods and radial fills (GEN.101).** `planetgen galaxy --center-sector ID --radius-pc R`, a random start's neighbourhood and a single address with a radius used to generate their sectors in the order the grid enumerated them, shell by shell, so the first tenth of a run held almost none of the nearest tenth. They now start with the sector at the centre and go outward in blocks of four sectors a side, the blocks nearest first and ordered along a Hilbert curve inside and between equal distances (`generation/fill_order.py`), so the order is deterministic whatever order the candidates come in. A ball cannot always be walked between neighbouring sectors with no jump (Boss kept the Hilbert order with a logged jump), so the run prints one line with the share of its steps that are between neighbours (about 97% in a 50 pc ball), the number of jumps and their mean and longest length in sector edges.

TODO items resolved: GEN.101. Added: none.
