# planetgen/datasets/catalogue.py

"""
The datasets the setup step downloads (GEN.548): one `Dataset` per source, each with the exact size and sha256 of
every file, the citation, the licence and the stand-in the readers use until the files are present.

The three of GEN.548 (isochrone grids, white-dwarf layers, the neutron-star equation of state) and the data of the
atmosphere group (AD1 to AD8) are added here, each by the item that reads it; every entry needs its verified source
(docs/design/gap-filling-plan.md 6.5) and a row in the stand-in registry (section 6).
"""

ALL = ()
"""tuple: Every `fetch.Dataset` the setup step downloads, in the order it downloads them."""
