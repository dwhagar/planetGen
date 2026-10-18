# planetgen/datasets/catalogue.py

"""
The datasets the setup step downloads (GEN.548): one `Dataset` per source, each with the exact size and sha256 of
every file, the citation, the licence and the stand-in the readers use until the files are present. Each module named
here holds its dataset and the functions that read it: the MIST main sequence (`isochrones`), the white-dwarf cooling
sequences with their layers (`white_dwarf_models`) and the neutron-star equations of state with the outer-crust
compositions (`neutron_star_matter`). The data of the atmosphere group (AD1 to AD8) is added the same way, each by the
item that reads it; every entry needs its verified source (docs/design/gap-filling-plan.md 6.5) and a stand-in row
(section 6.6).
"""

from planetgen.datasets import isochrones, neutron_star_matter, white_dwarf_models

ALL = (isochrones.DATASET, white_dwarf_models.DATASET, neutron_star_matter.DATASET)
"""tuple: Every `fetch.Dataset` the setup step downloads, in the order it downloads them."""
