# planetgen/datasets/catalogue.py

"""
The datasets the setup step downloads (GEN.548): one `Dataset` per source, each with the exact size and sha256 of
every file, the citation and the stand-in it replaced (named in the message when the files are missing). Each module named
here holds its dataset and the functions that read it: the MIST main sequence (`isochrones`), the white-dwarf cooling
sequences with their layers (`white_dwarf_models`) and the neutron-star equations of state with the outer-crust
compositions (`neutron_star_matter`), the ozone absorption cross sections of the radiation model (`ozone_uv`) and the infrared opacities
of the greenhouse model: the ExoMol correlated-k tables of water, carbon dioxide and methane (`ir_ktables`) and the HITRAN collision-induced absorption (`hitran_cia`), the Harris globular cluster catalog (`harris_clusters`), the NOAA solar proton events
(`noaa_proton_events`), the hydrogen-helium equation of state of giant planets (`hhe_eos`), the mean dwarf spectral-type sequence of Pecaut and Mamajek (`mamajek_sequence`) the Baraffe 2015 isochrones of low-mass stars and brown dwarfs (`baraffe_isochrones`) and the model-atmosphere spectra of stars (`stellar_spectra`). The data of the atmosphere group (AD1 to AD8) is added the same way, each by the
item that reads it; every entry needs its verified source (docs/design/research/gap-filling-plan.md 6.5) and a stand-in row
(section 6.6). The generator requires every dataset here (GEN.590): `missing` is the check bulk generation starts with.
"""

from planetgen.datasets import (baraffe_isochrones, fetch, harris_clusters, hhe_eos, hitran_cia, ir_ktables, isochrones, mamajek_sequence,
                                neutron_star_matter, noaa_proton_events, ozone_uv, stellar_spectra, white_dwarf_models)

ALL = (isochrones.DATASET, white_dwarf_models.DATASET, neutron_star_matter.DATASET, ozone_uv.DATASET, ir_ktables.DATASET, hitran_cia.DATASET,
       harris_clusters.DATASET, noaa_proton_events.DATASET, baraffe_isochrones.DATASET, stellar_spectra.DATASET, mamajek_sequence.DATASET, hhe_eos.DATASET)
"""tuple: Every `fetch.Dataset` the setup step downloads, in the order it downloads them."""


def missing():
    """
    The datasets whose files are not all downloaded and verified in the cache folder, each as the `fetch.Located` that says why
    (empty when everything is there). It opens no connection: setup downloads, this only looks.
    """
    return [located for located in (fetch.locate(dataset) for dataset in ALL) if located.missing]
