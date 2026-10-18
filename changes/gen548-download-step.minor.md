**Headline:** Static data files are downloaded once at setup and checked, with a labelled fallback when a host is down

### Changed

`update.sh` and `install.sh` now run one shared download step (`python3 -m planetgen.cli.datasets fetch`) that fetches the static data files the generator reads into `/var/cache/planetgen/datasets`, outside the checkout. Each request sends the standard `User-Agent` (software, release and API version), and each file must match the size and sha256 the catalogue states, otherwise it is not kept. A host that does not answer only prints a warning; the generator then uses a stand-in that says it is one. Nothing downloads while the site or the generator runs. The catalogue is empty until the isochrone grids, white-dwarf layers and neutron-star equation of state are added.
