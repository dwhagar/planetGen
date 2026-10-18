**Headline:** The generator's stellar isochrone grid, white-dwarf models and neutron-star matter tables are downloaded once at setup and checked

### Changed

`update.sh` and `install.sh` now run one shared download step (`python3 -m planetgen.cli.datasets fetch`) that fetches static data files into `/var/cache/planetgen/datasets`, outside the checkout. Each request sends the standard `User-Agent` (software, release and API version), is retried if the host does not answer, and each file must match the size and sha256 the catalogue states, otherwise it is not kept. A host that does not answer only prints a warning; the generator then uses a stand-in that says it is one. Nothing downloads while the site or the generator runs.

Three datasets use it: the MIST v1.2 main-sequence isochrones (luminosity, temperature and radius against mass, age and [Fe/H], and main-sequence lifetimes; a 221 MB archive reduced once to 5 MB), the Bédard et al. 2020 white-dwarf cooling sequences (radius and layer masses), and 21 neutron-star equations of state with the outer-crust composition for 20 nuclear mass models. Each has a reader in `planetgen.datasets`; the first download is about 225 MB and takes a minute or two.
