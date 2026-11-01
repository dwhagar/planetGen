# planetGen

**Version:** 8.0.3058 &middot; [Changelog](CHANGELOG.md) &middot; [Repository](https://github.com/dwhagar/planetGen) &middot; License: [CC0 1.0 Universal](LICENSE.md)

planetGen generates a galaxy: stars, star systems, planets, moons,
asteroid belts and exotic phenomena, placed in a physically modeled
spiral galaxy and saved to a MySQL database. A website on top of it lets
you browse and search the galaxy, fly around it on 3D maps, and copy any
system's write-up as Markdown. It was built for
the Molten Aether FFRP game.

**To install it, see [INSTALL.md](INSTALL.md).**

## Headline features

Every release in the [changelog](CHANGELOG.md) opens with its one best
feature. These are the ones that shaped the program, oldest first.

| Release | Date | Headline feature |
|---|---|---|
| 1.0.0 | 2024-06-25 | First fully working end-to-end planet generator validated against Earth |
| 3.0.0 | 2026-07-07 | Binary star systems |
| 5.0.0 | 2026-09-05 | Space sectors placing star systems in 3D positions with realistic spacing and save/load |
| 5.2.1 | 2026-09-06 | Web interface for browsing saved databases, with Apache deployment and one-shot installer |
| 5.4.8 | 2026-09-09 | Galaxy-wide density skeleton with lazy, visit-triggered sector generation |
| 5.5.0 | 2026-09-09 | MySQL backend replaces SQLite, plus a paginated, validated API |
| 5.8.0 | 2026-09-10 | NAV feature: course, distance, warp times, and optimal routes between systems |
| 5.29.0 | 2026-09-12 | Exotic phenomena: black holes, neutron stars, nebulae, supernova remnants, rogue planets, comets |
| 5.37.0 | 2026-09-18 | Publish systems and sectors to Wiki.js or MediaWiki from the web interface |
| 5.46.13 | 2026-09-23 | Interactive 3D Galaxy Map with live viewport queries and copyable generate commands |
| 5.47.0 | 2026-09-24 | Galaxy Map loads space as cached cube tiles, preventing server crashes |
| 6.1.0 | 2026-09-30 | Generate, plan and reset the galaxy from the web interface |
| 7.23.0 | 2026-10-01 | Stars drawn from real physics: mass function, age and stellar evolution |
| 7.44.0 | 2026-10-01 | Galaxy Map drill-down from whole galaxy to single sectors with linkable stages |
| 7.49.0 | 2026-10-01 | Dominant species, civilizations and polities that claim territories in 3D |
| 7.143.452 | 2026-10-02 | One 128-bit seed per galaxy gives reproducible generation at any worker count |
| 7.225.618 | 2026-10-07 | Milky Way galaxy model with thin disk, thick disk and boxy bar bulge |
| 7.339.675 | 2026-10-08 | 3D system view with orbits, belts, free camera and time control |
| 8.0.711 | 2026-10-09 | Phase 0 complete: saved sector paths, nearby-object search and re-centerable 3D maps |
| 8.0.912 | 2026-10-10 | Every planet and moon gets a habitability score and equipment-needed rating |

### Coming next

Each upcoming release has one headline feature, and the project plan
traces exactly which items each one needs (the plan pages list them in
build order).

| Release | Headline feature | Status |
|---|---|---|
| 8.0 | The orbital update system: the galaxy moves in real time under real gravity | In progress |
| 8.1 | Infinite zoom: fly from the whole galaxy down to a single moon on one map | Planned |
| 8.5 | Courses that bend around gravity wells, with waypoints, editable and saved courses | Planned |
| 9.0 | The gravity map: the pull, depth and tidal strength of a sector as a heat map | Planned |

## What it does today

- **Stars drawn as physics**: a mass from the Kroupa initial mass
  function, an age from the galaxy's star-formation history, then how
  the star has evolved by now (main sequence, subgiant, giant,
  supergiant or white dwarf). The result matches the real galaxy: about
  three in four stars are M dwarfs and supergiants are about one in a
  million. A star can also be pinned to a type such as `G2V`.
- **Binary systems**, both close (circumbinary) pairs and wide pairs where
  each star keeps its own planets. The two stars are born together, and
  orbits respect published stability limits.
- **Planets and moons** from 25 planet classes, placed in physically
  plausible orbits and named from their system (`Sol III`, `Sol IIIa`).
  They follow their star's history: a very young star has only belts,
  and a giant has swallowed its inner planets.
- **Habitable worlds and life**: which life chemistries a star supports,
  and a speculative evolutionary timeline for each living world.
- **A whole galaxy**: a density model of a spiral galaxy, with young,
  old and bulge star populations, cut into 4 pc (13 ly) sectors on a
  grid of spherical shells, bands and wedges around a round core. Planning the galaxy also places every bright star
  (500 times the Sun's luminosity or more, about 60 million in a Milky
  Way) so the spiral arms exist from the start. Sectors are then filled
  a neighborhood at a time, where and when you need them, with star
  ages that follow where each sector sits.
- **Exotic phenomena** at research-based interstellar rates: black
  holes, neutron stars, nebulae (classed A-Q, including planetary
  nebulae) and supernova remnants (R-W) with what fills them, rogue
  planets in four mass bins plus free-floating brown dwarfs,
  interstellar comets, and runaway and hypervelocity stars. Nebulae come
  with the stars that make them: H II regions around O and B stars, a
  white dwarf at the heart of each planetary nebula. Each object records
  the cloud it sits inside (which squeezes the heliosphere of systems
  within it), its sector octant and its three nearest systems. Every
  galaxy has a supermassive black hole at its center, or one time in ten
  an active quasar.
- **Facilities**: starbases, stations, colonies and outposts on planets
  and moons, in orbit, in belts or in open space, placed by an admin
  from a system's page or through the API.
- **Population (optional, off by default)**: a pass that names the
  dominant species of worlds with complex life, dates civilizations and
  their eras, and gives each spacefaring species a polity whose
  territory claims the systems around it.
- **Time passes**: a monthly update moves planets and moons along their
  orbits and every system and phenomenon along its galactic orbit,
  carrying anything that drifts into another sector over to it.
- **A website** to browse it all: sector and system pages, faceted
  search, 3D maps of the galaxy and of each sector, a System Map that
  measures routes around planets and stars, a NAV page that plots
  courses between systems, a JSON API, and admin pages to generate the
  galaxy from the browser.

## Requirements

- Python 3.12 or later, with numba (the compiled numeric kernels; `install.sh` and `update.sh` install it, and there is no fallback without it). This replaced the old Python 3.9 minimum with numba (OPS.42).
- MySQL 8.0.19+ (MariaDB is not supported).
- Redis (any recent version): the Generate page's jobs, the work queue behind `--workers` above 1, the rate limits and the login lockouts run on it. `install.sh` and `update.sh` install and start it on Debian and Ubuntu and check its settings.
- For the website, a web server: Apache2 with mod_wsgi on Debian or
  Ubuntu (the reference setup), nginx or Caddy on Linux, or nginx on
  macOS.
- Install and update scripts: `install.sh` and `update.sh` on Linux and
  macOS.
- Windows: see [docs/deployment/WINDOWS.md](docs/deployment/WINDOWS.md).

[INSTALL.md](INSTALL.md) lists the full requirements and walks through
the install.

## External files

The generator reads a few large published data sets that are not in this repository (raw data sets are never committed). `install.sh` and `update.sh` download them once at setup, as root, with `python3 -m planetgen.cli.datasets fetch`; nothing downloads while the site or the generator runs. Each file is requested over https, must have the exact size and sha256 the catalogue (`src/planetgen/datasets/catalogue.py`) states, and is kept by an atomic rename. They go in the cache folder `/var/cache/planetgen/datasets` (one sub-folder per set below), or in the folder `PLANETGEN_DATA_DIR` names. The licence and citation of each set are in [DATA-LICENSE.md](DATA-LICENSE.md).

- **Fetched at setup:** every set listed below. **Provided by hand:** none; the program reads no other external file.
- **See where each stands:** `python3 -m planetgen.cli.datasets status` fetches nothing, says which sets are present, prints the URL of every missing file, and exits 1 while one is missing. `fetch` prints one line per file with its URL and what happened to it (`fetched`, `kept` or `failed` with the reason), and `update.sh` and `install.sh` show those lines.
- **When a file is missing:** a bulk run (`galaxy`, `plan`, `population`, `sector` with more than one sector, the Generate page's jobs and the Sector page's neighbourhood button) refuses to start, writes nothing, and names each missing set, why, every file with its URL and the setup command (`python3 -m planetgen.cli.datasets fetch KEY`). Unless it is told, in so many words, to disregard missing files (`--disregard-missing-data` on any subcommand, the "Disregard missing data files" box on the Generate page, or `"disregard_missing_data": true` in the neighbourhood request): it then starts, uses the labelled fallback listed under each set below, and prints the fallbacks in the run header. One system, one sector or one phenomenon is not gated.
- **Installing a set by hand** (the host is down, or the machine has no network): download each URL listed under the set on a machine that can, check the file's size and sha256 against the catalogue entry, copy it into the folder `<cache folder>/<set key>/` under the name shown, make the folder readable for the web server's user, and run `python3 -m planetgen.cli.datasets fetch KEY`: a file already there with the right size and sha256 is kept and not requested again (one that fails the check is requested again), and a set that is reduced to a smaller file (the "kept after setup" file) is reduced from it.

The list below is written from the catalogue by `scripts/readme_external_files.py --write`, and a test fails until it matches, so a new set cannot skip it.

<!-- external-files:start (written by scripts/readme_external_files.py) -->

### `isochrones`: MIST v1.2 isochrones, main sequence (v/vcrit 0.4, solar-scaled)

- Source: Choi et al. 2016 (ApJ 823, 102) and Dotter 2016 (ApJS 222, 8), MIST v1.2, mist.science.
- Download: 1 file, 221.5 MB. Kept after setup: `main_sequence.npz`, `post_main_sequence.npz` (the download is deleted once it is reduced).
- Fetch: `python3 -m planetgen.cli.datasets fetch isochrones`.
- Replaces: the Baraffe et al. 2015 track at solar metallicity and the fitted relations of the stellar-composition note.
- With `--disregard-missing-data`: No MIST point, so no luminosity, temperature, radius or lifetime for a star (an open item on the star).
- Files (URL, then the file name it is kept under in the dataset's folder):
  - <https://mist.science/data/tarballs_v1.2/MIST_v1.2_vvcrit0.4_basic_isos.txz> as `MIST_v1.2_vvcrit0.4_basic_isos.txz`

### `white_dwarf_models`: Bédard et al. 2020 white-dwarf cooling sequences

- Source: Bédard, Bergeron, Brassard and Fontaine 2020, ApJ 901, 93 (doi:10.3847/1538-4357/abafbe).
- Download: 1 file, 2.5 MB. All kept as downloaded.
- Fetch: `python3 -m planetgen.cli.datasets fetch white_dwarf_models`.
- Replaces: the Chandrasekhar-type law of the stellar-composition note, an even carbon-oxygen split and a fixed layer split.
- With `--disregard-missing-data`: The radius of the zero-temperature Chandrasekhar table (Sirius B check, no temperature or layers), no layer fractions, and the Mestel power law of the cooling (the stand-in that came before the sequences: luminosity of the cooling age, one radius per mass).
- Files (URL, then the file name it is kept under in the dataset's folder):
  - <https://warwick.ac.uk/fac/sci/physics/research/astro/people/antoinebedard/models/sequences.tar> as `sequences.tar`

### `neutron_star_matter`: Neutron-star equations of state (Malik et al. 2024) and outer-crust compositions (Servais and Chamel 2026)

- Source: Malik, Pais and Providência 2024 (A&A 689, A242) and Servais and Chamel 2026 (A&A 710, A322), Zenodo.
- Download: 22 files, 2.0 MB. All kept as downloaded.
- Fetch: `python3 -m planetgen.cli.datasets fetch neutron_star_matter`.
- Replaces: the nickel-to-krypton crust table and the 10 km radius floor of the stellar-composition note.
- With `--disregard-missing-data`: The one outer-crust sequence of Chamel and Haensel 2008 Table 2, in place of the 20 nuclear mass models.
- Files (URL, then the file name it is kept under in the dataset's folder):
  - <https://zenodo.org/records/20085135/files/compositions.zip?download=1> as `compositions.zip`
  - <https://zenodo.org/records/10950421/files/EOS1.txt?download=1> as `EOS1.txt`
  - <https://zenodo.org/records/10950421/files/EOS2.txt?download=1> as `EOS2.txt`
  - <https://zenodo.org/records/10950421/files/EOS3.txt?download=1> as `EOS3.txt`
  - <https://zenodo.org/records/10950421/files/EOS4.txt?download=1> as `EOS4.txt`
  - <https://zenodo.org/records/10950421/files/EOS5.txt?download=1> as `EOS5.txt`
  - <https://zenodo.org/records/10950421/files/EOS6.txt?download=1> as `EOS6.txt`
  - <https://zenodo.org/records/10950421/files/EOS7.txt?download=1> as `EOS7.txt`
  - <https://zenodo.org/records/10950421/files/EOS8.txt?download=1> as `EOS8.txt`
  - <https://zenodo.org/records/10950421/files/EOS9.txt?download=1> as `EOS9.txt`
  - <https://zenodo.org/records/10950421/files/EOS10.txt?download=1> as `EOS10.txt`
  - <https://zenodo.org/records/10950421/files/EOS11.txt?download=1> as `EOS11.txt`
  - <https://zenodo.org/records/10950421/files/EOS12.txt?download=1> as `EOS12.txt`
  - <https://zenodo.org/records/10950421/files/EOS13.txt?download=1> as `EOS13.txt`
  - <https://zenodo.org/records/10950421/files/EOS14.txt?download=1> as `EOS14.txt`
  - <https://zenodo.org/records/10950421/files/EOS15.txt?download=1> as `EOS15.txt`
  - <https://zenodo.org/records/10950421/files/EOS16.txt?download=1> as `EOS16.txt`
  - <https://zenodo.org/records/10950421/files/EOS17.txt?download=1> as `EOS17.txt`
  - <https://zenodo.org/records/10950421/files/EOS18.txt?download=1> as `EOS18.txt`
  - <https://zenodo.org/records/10950421/files/EOS19.txt?download=1> as `EOS19.txt`
  - <https://zenodo.org/records/10950421/files/EOS20.txt?download=1> as `EOS20.txt`
  - <https://zenodo.org/records/10950421/files/EOS21.txt?download=1> as `EOS21.txt`

### `ozone_uv`: Ozone absorption cross sections, 195 to 830 nm at 218, 243 and 295 K, and the Chappuis band at 220 to 298 K (Brion, Daumont, Malicet; Burkholder and Talukdar; MPI-Mainz atlas)

- Source: Keller-Rudek et al. 2013 (ESSD 5, 365), MPI-Mainz UV/VIS Spectral Atlas snapshot 2013-06-24, Zenodo 10.5281/zenodo.6951, CC BY 4.0.
- Download: 1 file, 126.3 MB. Kept after setup: `ozone_uv.npz` (the download is deleted once it is reduced).
- Fetch: `python3 -m planetgen.cli.datasets fetch ozone_uv`.
- Replaces: the 20 nm band means of the gas table for ozone, which understate the transmission in the UV-B and misplace the UV index.
- With `--disregard-missing-data`: The 35 band means of the gas table (physics/data/gases.json) for ozone.
- Files (URL, then the file name it is kept under in the dataset's folder):
  - <https://zenodo.org/records/6951/files/uv-vis-spectral-atlas-2013-06-24.zip?download=1> as `uv-vis-spectral-atlas-2013-06-24.zip`

### `ir_ktables`: Infrared correlated-k tables of H2O, CO2 and CH4 (ExoMolOP, R1000, line lists POKAZATEL, UCL-4000 and YT34to10)

- Source: Chubb et al. 2021 (A&A 646, A21), ExoMolOP R1000 ARCiS k-tables from exomol.com.
- Download: 3 files, 1257.8 MB. Kept after setup: `ir_ktables.npz` (the download is deleted once it is reduced).
- Fetch: `python3 -m planetgen.cli.datasets fetch ir_ktables`.
- Replaces: the band-strength and wing stand-ins of the grey infrared opacity: plateau opacities per band centre of the gas table.
- With `--disregard-missing-data`: None at generation: the fitted greenhouse table (physics/greenhouse_table.py) carries the answer; the tables are read to fit and check it.
- Files (URL, then the file name it is kept under in the dataset's folder):
  - <https://exomol.com/db/H2O/1H2-16O/POKAZATEL/1H2-16O__POKAZATEL__R1000_0.3-50mu.ktable.ARCiS.fits.gz> as `ir_ktable_h2o.fits.gz`
  - <https://exomol.com/db/CO2/12C-16O2/UCL-4000/12C-16O2__UCL-4000.R1000_0.3-50mu.ktable.ARCiS.fits.gz> as `ir_ktable_co2.fits.gz`
  - <https://exomol.com/db/CH4/12C-1H4/YT34to10/12C-1H4__YT34to10.R1000_0.3-50mu.ktable.ARCiS.fits.gz> as `ir_ktable_ch4.fits.gz`

### `ir_cia`: Collision-induced absorption of N2, O2, CO2, H2, He and CH4 pairs (HITRAN CIA, 2011 to 2024 releases)

- Source: Karman et al. 2019 (Icarus 328, 160), HITRAN collision-induced absorption files from hitran.org.
- Download: 11 files, 195.4 MB. Kept after setup: `ir_cia.npz` (the download is deleted once it is reduced).
- Fetch: `python3 -m planetgen.cli.datasets fetch ir_cia`.
- Replaces: the coarse pair coefficients of the grey infrared opacity, read at six temperatures from a part-03 table.
- With `--disregard-missing-data`: None at generation: the fitted greenhouse table (physics/greenhouse_table.py) carries the answer; the pairs are read to fit and check it.
- Files (URL, then the file name it is kept under in the dataset's folder):
  - <https://hitran.org/data/CIA/main/N2-N2_2021.cia> as `N2-N2_2021.cia`
  - <https://hitran.org/data/CIA/main/O2-O2_2024.cia> as `O2-O2_2024.cia`
  - <https://hitran.org/data/CIA/main/O2-N2_2024.cia> as `O2-N2_2024.cia`
  - <https://hitran.org/data/CIA/main/CO2-CO2_2024.cia> as `CO2-CO2_2024.cia`
  - <https://hitran.org/data/CIA/main/N2-CH4_2024.cia> as `N2-CH4_2024.cia`
  - <https://hitran.org/data/CIA/main/N2-H2_2024.cia> as `N2-H2_2024.cia`
  - <https://hitran.org/data/CIA/main/CO2-CH4_2024.cia> as `CO2-CH4_2024.cia`
  - <https://hitran.org/data/CIA/main/CO2-H2_2024.cia> as `CO2-H2_2024.cia`
  - <https://hitran.org/data/CIA/main/H2-H2_2011.cia> as `H2-H2_2011.cia`
  - <https://hitran.org/data/CIA/main/H2-CH4_norm_2011.cia> as `H2-CH4_norm_2011.cia`
  - <https://hitran.org/data/CIA/main/H2-He_2011.cia> as `H2-He_2011.cia`

### `harris_clusters`: Harris catalog of parameters for Milky Way globular clusters (2010 edition)

- Source: Harris 1996, AJ 112, 1487 (2010 edition), McMaster University.
- Download: 1 file, 49 kB. All kept as downloaded.
- Fetch: `python3 -m planetgen.cli.datasets fetch harris_clusters`.
- Replaces: the globular cluster laws of tuning.py, recalled from Harris 1996.
- With `--disregard-missing-data`: The recalled medians of the cluster laws in tuning.py.
- Files (URL, then the file name it is kept under in the dataset's folder):
  - <https://physwww.mcmaster.ca/~harris/mwgc.dat> as `mwgc.dat`

### `noaa_proton_events`: NOAA solar proton events above 10 MeV, 1976 to 2017

- Source: NOAA Space Environment Services Center, Solar Proton Events Affecting the Earth Environment, via NASA SDAC.
- Download: 1 file, 45 kB. All kept as downloaded.
- Fetch: `python3 -m planetgen.cli.datasets fetch noaa_proton_events`.
- Replaces: the one solar particle dose and the exponential in rigidity of the radiation model.
- With `--disregard-missing-data`: The stand-in dose and rigidity law of physics/radiation.py.
- Files (URL, then the file name it is kept under in the dataset's folder):
  - <https://umbra.nascom.nasa.gov/SEP/seps.html> as `seps.html`

### `baraffe_isochrones`: Baraffe et al. 2015 (BHAC15) isochrones of low-mass stars and brown dwarfs

- Source: Baraffe, Homeier, Allard and Chabrier 2015, A&A 577, A42.
- Download: 1 file, 65 kB. All kept as downloaded.
- Fetch: `python3 -m planetgen.cli.datasets fetch baraffe_isochrones`.
- Replaces: the nine brown dwarf radius points and the power-law extension of the main sequence under 0.1 solar masses (the three pre-main-sequence constants of tuning.py stay as fitted to it).
- With `--disregard-missing-data`: The fitted pre-main-sequence constants of tuning.py; a brown dwarf's radius and the main sequence under 0.1 solar masses are left out.
- Files (URL, then the file name it is kept under in the dataset's folder):
  - <https://perso.ens-lyon.fr/isabelle.baraffe/BHAC15dir/BHAC15_iso.2mass> as `BHAC15_iso.2mass`

### `stellar_spectra`: Model-atmosphere spectra: Castelli and Kurucz ATLAS9 fluxes (3500 to 50000 K) and Husser et al. PHOENIX-ACES (2300 to 3500 K)

- Source: Castelli and Kurucz 2003 (IAU Symp. 210, A20), Husser et al. 2013 (A&A 553, A6).
- Download: 23 files, 156.4 MB. Kept after setup: `stellar_spectra.npz` (the download is deleted once it is reduced).
- Fetch: `python3 -m planetgen.cli.datasets fetch stellar_spectra`.
- Replaces: the blackbody of the stellar spectrum with a seeded UV-C share (R17) and the BT-Settl chromaticity offsets.
- With `--disregard-missing-data`: The blackbody of the effective temperature with the seeded UV shares of physics/stellar_spectrum.py.
- Files (URL, then the file name it is kept under in the dataset's folder):
  - <https://wwwuser.oats.inaf.it/fiorella.castelli/grids/gridp00k2odfnew/fp00k2odfnew.pck> as `fp00k2odfnew.pck`
  - <https://ftp.gwdg.de/pub/misc/phoenix/HiResFITS/WAVE_PHOENIX-ACES-AGSS-COND-2011.fits> as `WAVE_PHOENIX-ACES-AGSS-COND-2011.fits`
  - <https://ftp.gwdg.de/pub/misc/phoenix/HiResFITS/PHOENIX-ACES-AGSS-COND-2011/Z-0.0/lte02300-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits> as `lte02300-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits`
  - <https://ftp.gwdg.de/pub/misc/phoenix/HiResFITS/PHOENIX-ACES-AGSS-COND-2011/Z-0.0/lte02300-5.00-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits> as `lte02300-5.00-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits`
  - <https://ftp.gwdg.de/pub/misc/phoenix/HiResFITS/PHOENIX-ACES-AGSS-COND-2011/Z-0.0/lte02300-5.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits> as `lte02300-5.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits`
  - <https://ftp.gwdg.de/pub/misc/phoenix/HiResFITS/PHOENIX-ACES-AGSS-COND-2011/Z-0.0/lte02500-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits> as `lte02500-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits`
  - <https://ftp.gwdg.de/pub/misc/phoenix/HiResFITS/PHOENIX-ACES-AGSS-COND-2011/Z-0.0/lte02500-5.00-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits> as `lte02500-5.00-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits`
  - <https://ftp.gwdg.de/pub/misc/phoenix/HiResFITS/PHOENIX-ACES-AGSS-COND-2011/Z-0.0/lte02500-5.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits> as `lte02500-5.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits`
  - <https://ftp.gwdg.de/pub/misc/phoenix/HiResFITS/PHOENIX-ACES-AGSS-COND-2011/Z-0.0/lte02700-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits> as `lte02700-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits`
  - <https://ftp.gwdg.de/pub/misc/phoenix/HiResFITS/PHOENIX-ACES-AGSS-COND-2011/Z-0.0/lte02700-5.00-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits> as `lte02700-5.00-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits`
  - <https://ftp.gwdg.de/pub/misc/phoenix/HiResFITS/PHOENIX-ACES-AGSS-COND-2011/Z-0.0/lte02700-5.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits> as `lte02700-5.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits`
  - <https://ftp.gwdg.de/pub/misc/phoenix/HiResFITS/PHOENIX-ACES-AGSS-COND-2011/Z-0.0/lte02900-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits> as `lte02900-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits`
  - <https://ftp.gwdg.de/pub/misc/phoenix/HiResFITS/PHOENIX-ACES-AGSS-COND-2011/Z-0.0/lte02900-5.00-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits> as `lte02900-5.00-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits`
  - <https://ftp.gwdg.de/pub/misc/phoenix/HiResFITS/PHOENIX-ACES-AGSS-COND-2011/Z-0.0/lte02900-5.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits> as `lte02900-5.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits`
  - <https://ftp.gwdg.de/pub/misc/phoenix/HiResFITS/PHOENIX-ACES-AGSS-COND-2011/Z-0.0/lte03100-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits> as `lte03100-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits`
  - <https://ftp.gwdg.de/pub/misc/phoenix/HiResFITS/PHOENIX-ACES-AGSS-COND-2011/Z-0.0/lte03100-5.00-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits> as `lte03100-5.00-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits`
  - <https://ftp.gwdg.de/pub/misc/phoenix/HiResFITS/PHOENIX-ACES-AGSS-COND-2011/Z-0.0/lte03100-5.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits> as `lte03100-5.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits`
  - <https://ftp.gwdg.de/pub/misc/phoenix/HiResFITS/PHOENIX-ACES-AGSS-COND-2011/Z-0.0/lte03300-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits> as `lte03300-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits`
  - <https://ftp.gwdg.de/pub/misc/phoenix/HiResFITS/PHOENIX-ACES-AGSS-COND-2011/Z-0.0/lte03300-5.00-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits> as `lte03300-5.00-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits`
  - <https://ftp.gwdg.de/pub/misc/phoenix/HiResFITS/PHOENIX-ACES-AGSS-COND-2011/Z-0.0/lte03300-5.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits> as `lte03300-5.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits`
  - <https://ftp.gwdg.de/pub/misc/phoenix/HiResFITS/PHOENIX-ACES-AGSS-COND-2011/Z-0.0/lte03500-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits> as `lte03500-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits`
  - <https://ftp.gwdg.de/pub/misc/phoenix/HiResFITS/PHOENIX-ACES-AGSS-COND-2011/Z-0.0/lte03500-5.00-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits> as `lte03500-5.00-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits`
  - <https://ftp.gwdg.de/pub/misc/phoenix/HiResFITS/PHOENIX-ACES-AGSS-COND-2011/Z-0.0/lte03500-5.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits> as `lte03500-5.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits`

### `mamajek_dwarf_sequence`: Mamajek mean dwarf stellar colour and effective temperature sequence (2022.04.16)

- Source: Pecaut and Mamajek 2013, ApJS 208, 9 (doi:10.1088/0067-0049/208/1/9), version 2022.04.16 of E. Mamajek's table.
- Download: 1 file, 56 kB. All kept as downloaded.
- Fetch: `python3 -m planetgen.cli.datasets fetch mamajek_dwarf_sequence`.
- Replaces: the subclasses linear in temperature inside a letter, with a K/M boundary at 3,700 K.
- With `--disregard-missing-data`: nothing falls back, the run is still refused.
- Files (URL, then the file name it is kept under in the dataset's folder):
  - <https://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt> as `EEM_dwarf_UBVIJHK_colors_Teff.txt`

<!-- external-files:end -->

## Usage

### The website

Once installed, the site's front page leads to:

- **Galaxy Map**: the whole galaxy in 3D, with every bright star as a
  tiny point in a big soft glow and nebulae as translucent clouds. It
  opens on a drill-down: pick a slab of the galaxy, then a block, then
  smaller blocks down to single sectors, with a breadcrumb, keyboard
  keys, and a link for every stage so Back and Forward work. Type a
  sector, a `shell/band/wedge` address, coordinates or a name to fly
  there. A NAV course can be drawn on it, and once population has run,
  a Territories overlay shows who holds what. As an admin, click an
  empty sector to generate it, a neighborhood of a radius you choose,
  its column or its whole shell, or generate a whole small block or one
  of its layers.
- **Sectors** and **Systems**: every generated sector and every system,
  50 to a page. A sector has a 3D Sector Map (with the same Generate
  buttons for admins, and Nav from/to on any system) and a Contents
  list; a system has its stars, planets, moons and facilities, the cloud
  it sits in, its nearest systems, a System Map, and its write-up
  ready to copy as Markdown. Phenomena have pages too, with their class, what fills
  them, and a view that suits each one (a neutron star's spinning beams,
  a black hole's accretion disk).
- **Classes**: a reference page for every kind of class the generator
  gives out: stars, planets, nebulae, asteroid fields and more.
- **Species** and **Polities** (once the optional population pass has
  run): every species with its homeworld and era, and every polity with
  the systems it holds. They stay hidden until there is something to
  show.
- **NAV**: a course and distance between any two systems or phenomena,
  picked from their pages, the Sector Map or the Galaxy Map, with "Show
  on Galaxy Map" to see it drawn there.
- **Search**: filter by object type, star class, planet class, life
  chemistry, phenomenon and phenomenon class, or search by name from the
  box in every page header.
- **Admin** (after logging in): **Generate** plans the galaxy and
  generates sectors as background jobs you can watch and cancel;
  **Generate a one-off system** makes a single system with every
  command-line option and never saves it; plus your account, API keys
  and server stats. Times show in your own time zone.

[`docs/reference/html-interface.md`](docs/reference/html-interface.md) describes every page
and [`docs/reference/api.md`](docs/reference/api.md) the JSON API.

### The command line

`planetgen.cli.generate`, in the checkout, is the one entry point for
every generator. On Linux and macOS the installer adds a `planetgen`
command that runs it. Each subcommand saves what it makes to the database
unless told otherwise:

```bash
planetgen plan         # plan the galaxy's shape (once, before galaxy)
planetgen scatter      # draw one scatter on the plan: --mass, --luminosity, --phenomena (before galaxy)
planetgen galaxy       # generate sectors in the galaxy
planetgen sector       # one or more standalone sectors
planetgen system       # one star system
planetgen phenomenon   # one exotic phenomenon
planetgen population   # optional: species, civilizations and polities
planetgen territories  # optional: which polity owns which system (after population)
```

Common examples:

```bash
# A random neighborhood: a random start and every sector within 100 ly of it
planetgen galaxy

# One whole shell, one band of it, or one sector (shell 12, the equatorial
# band, wedge 5)
planetgen galaxy --shell 12 --limit 500
planetgen galaxy --shell 12 --band 0
planetgen galaxy --shell 12 --band 0 --wedge 5

# A span of shells and bands (negative bands lie south of the plane)
planetgen galaxy --shells 10:12 --bands=-2:2 --limit 500

# A quick test galaxy: scatter fewer pre-placed bright stars
planetgen scatter --mass --bright-star-min-luminosity 5000
planetgen scatter --luminosity

# A single system with the odds of a habitable world raised (it may come out without one), printed as Markdown on request, not saved
planetgen system +habitable_world --format markdown --output -

# A system built from a specification file
planetgen system --system-file examples/systems/solar_system.json
```

Most options use a `+name`/`-name` form: `+moons` forces moons,
`-moons` forbids them, and leaving it off leaves it to chance.
`planetgen <command> --help` lists every option, and
[`docs/reference/cli.md`](docs/reference/cli.md) is the full reference, including the
system specification file format and the incompatible option pairs.

### Other tools

Each runs as `python3 -m <module>` from planetGen's own virtual environment
(`/usr/local/planetgen/venv/bin/python -m <module>` on a server `install.sh`
set up; a bare `python3` there finds none of the libraries); `--help` shows
its options:

- **`planetgen.cli.orbits`** advances the galaxy by the time elapsed since
  the last run: planets, moons and orbiting facilities around their
  hosts, and every system and phenomenon along its galactic orbit. Run
  it monthly; the maintenance timers in
  [INSTALL.md](INSTALL.md#scheduled-maintenance) do.
- **`planetgen.cli.migrate`** brings the database schema up to date. The install
  and update scripts run it for you.
- **`planetgen.cli.reset`** deletes all generated galaxy data (asks for the
  database name first) so you can start a new galaxy.
- **`planetgen.cli.query`** searches the database from the terminal.

## Documentation

| Document | What it covers |
|---|---|
| [INSTALL.md](INSTALL.md) | Installing, updating and the example configs |
| [docs/reference/cli.md](docs/reference/cli.md) | Every command-line option |
| [docs/reference/config.md](docs/reference/config.md) | Every `config.json` setting |
| [docs/deployment/](docs/deployment/README.md) | Hosting guides per platform |
| [docs/deployment/server-checklist.md](docs/deployment/server-checklist.md) | Confirming the site loads after a deploy or update |
| [docs/deployment/WINDOWS.md](docs/deployment/WINDOWS.md) | Windows is unsupported; the simplest setup, WSL2 |
| [docs/reference/html-interface.md](docs/reference/html-interface.md) | How the website works |
| [docs/reference/api.md](docs/reference/api.md) | The JSON API |
| [docs/reference/atmosphere-colour.md](docs/reference/atmosphere-colour.md) | How a planet's atmosphere colour is worked out |
| [docs/reference/system-file-format.md](docs/reference/system-file-format.md) | The JSON file that describes one star system for `planetgen system` |
| [docs/reference/example-systems.md](docs/reference/example-systems.md) | The example system files and what each one tests |
| [docs/reference/database-schema.md](docs/reference/database-schema.md) | The database tables and migrations |
| [docs/design/architecture/architecture.md](docs/design/architecture/architecture.md) | How the program fits together: every file and the main flows |
| [docs/design/architecture/architecture.html](docs/design/architecture/architecture.html) | The same as an illustrated page to open in a browser: modules, flows, data and deployment (a snapshot of 2026-10-17, rebuilt from the architecture.md above) |
| [docs/generated/architecture-map.html](docs/generated/architecture-map.html) | The code map to open in a browser and drill down: the project's areas, their modules, each module's functions with what calls what, and the main flows traced from the code (generated by `scripts/build_architecture_map.py`) |
| [docs/design/](docs/design/) | Design notes: the galaxy model, coordinates, navigation, phenomena, and why each choice was made |
| [docs/TODO.md](docs/TODO.md) | Open work, and the index of the phase plans |
| [docs/plan/](docs/plan/) | The phase plans: order, dependencies and research notes for each phase of the open work |
| [docs/generated/version-plan.html](docs/generated/version-plan.html) | Every open item by version, phase, thread and dependency, as one page to open in a browser |
| [docs/generated/todo-tree.html](docs/generated/todo-tree.html) | Every open item as a tree, each under the prerequisite that finishes last, with search and phase and bug filters |
| [docs/generated/todo-reference.html](docs/generated/todo-reference.html) | Every TODO ID ever issued, explained in full, with search and filters (all three pages are rebuilt by `python scripts/build_todo_docs.py`) |
| [docs/reference/logging.md](docs/reference/logging.md) | The form of every log line (time, severity, component, stage markers) and where each log is found |
| [docs/reference/testing.md](docs/reference/testing.md) | Running the tests |
| [docs/reference/ci-runners.md](docs/reference/ci-runners.md) | Running the scheduled CI on your own computers (optional) |
| [CHANGELOG.md](CHANGELOG.md) | What changed in each release |

## A note on accuracy

planetGen aims for plausibility, not a simulation. Forcing extremes
(a large star with the most planets it can hold, say) can give systems
with worlds thousands of AU out; edit those by hand before using them.

## Contributing

Every pull request adds a release note under [`changes/`](changes/README.md)
instead of bumping the version; a workflow stamps the version after
merge. A change that alters what a seed generates also bumps the generator epoch
(`python -m tests.golden_update --reason "..." --affects fill`, [`changes/README.md`](changes/README.md)). Run the tests with `scripts/test.py changed` while you work and `scripts/test.py full --check` before
merging a finished piece of work; both need MySQL and Redis ([`docs/reference/testing.md`](docs/reference/testing.md)).
Plain `pytest -n auto` still works.

## License

[CC0 1.0 Universal](LICENSE.md).

The bundled list of common passwords (`src/planetgen/admin/common_passwords.txt.gz`)
is derived from [SecLists](https://github.com/danielmiessler/SecLists) under the
MIT licence; see `src/planetgen/admin/common_passwords.LICENSE`.

Where the data files in the repository come from, and what is known of their terms, is recorded in
[DATA-LICENSE.md](DATA-LICENSE.md).
