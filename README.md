# planetGen

**Version:** 7.85.177 &middot; [Changelog](CHANGELOG.md) &middot; [Repository](https://github.com/dwhagar/planetGen) &middot; License: [CC0 1.0 Universal](LICENSE.md)

planetGen generates a galaxy: stars, star systems, planets, moons,
asteroid belts and exotic phenomena, placed in a physically modeled
spiral galaxy and saved to a MySQL database. A website on top of it lets
you browse and search the galaxy, fly around it on 3D maps, and copy any
system's write-up as MediaWiki wikitext or Markdown. It was built for
the Molten Aether FFRP game, so a generated system can be pasted straight
into the [game's wiki](https://wiki.moltenaether.com).

**To install it, see [INSTALL.md](INSTALL.md).**

## What it does

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
  cylindrical grid. Planning the galaxy also places every bright star
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

- Python 3.9 or later (3.10+ on macOS).
- MySQL 8.0.16+ or MariaDB 10.4+.
- For the website, a web server: Apache2 with mod_wsgi on Debian or
  Ubuntu (the reference setup), nginx or Caddy on Linux, IIS, Caddy or
  Apache on Windows, or nginx on macOS.
- Install and update scripts for all three: `install.sh` and
  `update.sh` on Linux and macOS, `install.ps1` and `update.ps1` on
  Windows.

[INSTALL.md](INSTALL.md) lists the full requirements and walks through
the install.

## Usage

### The website

Once installed, the site's front page leads to:

- **Galaxy Map**: the whole galaxy in 3D, with every bright star as a
  tiny point in a big soft glow and nebulae as translucent clouds. It
  opens on a drill-down: pick a slab of the galaxy, then a block, then
  smaller blocks down to single sectors, with a breadcrumb, keyboard
  keys, and a link for every stage so Back and Forward work. Type a
  sector, a `ring/layer/slot` address, coordinates or a name to fly
  there. A NAV course can be drawn on it, and once population has run,
  a Territories overlay shows who holds what. As an admin, click an
  empty sector to generate it, a neighborhood of a radius you choose,
  its column or its whole shell, or generate a whole small block or one
  of its layers.
- **Sectors** and **Systems**: every generated sector and every system,
  50 to a page. A sector has a 3D Sector Map (with the same Generate
  buttons for admins, and Nav from/to on any system) and a Contents
  list; a system has its stars, planets, moons and facilities, the cloud
  it sits in, its nearest systems, a System Map, and its wiki write-up
  ready to copy. Phenomena have pages too, with their class, what fills
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

[`docs/html-interface.md`](docs/html-interface.md) describes every page
and [`docs/api.md`](docs/api.md) the JSON API.

### The command line

`generate.py`, in the checkout, is the one entry point for every
generator. On Linux and macOS the installer also adds a `planetgen`
command that runs it; on Windows run it with the venv's Python
(`C:\srv\planetgen-venv\Scripts\python.exe generate.py`). Each
subcommand saves what it makes to the database unless told otherwise:

```bash
python3 generate.py plan         # plan the galaxy's shape (once, before galaxy)
python3 generate.py galaxy       # generate sectors in the galaxy
python3 generate.py sector       # one or more standalone sectors
python3 generate.py system       # one star system
python3 generate.py phenomenon   # one exotic phenomenon
python3 generate.py population   # optional: species, civilizations and territories
```

Common examples:

```bash
# A random neighborhood: a random start and every sector within 100 ly of it
python3 generate.py galaxy

# One whole ring of the galactic plane, one sector of it, or one column
# of sectors through every layer
python3 generate.py galaxy --ring 12
python3 generate.py galaxy --ring 12 --layer 0 --slot 5
python3 generate.py galaxy --ring 12 --slot 5 --column

# A quick test galaxy: plan with fewer pre-placed bright stars
python3 generate.py plan --bright-star-min-luminosity 5000

# A single system with a habitable world, printed as Markdown, not saved
python3 generate.py system +habitable_world --markdown --output -

# A system built from a specification file
python3 generate.py system --system-file examples/systems/solar_system.json
```

Most options use a `+name`/`-name` form: `+moons` forces moons,
`-moons` forbids them, and leaving it off leaves it to chance.
`python3 generate.py <command> --help` lists every option, and
[`docs/cli.md`](docs/cli.md) is the full reference, including the
system specification file format and the incompatible option pairs.

### Other tools

All in `src/`; `python3 src/<tool>.py --help` shows each one's options:

- **`updateOrbits.py`** advances the galaxy by the time elapsed since
  the last run: planets, moons and orbiting facilities around their
  hosts, and every system and phenomenon along its galactic orbit. Run
  it monthly; the maintenance timers in
  [INSTALL.md](INSTALL.md#scheduled-maintenance) do.
- **`migrateDb.py`** brings the database schema up to date. The install
  and update scripts run it for you.
- **`resetDb.py`** deletes all generated galaxy data (asks for the
  database name first) so you can start a new galaxy.
- **`queryDb.py`** searches the database from the terminal.

## Documentation

| Document | What it covers |
|---|---|
| [INSTALL.md](INSTALL.md) | Installing, updating and the example configs |
| [docs/cli.md](docs/cli.md) | Every command-line option |
| [docs/config.md](docs/config.md) | Every `config.json` setting |
| [docs/deployment/](docs/deployment/README.md) | Hosting guides per platform |
| [docs/html-interface.md](docs/html-interface.md) | How the website works |
| [docs/api.md](docs/api.md) | The JSON API |
| [docs/database-schema.md](docs/database-schema.md) | The database tables and migrations |
| [docs/design/architecture.md](docs/design/architecture.md) | How the program fits together: every file and the main flows |
| [docs/design/](docs/design/) | Design notes: the galaxy model, coordinates, navigation, phenomena, and why each choice was made |
| [docs/testing.md](docs/testing.md) | Running the tests |
| [CHANGELOG.md](CHANGELOG.md) | What changed in each release |

## A note on accuracy

planetGen aims for plausibility, not a simulation. Forcing extremes
(a large star with the most planets it can hold, say) can give systems
with worlds thousands of AU out; edit those by hand before using them.

## Contributing

Every pull request adds a release note under [`changes/`](changes/README.md)
instead of bumping the version; a workflow stamps the version after
merge. Run the tests with `pytest` ([`docs/testing.md`](docs/testing.md)).

## License

[CC0 1.0 Universal](LICENSE.md).

The bundled list of common passwords (`src/stellarObjects/common_passwords.txt.gz`)
is derived from [SecLists](https://github.com/danielmiessler/SecLists) under the
MIT licence; see `src/stellarObjects/common_passwords.LICENSE`. The QR code
generator (`src/stellarObjects/qrcodegen.py`) is Project Nayuki's, under
the MIT licence in its header.
