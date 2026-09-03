# planetGen

A procedural planet and star system generator, designed for the Molten Aether FFRP game. The output is designed to be easily copied and pasted into the wiki.

## Features

*   **Star generation**: Stars are generated with physically-grounded mass, radius, temperature, and luminosity, either fully at random (weighted by realistic galactic spectral-class prevalence) or pinned to a specific spectral type and Yerkes luminosity class (e.g. `G2V`) via `--star-type`. `--force-large-star` biases generation toward hotter, more massive stars.
*   **Binary star systems**: `--binary-system` generates a P-type (circumbinary) binary pair — a primary and secondary star orbiting each other — represented as a single unified effective star (combined mass, luminosity, and habitable zone) for the purposes of planet placement, while still reporting each star's individual properties in the output.
*   **Planet and moon generation**: Planets are drawn from 25 distinct planet classes (terrestrial and gas giant), each with its own composition, atmosphere, and valid orbital zones (hot/ecosphere/cold). Planets can generate their own moons, with orbital placement, atmospheric conditions, and surface gravity calculated per body. Orbital spacing is validated against each object's Hill sphere to keep the system physically plausible.
*   **Asteroid belts**: Belts can appear between planets (or be forced with `--force-asteroid-belt`), each with a randomly generated density and mineral/gem composition.
*   **Stellar age and lifespan modeling**: Star age and lifespan are derived from spectral and Yerkes class data, then adjusted (via `--age young|old`) or extended as needed so any habitable planets' life stages remain consistent with how long the star has existed and how much longer it has left.
*   **Life chemistry and evolutionary timelines**: Habitable planets are evaluated against the star's spectral class to determine which of several life chemistries (e.g. Chlorophyll a, Melanin, Retinal) are viable, each with its own evolutionary pace (fast/normal/slow). A speculative evolutionary timeline (abiogenesis through technological civilization) is generated for habitable worlds, influenced by `--force-intelligent-life` / `--no-intelligent-life`.
*   **Flavor text**: Randomly-selected descriptive "sensor readings" flavor text can be appended to systems and planets, with limits and chances controllable via `--flavor-chance-system`, `--flavor-chance-planet`, and `--max-planet-flavor`.
*   **Dual output formatting**: Every generated system can be rendered as either MediaWiki wikitext templates (default, ready to paste into the wiki) or Markdown (`--markdown`).

## Setup

This project uses the `nltk` library to generate phonetically pleasing names. To install the necessary dependencies, you will need to have `setuptools` installed.

First, ensure you have `setuptools` installed:

```bash
pip install setuptools
```

Then, run the setup script to install the project dependencies:

```bash
python setup.py install
```

This will install the `nltk` library and download the 'words' corpus, which is required for the name generation.

## Usage

To generate a new star system, run the `planetGen.py` script from the root of the project:

```bash
python planetGen.py [options]
```

### Options

*   `--force-habitable-world`, `-fhw`: Force the generation of a habitable world in the system.
*   `--force-asteroid-belt`, `-fab`: Force the generation of an asteroid belt in the system.
*   `--force-large-star`, `-fls`: Force the generation of a large star.
*   `--force-moons`, `-fm`: Force the generation of lots of moons in a system.
*   `--force-max-planets`, `-fmp`: Force the system to maximized the number of planets.
*   `--absurd`: Force the system to generate the largest star possible with the maximum number of planets and moons.
*   `--output`, `-o <file_path>`: Specifies a file path to write the output to.
*   `--markdown`, `-m`: Formats the output in Markdown instead of the default wikitext.
*   `--no-planets`, `-np`: Skip planet generation and create only a star.
*   `--star-type <type>`: Force the generation of a specific star type (e.g., G2V).
*   `--name <name>`: Specifies a name for the star system, overriding the default random generation.
*   `--age <young|old>`: Specifies the age of the star system (young or old).
*   `--force-intelligent-life`, `-fil`: Ensures at least one planet with intelligent life is generated. Implies `--force-habitable-world`.
*   `--no-intelligent-life`, `-nil`: Ensures no planet with intelligent life is generated. Implies `--force-habitable-world`.
*   `--no-habitable-world`, `-nhw`: Ensures no habitable world is generated in the system.
*   `--binary-system`, `-bs`: Generates a binary star system (P-type).
*   `--flavor-chance-system <float>`: Overrides the default system-level flavor text chance (0.0 to 1.0).
*   `--flavor-chance-planet <float>`: Overrides the default planet-level flavor text chance (0.0 to 1.0).
*   `--max-planet-flavor`: Sets the maximum flavor text total for planets to 99.

**Note on Incompatible Options:**

*   `--no-planets` cannot be combined with `--force-moons`, `--force-max-planets`, `--absurd`, or `--force-habitable-world`.
*   `--star-type` cannot be combined with `--force-large-star`.
*   `--force-intelligent-life` and `--no-intelligent-life` cannot be used together.
*   `--force-habitable-world` and `--no-habitable-world` cannot be used together.
*   `--flavor-chance-system` must be a float between 0.0 and 1.0.
*   `--flavor-chance-planet` must be a float between 0.0 and 1.0.

### Additional Information

This tool is designed as a personal tool for the Molten Aether FFRP game. The output is designed to be simply cut and paste from the program output into the wiki. See https://wiki.moltenaether.com for wiki and game information.

Using commands to force a habitable world and an asteroid belt will automatically force a large star to ensure there is room for both objects. The most common stars are small dwarf stars which make a smaller star system. Forcing a large star as well as maximized planets will cause the generated system to be very large with an extremely high and sometimes absurd number of planets. Do not assume that just because it is generated here, it is accurate or possible, such as large systems may require editing as some worlds may end up saying they are several thousand AU's from the central star.

## License

This project is licensed under the [CC0 1.0 Universal](LICENSE.md) license.