# planetGen

A procedural planet and star system generator, designed for the Molten Aether FFRP game. The output is designed to be easily copied and pasted into the wiki.

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
*   `--age <young|old>`: Specifies the age of the star system, either "young" or "old".
*   `--force-intelligent-life`, `-fil`: Ensures at least one planet with intelligent life is generated. Implies `--force-habitable-world`.
*   `--no-intelligent-life`, `-nil`: Ensures no planet with intelligent life is generated. Implies `--force-habitable-world`.
*   `--no-habitable-world`, `-nhw`: Ensures no habitable world is generated in the system.

**Note on Incompatible Options:**

*   `--no-planets` cannot be combined with `--force-moons`, `--force-max-planets`, `--absurd`, or `--force-habitable-world`.
*   `--star-type` cannot be combined with `--force-large-star`.
*   `--force-intelligent-life` and `--no-intelligent-life` cannot be used together.
*   `--force-habitable-world` and `--no-habitable-world` cannot be used together.

### Additional Information

This tool is designed as a personal tool for the Molten Aether FFRP game. The output is designed to be simply cut and paste from the program output into the wiki. See https://wiki.moltenaether.com for wiki and game information.

Using commands to force a habitable world and an asteroid belt will automatically force a large star to ensure there is room for both objects. The most common stars are small dwarf stars which make a smaller star system. Forcing a large star as well as maximized planets will cause the generated system to be very large with an extremely high and sometimes absurd number of planets. Do not assume that just because it is generated here, it is accurate or possible, such as large systems may require editing as some worlds may end up saying they are several thousand AU's from the central star.

## License

This project is licensed under the [CC0 1.0 Universal](LICENSE.md) license.