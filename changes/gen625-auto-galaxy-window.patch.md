**Headline:** A default galaxy no longer raises the plan's coverage warning, and the average default galaxy is about half the size of a Milky Way disc.

### Changed
- **The all-Auto size and flatness are redrawn until the plan can hold the shape (GEN.625).** When the size, flatness, bulge share and bulge size are all Auto, `shape_draw._auto_size` takes the seed's draw of the disc size (Shen et al. 2003) and flatness (Kregel et al. 2002) and, when the shape would raise a grid warning (the thin disc's crest under 20 times the one-star density, a scale height under 8 pc, a bulge under 4 pc) or the flatness would leave the Structure dial's 3 to 20, draws the next pair from the same two streams, up to 32 times. A seed whose first draw is fine keeps its galaxy; about one in five gets a later draw. Over 3,000 seeds the share of default plans that warn goes from 18 percent to none, and the mean disc scale length from 1,515 pc to 1,235 pc (47.5 percent of a 2.6 kpc Milky Way disc, a test yardstick only). Pinned dials keep the seed's own first draws, and the warnings still fire for them. The Default (Auto) column of the dial table (`docs/cli.md`) says so.
- `tests/test_shape_draw.py` draws 1,000 all-Auto galaxies and finds no grid warning, a mean size and mass at most half a Milky Way's, and every draw inside the dial ranges; the law tests read size and flatness draws with a neighbour pinned, which leaves the window out of the way.
- Generator epoch 51 (numbered at merge): the default size and flatness of about one seed in five change.

generation-output: changed

TODO items resolved: GEN.625. Added: none.
