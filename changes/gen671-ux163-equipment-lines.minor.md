**Headline:** Cosmic-ray worlds like Mars and the Moon now require a hardened suit, and every planet and moon page shows a Required line, a Recommended line (each left off when empty) and the per-item specs.

### Changed
- **A hardened suit against cosmic rays (GEN.671).** Where no ordinary garment or vest brings the yearly dose under the Green edge, the hardened suit is required and the hardened habitat is recommended. Mars is now "Required: hardened pressure suit + dust-tight seals. Recommended: hardened habitat." Run `planetgen rescore-habitability` to update stored worlds.
- **A Required line and a Recommended line (UX.163).** A planet or moon's details head the habitability rating with a Required line and a Recommended line (a line with nothing in it is left off), followed by a Survival equipment table with each item's severity, category and tech level, and the shelter's shielding and hours outside. The API's `equipment` object gains `required_line` and `recommended_line`.
