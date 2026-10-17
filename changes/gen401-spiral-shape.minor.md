**Headline:** A new galaxy now draws its size, thickness, bulge, bar and spiral arms from its seed instead of every galaxy being built to the Milky Way's measurements.

### Changed
- **The spiral shape is drawn (GEN.401).** `planetgen plan` and the Generate page leave `--disk-scale-length-pc`, `--disk-scale-height-pc`, `--bulge-scale-radius-pc`, `--bulge-amplitude`, `--arm-count` and `--pitch-angle-deg` blank by default; blank means drawn from the galaxy seed (`galaxy/shape_draw.py`). The disk length follows the size-mass relation of late-type galaxies from the seed's stellar mass (GEN.400), the height follows the length, the bulge's share and size follow the galaxy's stage, and the arm count and pitch angle follow the arm statistics. An option given still replaces its draw, and the bulge's density is solved for the disk used. The plan prints what it drew.
- The old defaults (2,600 pc, 300 pc, 1,580 pc, 3.11, two arms at 15 degrees) are gone from the plan, the CLI and the page. They stay only as the reference shape the math check measures the model against.
- Bar presence and the bar's axes and angle stay the model's constants until the shape columns of DB.66; the stage rule is a stand-in (design note 5.13.8).

TODO items resolved: GEN.401. Added: none.
