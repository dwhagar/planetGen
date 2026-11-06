**Headline:** A system asked to have no habitable world now has no planet or moon of a life-carrying class anywhere, where about a third of them still held one before.

### Changed
- **The habitable-world bar covers the whole system (GEN.1275).** `HABITABLE_WORLD` set to False barred only an ecosphere planet that met the request, so about 0.35 of such systems still held a body of a class that can carry life: a moon, an outer planet, or a body whose class changed after the planet was drawn. The bar of GEN.1250 (`StarSystem._bar_habitable_worlds`, which judges the finished bodies and makes a barred one's end stages again, up to 30 attempts, each from a stream of its own so the run's stream does not move) now keeps out a body that meets the request's conditions or whose finished class can carry life. A planet that comes out in a life-carrying class at its first draw is drawn again too. A body that never leaves those classes is dropped from the system, its orbit and moons with it. In 60 systems drawn this way none holds a body of those classes, and none reports the bar as unmet. A barred system costs about 12 seconds to draw against about 3, so the habitable-world prevalence cases draw 50 systems.

generation-output: changed

TODO items resolved: GEN.1275.
