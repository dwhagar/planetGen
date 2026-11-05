**Headline:** A system asked to have no habitable world now has no planet or moon of a life-carrying class anywhere, where about a third of them still held one before.

### Changed
- **The habitable-world bar covers the whole system (GEN.1275).** `HABITABLE_WORLD` set to False barred only an ecosphere planet that met the request, so about 0.35 of such systems still held a body of a class that can carry life: a moon, an outer planet, or a body whose class changed after the planet was drawn. The bar now looks at the finished class of every planet and moon. A planet that comes out in a life-carrying class is drawn again, and a body whose class is settled later in the life pass is remade from the same conditions (up to the same 30 attempts, from a stream of its own so the run's stream does not move). A body that never leaves those classes is dropped from the system, its orbit and moons with it. In 60 systems drawn this way none holds a body of those classes, and none reports the bar as unmet.

generation-output: changed

TODO items resolved: GEN.1275.
