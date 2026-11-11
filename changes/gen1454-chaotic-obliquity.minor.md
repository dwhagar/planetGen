**Headline:** A planet in a crowded system whose tilt can wander chaotically under the pull of its neighbours now gets a tilt from its whole wander range, not the one it was made with.

### Changed

- GEN.1454 (part 2, a lock item): at the end of the settle pass each free-spinning planet of a star with other planets is checked against the spin-orbit resonances (`obliquity_variation`: the star's and moons' torque on its bulge against the nodal modes of its neighbours' pull on its orbit). A planet whose tilt lies in a chaotic zone, where two or more resonances overlap, takes a tilt drawn from that zone from a stream of its own; the axis keeps its angle round the orbit's normal, the sky of the other suns, the day and night shares and the planet's moons follow. Planets with a large moon, locked or resonant planets, lone planets and planets whose tilt is regular are not touched. In 100 sampled systems of four star types 80 planets changed tilt, almost all of them in crowded chains of 15 to 24 small planets whose slow precession puts many resonances inside their range. Registry section 6.271 (OB6). The resonance and chaos picture follows [R2571, R5220].

generation-output: changed
generation-affects: fill
