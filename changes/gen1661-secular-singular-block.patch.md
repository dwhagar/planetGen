**Headline:** a system with two close-in planets and a distant belt no longer fails to generate in the secular exchange ("Singular matrix", about 4 percent of systems asked for a habitable world and a belt).

### Changed

- `physics/secular_chain.py` (GEN.1661): the fixed point the held bodies set for the free ones (`_drive`) is zero when the held bodies are at rest (the flat rings of a belt), and the shortest solution when the free planets' block has a zero mode (the invariable plane) that the held bodies barely touch. Before, `np.linalg.solve` was called on the nodal block of the planets, whose rows sum to zero but for a belt's pull of 1e-17 of their own, and raised. A system whose block was not singular is unchanged.
- Census: the 600-system world-requested redraw group had 23 crashes (seeds 40, 53, 83, 97, 124 among them) and now has none.

generation-output: unchanged
