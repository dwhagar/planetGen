**Headline:** a moon whose orbit the star's pull reshapes now sits on its new orbit, not on the one it had before.

### Changed

- GEN.1540: the Kozai-Lidov pass (GEN.1447) sets a cycling moon's eccentricity to its place in the cycle at the star's age, but the moon's position and speed stayed on the ellipse it was drawn on, so a captured moon pumped down from a high eccentricity could sit far outside `a (1 + e)` (seen once as 0.507 AU against 0.233 to 0.313 AU; seed 1 of a G2V scan shows it every time). The moon is placed again on its new ellipse, in its place in its plane, as the tide settling does. Generation changes for the moons the pass changes: their stored position and speed.
