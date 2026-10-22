**Headline:** The galaxy density model can now be read as a list of components, each with an amplitude and a mass.

### Added
- `galaxy.density.density_components(shape)` lists the thin disk, thick disk, bulge, halo and the two nuclear parts, each with its population, amplitude, mass (in closed form and in solar masses) and the lengths, heights and slopes that shape it. A component with amplitude 0 has mass exactly 0. The density itself does not change, so every galaxy generates as before; the galactic potential (GEN.338) and the Sersic spheroid (GEN.381) read and extend the list.
