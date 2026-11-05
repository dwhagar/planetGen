**Headline:** The generator can now work out how bright any star would look from another system, with the dust of the galaxy in the way.

### Changed

- `physics/sky_photometry.py` (VIEW.20, the photometry of the system sky on arrays): `apparent_magnitude` and `absolute_visual_magnitude` for arrays of luminosity, temperature and distance; `bolometric_correction` for any temperature (Flower's fit from 4,100 K to 50,000 K, the main-sequence table of Pecaut and Mamajek under 4,100 K, which answers what VIEW.9 left open for a star cooler than 3,500 K, and the blackbody tail over 50,000 K); and `galactic_extinction_av`, the dust of the smooth disc of Drimmel and Spergel integrated along a line of sight. The Sun comes out at -26.75 at 1 au and Sirius A, Alpha Centauri A and Procyon A within 0.03 mag of their catalogue values. Nothing reads it yet: the sky catalogue (VIEW.23) and the glow band (VIEW.27) do. Four stand-ins are registered (6.171, SKP1 to SKP4). generation-output: unchanged
