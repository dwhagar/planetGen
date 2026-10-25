**Headline:** The impact history, the thermosphere and the planet-survey modules take G, the Sun's mass and Earth's mass and radius from the shared constants.

### Changed
- **Typed copies of G, solar mass, Earth mass and Earth radius are gone (GEN.587, a further file group).** Impact history, the
  thermosphere, the Kepler recovery and the radius-valley survey read `physics/constants.py` instead of their own rounded copies
  (G 6.674e-11 becomes 6.6743e-11, the Sun 1.989e30 kg becomes 1.98841e30 kg, Earth 5.972e24 kg becomes 5.97217e24 kg). The last
  digits of the impact-history escape speeds move, so this rides the generator epoch of the moon and rogue ring change.

generation-output: changed
