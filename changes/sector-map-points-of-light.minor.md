### Changed
- **Stars and glowing phenomena are points of light on the Sector Map
  (MAP.15).** Every star is now drawn the way the Galaxy Map draws its
  bright stars: a tiny bright core in a soft aura a fixed number of
  pixels across, the core sized by the star's radius, the aura's width
  and brightness by its luminosity, and the color by its temperature,
  instead of a textured ball with a glow shell. Quasars, neutron stars
  and accreting black holes are points of light too; quiescent black
  holes, rogue planets (still ringed by "Mark rogue planets") and
  interstellar comets keep their spheres, and nebulae, supernova
  remnants and asteroid fields stay clouds. Points grow a little (up to
  1.5 times) as the camera closes in, a click within a few pixels of one
  still selects it and shows its details, and on the light theme a
  point's core is a darker shade of its color so it stays visible.
  `lib/starmap.py` now sends each star's point (`light`: color, core
  and aura size in pixels, aura strength, core brightness) in place of
  the old glow-shell numbers.
