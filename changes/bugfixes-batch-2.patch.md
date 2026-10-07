### Fixed
- **A black hole without a disk shows its Hawking temperature and
  luminosity** (GEN.82) instead of zero, from its mass.
- **The sector summary names star kinds plainly and logs one line per
  entry** (UX.34, OPS.9): white dwarfs, neutron stars, black holes,
  giants and supergiants instead of raw spectral codes.
- **Species are kept only for worlds with a technological civilization**
  (GEN.80); a population pass removes any stored species without one.
- **Changing a planet's or moon's class re-generates it as that class**
  (ADM.27): radius, density, composition, atmosphere, temperature,
  pressure and life, keeping its orbit, mass and name. A class its mass
  can't fit (a gas giant class for a small rocky world) is refused with
  a message and nothing changes.
