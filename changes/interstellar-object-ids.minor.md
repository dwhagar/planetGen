### Changed

- Every interstellar object in a galaxy-placed sector is now named by a
  76-bit position ID, shown as 19 hex digits, instead of a generated name
  or designation (GEN.64). The ID packs the object's type, distance from
  the galactic center (with its unit, mpc to Gpc), bearing, mark and a
  collision number that tells apart up to sixteen objects at the same spot.
  This covers rogue planets, standalone black holes and neutron stars,
  supernova remnants and their collapsed cores, nebulae, quasars,
  interstellar comets, asteroid fields, and star systems built around a
  bright-sweep star, whose stars, planets and moons are named from the
  ID. Ordinary star systems keep their names, and names given with
  `--name` are kept. These objects no longer go through the name
  registry, which made dense sectors about 1.8 times faster to generate
  with one worker and 2.3 times with two. Objects saved before this keep
  their names. See `docs/design/object-ids.md`.
