# Decisions, 2026-09

20 sentences in the documents (outside TODO.md and the number map) that name Boss and a date in 2026-09; the choices as they were made, oldest first, with the file and line they are written on.

- 2026-09-24 ([docs/design/galaxy/spherical-shell-sectors.md:102](../design/galaxy/spherical-shell-sectors.md)): The no-backward-compatibility rule makes this acceptable; it is the same kind of reset Boss accepted on 2026-09-24, on a codebase many times larger.
- 2026-09-30 ([docs/archive/galaxy-coordinate-system-shells.md:49](../archive/galaxy-coordinate-system-shells.md)): master-wedge rule (v35, Boss's choice of 2026-09-30):
- 2026-09-30 ([docs/design/architecture/design-decisions.md:46](../design/architecture/design-decisions.md)): Boss chose this hybrid on 2026-09-30.
- 2026-09-30 ([docs/design/architecture/designation-history.md:19](../design/architecture/designation-history.md)): - **Names change today, with no memory.** `store.rename_star_system`, `rename_star` and `rename_body` (reached by `PATCH /api/systems\|stars\|planets\|moons/&lt;id>`) overwrite `name`; planets and moons named after their system are renamed with it (`_rename_bodies_with_prefix`). `store.name_in_use` keeps sector, system, star and the uniquely named phenomena unique (`UNIQUE_NAME_TABLES`); a planet or moon renamed by hand ma
- 2026-09-30 ([docs/design/architecture/units-and-number-formatting.md:44](../design/architecture/units-and-number-formatting.md)): - Distances (Boss, 2026-09-30): a parsec-family value carries one familiar unit in parentheses, "4.2 pc (13.7 ly)", "2.4 mpc (495 AU)" (`DISTANCE_PAREN_MIN_LY`, `DISTANCE_PAREN_MIN_AU` in `util/format.py`).
- 2026-09-30 ([docs/design/architecture/units-and-number-formatting.md:45](../design/architecture/units-and-number-formatting.md)): - Body radii (Boss, 2026-09-30): a planet, moon or star radius is always km in scientific notation, the one exception to the distance ladder (`format_body_radius_km`).
- 2026-09-30 ([docs/design/architecture/units-and-number-formatting.md:92](../design/architecture/units-and-number-formatting.md)): \| Body radius \| m \| Stays km in scientific notation with "(x Earth radii)" (Boss, 2026-09-30) \| Earth radii for planets; add Jupiter radii for gas giants \| `format_body_radius_km` is a recorded decision.
- 2026-09-30 ([docs/design/galaxy/multiple-galaxies.md:22](../design/galaxy/multiple-galaxies.md)): - Navigation frames are Boss's design of 2026-09-30 (navigation-frames.md): "North" points toward the local dominant center of mass.
- 2026-09-30 ([docs/design/orbits/course-routing.md:661](../design/orbits/course-routing.md)): 15, 21, 27, 36, ...; Boss, 2026-09-30); lookups use `sector_address_at` and
- 2026-09-30 ([docs/design/orbits/navigation-frames.md:5](../design/orbits/navigation-frames.md)): Boss's design for nested navigation frames, recorded 2026-09-30 (NAV.1
- 2026-09-30 ([docs/design/orbits/navigation-frames.md:26](../design/orbits/navigation-frames.md)): Decisions Boss approved on 2026-09-30:
- 2026-09-30 ([docs/design/orbits/navigation-frames.md:275](../design/orbits/navigation-frames.md)): - **Mark as elevation mod 360.** Approved by Boss on 2026-09-30 (the list at
- 2026-09-30 ([docs/design/orbits/sky-view.md:27](../design/orbits/sky-view.md)): - Boss, 2026-09-30 (navigation-frames.md): east is North x Up.
- 2026-09-30 ([docs/design/systems/interstellar-object-rates.md:6](../design/systems/interstellar-object-rates.md)): from the research Boss supplied on 2026-09-30.
- 2026-09-30 ([docs/design/systems/nebula-and-asteroid-field-classes.md:3](../design/systems/nebula-and-asteroid-field-classes.md)): The design for `docs/TODO.md` items GEN.10 to GEN.14, recorded 2026-09-30 from Boss's
- 2026-09-30 ([docs/design/systems/nebula-and-asteroid-field-classes.md:124](../design/systems/nebula-and-asteroid-field-classes.md)): Boss (2026-09-30): "stellar remnants" means supernova remnants only; the
- 2026-09-30 ([docs/design/systems/nebula-and-asteroid-field-classes.md:210](../design/systems/nebula-and-asteroid-field-classes.md)): The 5e-6 is Boss's chosen research value (2026-09-30), so this stays a
- 2026-09-30 ([docs/design/systems/nebula-and-asteroid-field-classes.md:256](../design/systems/nebula-and-asteroid-field-classes.md)): Boss (2026-09-30): size is "a digit and part of the class", so a full
- 2026-09-30 ([docs/design/systems/nebula-and-asteroid-field-classes.md:270](../design/systems/nebula-and-asteroid-field-classes.md)): - **Supernova remnants only.** Boss said on 2026-09-30 that "stellar
- 2026-09-30 ([docs/design/systems/nebula-and-asteroid-field-classes.md:292](../design/systems/nebula-and-asteroid-field-classes.md)): - Letters for black holes and neutron stars: not wanted (Boss, 2026-09-30).
