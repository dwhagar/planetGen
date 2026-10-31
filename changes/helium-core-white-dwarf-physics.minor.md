**Headline:** A white dwarf in a close pair can now be a helium core of 0.18 to 0.45 solar masses, the one its giant had when it filled its lobe, instead of always 0.5 solar masses or more.

### Changed
- **The helium-core white dwarf of a close binary (GEN.996).** A star under 2.3 solar masses that fills its Roche lobe on the red giant branch leaves the helium core it then had, 0.18 to 0.45 solar masses, and not the 0.55 or more the single-star relation gives. `physics/helium_core_white_dwarf.py` holds Rappaport et al. 1995's core mass and radius relation (read second hand), the helium flash at a core of 0.45 solar masses and the lobe of Eggleton 1983. The tree builder (`multistar._Builder.helium_core`) now reads it: a companion white dwarf in a post-common-envelope pair takes the core its giant had when its lobe reached the pair's periastron, so the pair is lighter, its period longer and the node masses above it follow. No rate is set. In pairs of a solar-type star and a white dwarf companion at 9 gigayears about 9 percent of the white dwarfs come out helium-core, from the generator's own period law and mass-ratio law (the 20 pc sample has 7.4 percent of all white dwarfs).
- **The stored member needs no flag.** A white dwarf lighter than 0.45 solar masses can only be a helium core (the single-star relation gives none under 0.5), so the rebuild of a stored tree makes it with the stored mass (`stored_trees.helium_core_mass_sol`, `member_record.make(white_dwarf_mass_sol=)`, `evolve_star(white_dwarf_mass_sol=)`). No migration.
- **Left for a follow-up.** The placed white dwarf of a tree (member 0) is built before its tree and stays as its own star draw made it; only companions are helium cores. The white dwarf cools from the end of its progenitor's whole life, as every white dwarf does, though its envelope left earlier (stand-in HC5, section 6.160).

generation-output: changed

TODO items advanced: GEN.996 (wired, lock window).
