**Headline:** A star or system told to have no binaries (or exactly one star) no longer gets a white dwarf or other compact companion anyway.

### Changed
- **A forced single object stays single (GEN.719, bug).** `draw_member_tree` with `n_stars=1` (what a binary prevalence of -100 or a forced single star asks for) drew the compact companion of a living star, or the companion of a placed white dwarf, neutron star or black hole, from the companion table regardless; an old star that had become a white dwarf therefore still came out as a wide double white dwarf. One star is now one member for every placed object, with no companion of any kind. Systems built without a directive draw exactly as before.
- **Tests.** A forced single object has one member and no orbit for a living star, a star that is now a white dwarf, and a placed white dwarf, neutron star and black hole, over 150 seeds each.
