**Headline:** A neutron star or a black hole with a white dwarf companion no longer comes out as a pair that would already have merged by gravitational waves.

### Fixed
- GEN.1851: the Peters screen of `multistar._Builder` (a pair whose gravitational-wave merger time is shorter than the system's age is redrawn) applied to two neutron stars or black holes and to two white dwarfs, but a neutron star or a black hole with a white dwarf took the law of a compact object with a star and skipped it. In a seeded scatter of a million neutron star rows, 2 of 20,269 trees held such a pair; with the companion forced to a white dwarf and a 10 Gy system, 3 in 100 pairs with a neutron star and 6 in 100 with a black hole were tighter than a merged pair can be (a neutron star of 1.4 and a white dwarf of 1.03 solar masses at 0.015 AU, merger time 4.9 Gy against an age of 7.4 Gy). Every pair of two degenerate objects (white dwarf, neutron star, black hole, any mix) now meets the one rule, `degenerate_pair` replacing `double_white_dwarf`; the redraw kind `pair.gravitational_waves` is the same one. The multiplicity audit of TEST.162 checks every degenerate pair of every stored tree.

generation-output: changed
generation-affects: scatter, fill
