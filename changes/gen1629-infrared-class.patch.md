**Headline:** a gas the gas table gave no infrared flag (52 of 124: the atoms, radicals and rock vapours) is now treated one way by the greenhouse table and the line opacity, by its structure.

### Fixed
- GEN.1629: 52 gases had `ir` null and three readers treated the null three ways (the greenhouse table's features counted it as no infrared share, its thin-air scaling counted it as absorbing, the line opacity gave it no band). `Gas.infrared_class` returns the row's flag, or for a row with none the gas table's own rule (monatomic none, homonuclear diatomics collision-induced only, the rest infrared): an atom `none`, a diatomic of one element `cia`, any other molecule `ir`, which takes the generic stand-in band. The table, the thin-air scaling and the line opacity read it (registry row IR21, tier 7, tweak L0). OH, FeH, VO, MgO, SiO2, TiO2, NaOH and HCOOH (and 40 more) now count as infrared-active where they are in an air; Cl, H, O and the other atoms do not.

generation-output: changed (lock item)
generation-affects: fill
