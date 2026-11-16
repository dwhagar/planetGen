**Headline:** A carbon-rock world now holds the phosphorus of its star, so a living sea on one has phosphate to live on.

### Fixed

- GEN.1749 (every galaxy): `bulk_draw._partition` gave a carbon rock no minor oxide at all, so `rock_p2o5_fraction` was NULL, `phosphorus_budget.generate` returned early, the ocean's phosphate was None and the life law's chemistry term was 0 on every carbon-rock world (Oleodzo II, seed 30391, a sea world in light, was one). The carbon rock now keeps the star's phosphorus as P2O5 in its rock and the share the core takes in its metal, by `mantle_oxides.oxides` as a silicate rock does, and holds none of the other six oxides (they stay NULL, so the norm and ground types, which read no carbon rock, are as they were).

generation-output: changed
generation-affects: fill
