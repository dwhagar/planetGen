**Headline:** Rocky planets and moons now carry a crust thickness, a craton share, a crustal sulfur index, and a record of which magmatic ore types (nickel-copper-PGE sulfide, chromitite, iron-titanium-vanadium oxide) they hold, how much and how rich.

### Added
- **Crust fields (GEN.1380, `physics/ore_deposits.py`).** `crust_thickness_km` (thick on a small stagnant body, thin on a hot mantle or a mobile lid), `craton_fraction` (a mobile lid, 2 Gyr or older, with granitoid ground) and `crust_sulfur_index` (the star's S/Fe, liquid water and age), on `planets`, `moons` and `rogue_planets`.
- **Magmatic ore types.** The new child tables `planet_deposit_types`, `moon_deposit_types` and `rogue_planet_deposit_types` hold one row per type a rocky body was assessed for: `status` (present, absent by state or absent by chance), the ore mass, the grade of the key element capped by what the host rock supplies, the age window and a capped flag. The gate is the magma volume (the mantle processed through melting), the sulfur index for the sulfide, and the mantle redox. A metal-poor star gives leaner ore. Stand-ins ORE1 to ORE7 (section 6.215). Later slices (GEN.1381 to GEN.1383, GEN.1378) add the other types to the same tables.
- Schema v228 (migration 0228; the train renumbers it).

### Changed
- Every rocky planet and moon made again has the new fields and rows, so the next epoch is needed (the stored bodies move; no existing draw of the run's stream moves, the ore uses the body's own streams).

generation-output: changed
TODO items resolved: GEN.1380
