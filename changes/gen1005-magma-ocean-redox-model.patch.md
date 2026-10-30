**Headline:** The oxygen fugacity of a magma ocean now has a model of its own (the pressure at its base, the ferric share of its iron), ready for the rows that draw it.

### Changed

- GEN.1005 (part 1, lock-free): `physics/magma_ocean_redox.py` gives the magma-ocean surface redox of a rocky body from its mass through the pressure at the base of its ocean, the ferric share of the iron and the fourth-power relation of the ratio to the oxygen fugacity, and the later upper mantle's value; the Moon, Mars and the Earth are tests of it and no longer inputs. The recipe expressions can call it (`magma_ocean_diw`, `magma_ocean_diw_sigma`, `later_mantle_diw`). Registry rows MR1 to MR7 (section 6.57 of the gap-filling plan) carry tiers and tweaks. Rows a1 and a2 still draw the old line until part 2, so nothing the generator draws changes.
