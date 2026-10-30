**Headline:** Every rocky or icy planet and moon now has a porosity: the mean pore space of its upper crust, higher where the gravity is lower.

### Changed
- **GEN.780: porosity is filled for every planet and moon that has a crust.** A body of 1e21 kg or more (a planet, or a moon that large) has no macroporosity, so it now stores the mean porosity of its upper 10 km (`physics/crust_porosity`): the pores close at a fixed pressure, so the depth that stays porous is the Earth's scale depth (2.82 km) times the Earth's gravity over its own, and the surface value is 20 percent for rock and 30 percent for ice, by the share of water and ice. It comes out near 5.5 percent for the Earth, 11 for Mars and 15 for the Moon, and higher for a small, low-gravity or icy body. A moon lighter than that, which is a heap of grains, keeps its macroporosity from the compaction law; one with no drawn value now takes the tapered law of Bierson and Nimmo instead of an empty cell. A gas giant has no crust and none. The sheet rule no longer excuses the column, and `planets.porosity` is written (TEST.461's exclusion is removed).

generation-output: changed
