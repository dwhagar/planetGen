**Headline:** Worlds with life now hold the gases their kind of life makes at their stage: methane and hydrogen drawdown for the hydrogen, donor and anoxygenic classes, a methane haze that follows it, and no methane once the oxygen has risen.

### Changed
- **The class gas ledger (GEN.891, GEN.857, GEN.858, GEN.859).** `physics/data/life_ledger.json` holds one row per life class and `physics/life_air.py` applies it after the oxygen. A built class (1, 2, 3, 6 and 8) that has producers by the star's age puts methane into an air with no free oxygen: its flux (the carbon the oxygenic class buries, or the donor and hydrogen supplies of the other classes over Earth's) times the 30,000 year life of methane, at most 3.5 percent of the air's molecules, with the carbon taken from the air's own carbon dioxide (a tenth of the carbon stays). The hydrogen eaters draw hydrogen down to 1e-4 of the molecules. An oxygenic class makes no methane once its oxygen has risen. The haze, ultraviolet shield, colour and habitability score are redone from the air life leaves. Classes 4, 5, 7 and 9 are `pending` rows and change nothing yet. Every number has a tier and tweak level (row LF60). The haze's cooling of the ground and the hydrosphere of a world the new air moves across freezing are not modelled.
- **Icy-moon tide (GEN.906).** An icy body on an orbit so round that its tide underflows to 0.0 now stores no tide (NULL), as a rocky body does, instead of 0.0.
- **Magnetopause in planet radii (GEN.843).** `magnetopause_rp` was stored in Earth radii but labelled planet radii, so any body not Earth-sized read wrong by the ratio of radii. It now divides the standoff by the body's own radius.

generation-output: changed

TODO items resolved: GEN.891, GEN.857, GEN.858, GEN.859, GEN.906, GEN.843. Added: stand-in registry row LF60.
