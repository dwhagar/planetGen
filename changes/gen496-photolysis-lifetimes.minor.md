**Headline:** The photolysis lifetime of a molecule can be worked out for its star and its distance, species by species.

### Changed
- **The photolysis lifetime table (GEN.496).** The new `physics/photolysis.py` and `physics/data/photolysis_lifetimes.json` hold, for H2O, CO, CO2, CH4, C2H6, NH3, HCN, CH3OH, CH2O, H2S and N2, the lifetime at 1 AU from the quiet Sun, the range the sources give, the active-Sun factor where one was read (H2O 1.68, CO2 2.35, NH3 1.075) and the rate ratio around a 4000 K and a 10000 K star to the Sun's at equal ultraviolet flux. The rate follows the flux and, species by species, the star's spectrum: around a 4000 K star CO photolyses 3.5e-4 as fast as around the Sun and N2 2.9e-5, while NH3 hardly changes, so one scalar for "stellar UV" would be wrong by orders of magnitude. NH3 is the quiet-Sun 6.7e3 s, not the 5.6e3 s mean of quiet and active Sun. Stars outside 4000 to 10000 K take the end template and are flagged; shielding is not included. Nothing is stored and no generated value changes.

generation-output: unchanged
