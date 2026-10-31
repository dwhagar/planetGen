**Headline:** The oxygen fugacity, dry carbon and comet species rows of the atmosphere recipes now come from models, not from Solar System bodies.

### Changed

- GEN.1005 (the flip): rows a1 and a2 of `atmosphere-a` draw the magma-ocean surface redox and the later upper mantle from `physics/magma_ocean_redox.py`; row b1's dry carbon (`mixture_rocky.INVENTORY['co2_dry_wt_percent']`) is `physics/dry_carbon.py`'s bracket, 0.029 to 0.108 weight percent CO2; row k7's N2, Ar, OCS and SO2 follow the host's own nitrogen, argon and sulfur per carbon (new subject quantities `n_over_c`, `ar_over_c`, `s_over_c`) through `physics/comet_species.py`, with a width of 0.4 dex in place of the design 0.6. The Moon, Mars, the Earth and 67P are tests of the models, except that the comet species' trapping shares and O2, Kr and Xe are tuned to 67P and stay flagged stand-ins (rows CK1 to CK6 of section 6.59, replaced by GEN.1021). a1, a2 and b1 are `partly verified`; SSR1 to SSR3 are retired; k7 stays a stand-in on GEN.1021.
- generation-output: changed
