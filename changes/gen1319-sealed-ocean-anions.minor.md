**Headline:** A sealed ocean now holds the chloride and sulfate its volcanoes gave it, so its charge balance is no longer sodium and calcium against hydroxide alone; the pH of such oceans comes down and a few more of them keep carbon dioxide in the air (GEN.1319).

### Changed

- `waterworld_carbon.ocean_pco2` takes `anion_charge_mol_m2` (chloride once, sulfate twice) into the charge balance; with no anions it gives the old result.
- `surface_air._shares` passes the chloride the volcanoes degassed (which the air still loses to the ocean) and the sulfate of `sulfur_balance` to it (stand-ins AN1 and AN2, section 6.197).
- Measured on 100 habitable-zone G2V and 100 K3V stars: about 1 to 3 percent more sealed-ocean worlds hold 5 Pa of carbon dioxide; the median sealed-ocean pH falls from 11.0 to 10.7 (G2V) and stays near 12 (K3V). The leached calcium, not the missing anions, locks most of the carbon in calcite.
