**Headline:** The generator now has one oxygen-fugacity scale, so a delta-IW means the same thing everywhere.

### Changed
- **One oxygen-fugacity scale (GEN.498).** The new `physics/oxygen_fugacity.py` fixes the reference as the Hirschmann 2021 iron-wustite buffer; every stored delta-IW, the mantle redox of a planet included, is on it. It converts the Frost 1991 buffers (QFM stands 3.1 to 3.6 above IW at 1 bar from 2200 to 1500 K) and the Sossi 2020 scale (his IW+0.5 for Earth is +0.0 on the reference) to it, gives the CO2/CO and H2O/H2 gas ratios at any delta-IW from the NASA Glenn fits of the five gases (CO2 equals CO at IW+1.37 and H2O equals H2 at IW-0.06 at 2173 K), and Sossi's melt Fe3+ law both ways. The coefficients' signs of the reference fit and the Frost constants are second hand and marked so in the module. No stored value changes.

generation-output: unchanged
