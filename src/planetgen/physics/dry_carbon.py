# planetgen/physics/dry_carbon.py

"""
The carbon a dry rocky body keeps in its mantle, from the carbon of its building blocks and what the core leaves (GEN.1005, part 2 of the redox, dry-carbon and comet rows;
docs/design/composition/dry-carbon-model.md).

Row b1 of the atmosphere possibility table (a body that formed inside the snow line and took in little water) drew the carbon of its outgassing mantle as a log bracket of 0.018 to
0.1 weight percent CO2, tuned to the bulk silicate Earth's 100 ppm of carbon (GEN.813). This module derives the bracket from three steps, none of which names a planet:

1. **What the building blocks hold.** Blanchard, Siebert, Badro, Ghiorso and Hirschmann 2022 (EPSL, arXiv:2202.06809, Table 2, READ) give the carbon of undifferentiated
   non-carbonaceous bodies as 0.16 weight percent and of differentiated non-carbonaceous bodies as 0.004 (the carbon left the planetesimal with its core): `NC_UNDIFFERENTIATED_C_PPM`,
   `NC_DIFFERENTIATED_C_PPM`. [Stand-in DC1, tier T1, tweak L0.]
2. **How much of the material was already differentiated.** In their accretion model a share of the planetesimals that built the planet had cored before they were accreted: 45 +/- 15
   percent reproduces the mantle's carbon when the carbonaceous bodies carry 3.35 weight percent, and 25 +/- 15 percent with 2.2. The bracket here is 30 to 60 percent (`DIFFERENTIATED_SHARE`),
   the same width on the better-fit centre. [Stand-in DC2, tier T7, tweak L1.]
3. **How much of the carbon the core leaves in the mantle.** Carbon is siderophile (Eq. 4 of the paper: log K_D = 0.3 + 3822 / T for carbon as C2+, D_C = K_D x D_Fe, which gives D_C of 20 to 140
   at 3600 to 4000 K, the range of their experiments), but "C is not extracted efficiently from the mantle during core formation": only a few percent of the mantle equilibrates with each
   impactor's core (0.16 to 2.5 percent for planetesimal impacts, 2.6 to 9.9 for embryos), and the model leaves 140 +/- 40 ppm in the mantle of 495 +/- 125 ppm of the whole. The share of the body's carbon that
   stays in the mantle and the air is about 19 percent, with a spread of 12 to 26 (`MANTLE_SHARE`). [Stand-in DC3, tier T7, tweak L2.]

`co2_wt_percent_bracket` multiplies the three at their extremes and writes it as the weight percent of CO2 in the body (`CO2_PER_C`): 0.029 to 0.108, against the 0.018 to 0.1 that row b1 draws.

Limits. The three numbers of steps 2 and 3 are results of one accretion model of one planet (it is run with Earth's mass and its three-stage core formation), so they are the centres of a
bracket for rocky bodies whose accretion is like that planet's; a body that accretes more of its metal late, or from carbonaceous bodies, or round a carbon-richer star is outside it. The star's own
carbon (the share of carbon per silicon of the disc) scales step 1 and is not in the model. The core's size, the depth of the ocean and the oxygen fugacity enter step 3 through the partitioning of
Eq. 4 and are not drawn: the model carries one share for all of them (Blanchard et al.: carbon becomes less siderophile with depth along a magma-ocean geotherm, which only raises the mantle's share).
"""

CO2_PER_C = 44.009 / 12.011
"""float: Mass of CO2 per mass of carbon."""

NC_UNDIFFERENTIATED_C_PPM = 1600.0
"""float: Carbon of an undifferentiated non-carbonaceous building block, ppm by weight (Blanchard et al. 2022 Table 2: 0.16 weight percent). [Stand-in DC1, tier T1, tweak L0.]"""

NC_DIFFERENTIATED_C_PPM = 40.0
"""float: Carbon of a non-carbonaceous planetesimal that cored before it was accreted, ppm by weight (the same table: 0.004 weight percent). [Stand-in DC1, tier T1, tweak L0.]"""

DIFFERENTIATED_SHARE = (0.30, 0.60)
"""tuple: The share of the accreted mass that was already differentiated, low and high (Blanchard et al. 2022: 45 +/- 15 percent for 3.35 weight percent carbonaceous carbon). [Stand-in DC2, tier T7, tweak L1.]"""

MANTLE_SHARE = (0.12, 0.26)
"""tuple: The share of a body's carbon that stays in its mantle and air after core formation, low and high (140 +/- 40 ppm of the mantle in 495 +/- 125 ppm of the whole, mantle 0.68 of
the mass: 0.19, with the two errors carried: 0.12 to 0.26). [Stand-in DC3, tier T7, tweak L2.]"""


def bulk_carbon_ppm(differentiated_share):
    """The carbon of the accreted material, ppm by weight, when `differentiated_share` of it had cored before it was accreted."""
    if not 0.0 <= differentiated_share <= 1.0:
        raise ValueError(f"a share of the accreted mass is between 0 and 1, not {differentiated_share}")
    return (1.0 - differentiated_share) * NC_UNDIFFERENTIATED_C_PPM + differentiated_share * NC_DIFFERENTIATED_C_PPM


def mantle_carbon_ppm(differentiated_share, mantle_share):
    """The carbon of the mantle and air of the whole body, ppm by weight."""
    if not 0.0 <= mantle_share <= 1.0:
        raise ValueError(f"the share of the carbon that stays in the mantle is between 0 and 1, not {mantle_share}")
    return bulk_carbon_ppm(differentiated_share) * mantle_share


def co2_wt_percent_bracket():
    """The `(low, high)` weight percent of the body that is CO2 in a dry body's mantle and air: the least differentiated and the lowest mantle share against the most and the highest."""
    low = mantle_carbon_ppm(max(DIFFERENTIATED_SHARE), min(MANTLE_SHARE))
    high = mantle_carbon_ppm(min(DIFFERENTIATED_SHARE), max(MANTLE_SHARE))
    return low * CO2_PER_C / 1e4, high * CO2_PER_C / 1e4
