# The dry-body carbon model (GEN.1005, part 1 of row b1)

Status: module and registry rows, written 2026-10-30 (as the lane reported) by the Aux 3 lane. `src/planetgen/physics/dry_carbon.py`; tests in `src/tests/test_dry_carbon.py`; registry rows DC1 to DC3 (section 6.58 of the gap-filling plan). Part 1 changed nothing the generator draws. Part 2 (the digest flip, 2026-10-30) has `mixture_rocky.INVENTORY['co2_dry_wt_percent']` read `co2_wt_percent_bracket()`; SSR3 is retired and row b1 is `partly verified`.

## 1. What replaces what

Row b1 (a body that formed inside the snow line and took in little water) carried a log bracket of 0.018 to 0.1 weight percent CO2: the stored dry value of [R3378] to its high case, chosen so that an Earth-distance draw's median sat at the bulk silicate Earth's 100 ppm of carbon (`volatile-inventory-against-planet-totals.md`, GEN.813). The Earth was the anchor. This model derives the bracket from the carbon of the building blocks and what core formation leaves, and uses the Earth as a test.

## 2. The three steps (all from Blanchard, Siebert, Badro, EPSL, [R0419]; read in full text)

| Step | Number | Where | Registry |
|---|---|---|---|
| Carbon of a non-carbonaceous building block | 1600 ppm undifferentiated, 40 ppm differentiated (it left with the core) | Table 2: 0.16 and 0.004 weight percent | DC1, T1, L0 |
| Share of the accreted mass that had cored before it was accreted | 0.30 to 0.60 | the accretion model: 45 +/- 15 percent with 3.35 weight percent carbonaceous carbon, 25 +/- 15 with 2.2 | DC2, T7, L1 |
| Share of the carbon left in the mantle and air | 0.12 to 0.26 | 140 +/- 40 ppm in the mantle against 495 +/- 125 ppm of the whole, mantle 0.68 of the mass: 0.19 with the errors carried | DC3, T7, L2 |

Multiplying the extremes and converting carbon to CO2 (3.664) gives 0.029 to 0.108 weight percent of the body, 0.6 dex inside the old ends on both sides. The old bracket was the Earth; this one comes from the planetesimals' carbon and the incomplete equilibration of the core.

Why the mantle keeps a fifth: Eq. 4 of the paper (carbon as C2+) is log K_D = 0.3 + 3822 / T, and D_C = K_D x D_Fe with D_Fe = 10^(-dIW/2); at 3800 K and dIW -1.06 that is about 70 (the paper's runs span 23 to 157), so a core in full equilibrium with the whole mantle would take 97 percent of the carbon. The paper's own answer is that only 0.16 to 2.5 percent of the mantle equilibrates with a planetesimal's core and 2.6 to 9.9 percent with an embryo's, "C is not extracted efficiently from the mantle during core formation in spite of its siderophile behaviour". `test_dry_carbon.py` holds the partitioning arithmetic as a test of the share (full equilibrium would leave less than the bracket's low end).

## 3. What it does not do

- The three numbers of steps 2 and 3 are results of one accretion model of one planet (Earth's mass, three-stage core formation). They are the centres of a bracket for bodies whose accretion is like that planet's. A body that accretes its metal late, from carbonaceous bodies, or round a carbon-richer star is outside it. DC2 and DC3 are T7 for that reason.
- The star's carbon (carbon per silicon of the disc) scales step 1 and is not in the model; the subject of the rocky mixture does not carry it yet.
- The core's size, the ocean's depth and the oxygen fugacity enter step 3 through Eq. 4 and are not drawn; carbon becomes less siderophile with depth along a magma-ocean geotherm, which only raises the mantle's share.
- Carbonaceous (wet) material is not handled here; rows b2 to b4 carry that.
- Li, [R2664] (Science, [R2664]; read) puts the carbon of the bulk Earth at 0.4 +/- 0.2 weight percent (upper bound 0.2 to 1.9) and has the refractory carbon of the disc sublimate near 500 K (the soot line) with a drop of more than a decade in condensed carbon inside it. It supports step 1's low carbon for material that formed hot but gives no formula for a planet's carbon against formation distance, so the soot line is not drawn.
