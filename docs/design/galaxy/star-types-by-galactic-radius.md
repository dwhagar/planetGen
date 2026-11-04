# Star types against distance from the galactic core (GEN.133)

**Status (2026-10-17):** analysis (GEN.133, #807); proposals 1 and 3 built derived as GEN.134 (see "Built: GEN.134" at the end), proposal 2 not met (see there), proposal 4 already built by GEN.158.

Analysis of 2026-10-09 for Boss: how the generator's star-type rates change with
galactocentric radius, set against what is observed. Method: `scripts/star_type_by_radius.py`
draws 30,000 stars per population from the generator's own model
(`stellar_evolution.sample_living_star`) and mixes them with the default galaxy's
population shares at each radius (`density.population_densities`). Nothing was changed
by this note; the proposals at the end are for Boss to choose from.

## How the generator works today

A star's type comes only from its population (young, intermediate, old, bulge: an age
range each) through a Kroupa mass function and stellar evolution. Position enters only
by changing the population mix: arms and the thin plane are young, the bulge and the thick
disk and halo are old. There is no metallicity, no radial star-formation history and no
radial change in the mass function. So the class shares move only where the mix moves:
fractions of B, A and white dwarfs change, and M, K and G stay flat (M 74-75%, K 13-14%, G
2.8-3.4% everywhere).

## Results

### Class shares inside each population (main-sequence dwarfs by class)

| population | O | B | A | F | G | K | M | subgiant | giant | wd |
|---|---|---|---|---|---|---|---|---|---|---|
| young | 0.0133% | 2.99% | 2.65% | 2.86% | 3.22% | 13.9% | 74.3% | 0.02% | 0% | 0.05% |
| intermediate | 0% | 0.51% | 2.09% | 3.12% | 3.34% | 13.8% | 74.2% | 0.197% | 0.143% | 2.6% |
| old | 0% | 0% | 0.0667% | 1.44% | 3.37% | 13.3% | 74.7% | 0.347% | 0.357% | 6.42% |
| bulge | 0% | 0% | 0% | 0.06% | 2.84% | 14% | 74% | 0.55% | 0.417% | 8.12% |

### spiral-arm crest, z = 0 pc

| R (kpc) | relative density | young | intermediate | old | bulge | O | B | A | F | G | K | M | subgiant | giant | wd |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.25 | 160 | 2% | 15% | 15% | 69% | 0.000273% | 0.135% | 0.369% | 0.765% | 3% | 13.8% | 74.2% | 0.458% | 0.36% | 6.9% |
| 0.5 | 153 | 2% | 14% | 14% | 71% | 0.000258% | 0.128% | 0.349% | 0.726% | 2.99% | 13.8% | 74.2% | 0.463% | 0.363% | 6.97% |
| 1 | 112 | 2% | 16% | 16% | 67% | 0.000292% | 0.145% | 0.396% | 0.814% | 3.01% | 13.8% | 74.2% | 0.452% | 0.356% | 6.82% |
| 2 | 27.2 | 6% | 44% | 43% | 8% | 0.000817% | 0.406% | 1.1% | 2.16% | 3.31% | 13.6% | 74.4% | 0.277% | 0.246% | 4.49% |
| 3 | 19.1 | 6% | 42% | 41% | 11% | 0.00079% | 0.392% | 1.07% | 2.08% | 3.29% | 13.6% | 74.4% | 0.287% | 0.252% | 4.62% |
| 4 | 11.5 | 7% | 48% | 46% | 0% | 0.000893% | 0.444% | 1.21% | 2.34% | 3.34% | 13.6% | 74.4% | 0.253% | 0.231% | 4.17% |
| 5 | 7.8 | 7% | 48% | 45% | 0% | 0.000898% | 0.446% | 1.21% | 2.34% | 3.34% | 13.6% | 74.4% | 0.253% | 0.23% | 4.16% |
| 6 | 5.36 | 7% | 48% | 44% | 1% | 0.00089% | 0.442% | 1.2% | 2.32% | 3.34% | 13.6% | 74.4% | 0.256% | 0.232% | 4.2% |
| 8 | 2.43 | 7% | 49% | 45% | 0% | 0.000908% | 0.451% | 1.23% | 2.35% | 3.34% | 13.6% | 74.4% | 0.252% | 0.229% | 4.13% |
| 10 | 1.12 | 7% | 49% | 44% | 0% | 0.000914% | 0.454% | 1.23% | 2.36% | 3.34% | 13.6% | 74.4% | 0.251% | 0.228% | 4.12% |
| 12 | 0.517 | 7% | 49% | 44% | 0% | 0.000918% | 0.456% | 1.24% | 2.36% | 3.34% | 13.6% | 74.4% | 0.251% | 0.228% | 4.11% |
| 15 | 0.162 | 7% | 49% | 44% | 0% | 0.000923% | 0.458% | 1.24% | 2.37% | 3.34% | 13.6% | 74.4% | 0.25% | 0.227% | 4.1% |

### midway between arms, z = 0 pc

| R (kpc) | relative density | young | intermediate | old | bulge | O | B | A | F | G | K | M | subgiant | giant | wd |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.25 | 131 | 0% | 5% | 13% | 82% | 1.16e-05% | 0.0283% | 0.116% | 0.391% | 2.93% | 13.9% | 74.1% | 0.506% | 0.395% | 7.62% |
| 0.5 | 103 | 0% | 6% | 14% | 80% | 1.33e-05% | 0.0325% | 0.133% | 0.439% | 2.94% | 13.9% | 74.1% | 0.5% | 0.392% | 7.55% |
| 1 | 56.1 | 0% | 9% | 22% | 69% | 2.02e-05% | 0.0494% | 0.203% | 0.634% | 3% | 13.8% | 74.2% | 0.474% | 0.379% | 7.25% |
| 2 | 27.2 | 0% | 12% | 30% | 57% | 2.84e-05% | 0.0695% | 0.285% | 0.859% | 3.06% | 13.8% | 74.3% | 0.444% | 0.364% | 6.91% |
| 3 | 7.78 | 1% | 29% | 70% | 0% | 6.76e-05% | 0.165% | 0.676% | 1.94% | 3.36% | 13.5% | 74.5% | 0.301% | 0.292% | 5.26% |
| 4 | 9.84 | 0% | 16% | 37% | 47% | 3.64e-05% | 0.0889% | 0.363% | 1.06% | 3.11% | 13.7% | 74.3% | 0.417% | 0.35% | 6.6% |
| 5 | 3.52 | 1% | 30% | 69% | 0% | 6.92e-05% | 0.169% | 0.69% | 1.95% | 3.36% | 13.5% | 74.5% | 0.3% | 0.291% | 5.24% |
| 6 | 2.37 | 1% | 30% | 69% | 0% | 6.99e-05% | 0.171% | 0.696% | 1.96% | 3.36% | 13.5% | 74.5% | 0.299% | 0.29% | 5.22% |
| 8 | 1.08 | 1% | 31% | 69% | 0% | 7.11e-05% | 0.174% | 0.707% | 1.97% | 3.36% | 13.5% | 74.5% | 0.299% | 0.289% | 5.2% |
| 10 | 0.495 | 1% | 31% | 68% | 0% | 7.2e-05% | 0.176% | 0.715% | 1.97% | 3.36% | 13.5% | 74.5% | 0.298% | 0.288% | 5.19% |
| 12 | 0.227 | 1% | 32% | 68% | 0% | 7.28e-05% | 0.178% | 0.722% | 1.98% | 3.36% | 13.5% | 74.5% | 0.297% | 0.287% | 5.17% |
| 15 | 0.0707 | 1% | 32% | 67% | 0% | 7.37e-05% | 0.18% | 0.73% | 1.99% | 3.36% | 13.5% | 74.5% | 0.297% | 0.286% | 5.16% |

(The full output, including z = 500 pc, is what the script prints.)

## Against the observed universe

| Quantity | Generator | Observed | Source |
|---|---|---|---|
| M dwarfs, solar neighborhood | 74-75% | 69% (10 pc census) to 75-78% | [R1112], [Mamajek](https://www.pas.rochester.edu/~emamajek/memo_star_dens.html) |
| B stars, solar neighborhood | 0.13% between arms, 0.37% in an arm | about 0.04% | [Mamajek](https://www.pas.rochester.edu/~emamajek/memo_star_dens.html) |
| A stars | 0.8-1.3% | about 0.6% | same |
| O stars | 7e-5% to 9e-4% | about 5e-5% | same |
| Bulge stars younger than 5 Gyr | 0% (bulge ages 8-12 Gyr) | 3% (HST) to 16-23% (microlensing, disputed) | [R0346, R4483] |
| Star formation across the disk | young share flat from 2 to 15 kpc | rate peaks near 5 kpc, falls about 0.28 dex per kpc beyond; 84% inside the solar circle, 1% beyond 13.5 kpc | [Herschel Hi-GAL study][R1241] |
| Galactic center star formation | none (bulge is old) | Central Molecular Zone holds ~5% of the Galaxy's rate | same, [R0141] |
| Iron gradient | none | -0.06 dex/kpc for young stars, flatter (-0.03) for 6-10 Gyr | [R0101] |

## Findings

1. Local class shares of M, K, G, A and O are inside the observed spread. B stars come out
   three to ten times too common at the Sun's radius (the young population's 3% B share times
   its arm weight).
2. Nothing varies with radius beyond the population mix. The young share stays 7% on an arm
   from 2 kpc to 15 kpc, so the outer galaxy has as many O, B and A stars per star as the Sun's
   neighborhood, when the real star formation falls by about a factor of 20 out to 15 kpc.
3. The inner galaxy (under 1 kpc) is all bulge and old disk: no O or B stars at all in the
   plane, which matches the bulge but leaves out the Central Molecular Zone's young stars (about
   5% of the Galaxy's star formation, a few hundred parsecs across).
4. White dwarfs run 4-8% (higher where it is older), plausible against the roughly 5% local figure (not separately sourced here).
5. The bulge has no stars under 8 Gyr; the observed fraction is 3% to about 20% depending on the method.
6. No metallicity: the observed -0.06 dex/kpc gradient changes planet occurrence and the giant
   branch, but not the main-sequence class shares to first order.

## Proposals (not built)

1. Weight the young and intermediate populations by the observed star-formation profile (peak at
   5 kpc, about -0.28 dex per kpc beyond, a dip inside 3 kpc, the Central Molecular Zone as its own
   small young region).
2. Lower the young population's B share, or cut its weight, to bring local B stars to about 0.04%.
3. Give the bulge a small young tail (about 10% under 5 Gyr, between the HST and microlensing figures).
4. Add a metallicity gradient only if planet occurrence is later tied to it. Globular clusters are a second reason to add a metallicity value: their giant-planet cull is `10^(2 [Fe/H])` ([globular-clusters.md](globular-clusters.md) section 7).

## How the proposals are built (Boss, 2026-10-16)

Boss's NO FIXED CONTENT ruling (2026-10-16, 4:11 and 4:14 PM EDT) applies to all four proposals above. The numbers in them
(a peak at 5 kpc, a slope of -0.28 dex per kpc, a dip inside 3 kpc, B stars at 0.04%, a bulge tail of 10% under 5 Gyr) are
the Milky Way's observed values. They are test targets that show the derived result is plausible, never inputs to the
generator and never a preset. GEN.134 is built derived: the star-formation profile by radius follows from the galaxy's own
shape (its disk scale length, its bulge and bar, and the gas supply those imply), the B share follows from the mass function
and the star-formation history, and the bulge's young tail follows from the bulge's own formation model. Each value is
calculated or drawn from the existing seeds, with real science named for every link, and the Milky Way figures above are only
what the plausibility test compares the result with.

Boss chose Core 4 to take GEN.134 first and to research the inner-disc dip (2026-10-16, 8:40 PM EDT, card "Core 4 first").
Core 4 asks Boss only for what it cannot look up or verify.

## Built: GEN.134 (2026-10-17)

`galaxy/star_formation.py` gives the young and intermediate populations a radial factor, applied in
`density.population_densities` and `population_densities_array`. It is derived from the shape's own disk scale length
and bulge scale radius. The Milky Way's figures are test targets in `tests/test_star_formation_profile.py`, never inputs.

**The chain.**

1. *The star-formation law.* Rate per area goes as the gas per area to the power n = 1.4 ([R2321]: 1.4 +/- 0.15
   over 61 spirals and 36 starbursts, gas 10 to 1e5 solar masses per pc2; other fits 1.3 +/- 0.3 and 1.65 +/- 0.16;
   resolved data give 1.5 to 2.7 for total gas). The thin disk's stars are the sum of its past rates, so their surface
   density is `exp(-R/L)`. A gas disk of scale length `g L` makes the present rate fall as `exp(-n R / (g L))` and the
   *share* of the stars that are recent as `exp(-(n/g - 1) R / L)`. `g` is 1 (`GAS_TO_STELLAR_SCALE_LENGTH_RATIO`, an
   assumption: the observed gas disk is more extended, which would make the outer disk younger). The Milky Way's own Kennicutt-Schmidt index
   is lower than the classic one: 1.2 +-0.2 ([R2770], molecular gas) and 1.12 +-0.05 [R4140], as the
   Herschel paper below reports them [P: as that paper quotes them]; with them the outer fall would be -0.19 dex per kpc, and the observed -0.28 needs an index of 1.65.
   *The index rises where the gas thins (GEN.450, Boss's answer to the card of 12:27 PM EDT 2026-10-23).* Instead of one index of 1.65, which is
   chosen to hit the observed figure (a test target, never an input), the index steepens where the disk's own gas is thin: [R0392]
   find that molecular gas forms stars at a constant efficiency (index 1.0) while atomic gas, which dominates the outskirts, saturates near 9
   solar masses per pc2 and shows a weak total-gas relation; [R2651] find the efficiency in the atomic-dominated
   disk falls exponentially with radius; [R3951] gives a cold-phase threshold. The rate is
   `A gas^n (gas / (gas + gas_s))^(m - n)`: `n` stays 1.4, `gas_s` is 9 solar masses per pc2 (`tuning.STAR_FORMATION_SATURATION_GAS_MSUN_PC2`) and `m`
   is 2.1 (`tuning.STAR_FORMATION_LOW_DENSITY_EXPONENT`), the middle of the 1.5 to 2.7 that resolved total-gas fits give. The gas is the rings' gas now
   (`disk_history.history_of`), about 4 solar masses per pc2 at 8 kpc for a Milky Way-like galaxy (the observed total is higher, about 13, so the
   threshold bites a little earlier in the model than in the Galaxy), and the extra factor is normalised away like the rest, so the disk still makes
   as many stars. The result for the Milky Way-like shape is a fall of -0.27 dex per kpc beyond the peak. What is open: `m` and `gas_s` are not read from
   one galaxy (no source gives the index of a single galaxy's outer disk), the smooth form `gas / (gas + gas_s)` is a choice, and the ring gas is
   lower than the observed total.
2. *The bar gap.* A bar drives the gas inside its corotation inward, and what is left is sheared and made turbulent by the bar's
   shocks, so the star-forming gas in the region the bar sweeps is depleted. What the papers say [P: each read in the paper, 2026-10-17]:
   in six zoom-in simulations the star formation "desert" either side of the bar lacks stars younger than 10 Myr because, after
   the bar forms, "gas in these regions is removed on 1 Gyr timescales" [R1142];
   bars drive gas inflow inside the corotation radius, which raises the central rate, yet disc galaxies with a suppressed centre often host
   a bar, and the causes are discussed, not settled [R4610]; galaxies with bar-driven radial gas flows have
   suppressed star-formation efficiency with "no indication of depleted gas mass fractions", so gas exhaustion is not the cause there
   [R1973]. Observed in the Milky Way: the molecular gas shows "a pronounced deficit at R about 0.5 to 3.5 kpc (a bar gap)", an outer boundary at
   3.5 kpc that "is consistent with the dynamical estimate of the Galactic bar half-length", with the caveat that axisymmetric distance
   inversions can partly fill such a gap [R0215]; the long bar's half-length is
   5.0 +-0.2 kpc in [R4650], so the estimates differ. The Central Molecular Zone forms stars at
   a rate at least ten times below what its dense gas predicts ([R2505], after [Longmore et al.
   2013][R2737]): episodic, with a raised critical density from the turbulent pressure. The generator's rate is scaled by
   `S(R) = f + (1 - f) u^p / (1 + u^p)`, `u = R / R_gap`. **The gap radius** is 2.2 bulge scale radii (`BAR_GAP_RADIUS_RATIO`), 3.5 kpc for the default shape, which
   is the gap's outer edge over the bulge's 1.58 kpc [P for the 3.5 kpc]; **the floor** `f` (0.03) is what 1 Gyr removal leaves after 3.5 Gyr
   (`exp(-3.5)`), the young and intermediate stars being under 3 Gyr old; that bar epoch is an assumption and the 1 Gyr is a simulated galaxy's, so the floor stays an open default;
   **the edge power** `p` (4) is a model default (2 never leaves a ring). Fast bars have corotation over bar length of 1.0 to 1.4 (Debattista and Sellwood
   2000; a search summary, the paper was not opened).
3. *Normalisation.* The factor is scaled so the stars of the disk, weighted by their own surface density, average 1: the
   disk makes as many stars over its life as the uniform-in-time populations already assumed. The three disk populations
   still share the same disk density, so no cell's total changes.

**Result for the default shape** (disk scale length 2.6 kpc, bulge scale radius 1.58 kpc; rate per area, test values):

| Quantity | Derived | Milky Way target |
|---|---|---|
| Peak of the rate per area | 3.3 kpc (the gap's edge) | 4 to 5 kpc (molecular ring) |
| Share of the star formation inside the solar circle | 91% (86% before GEN.450) | 84% |
| Share beyond 13.5 kpc | 0.4% (1.3% before GEN.450) | 1% |
| Fall beyond the peak | -0.27 dex per kpc (the plain law gives -0.23 = 1.4 x 0.434 / 2.6; the index rising where the gas thins, GEN.450, adds -0.04) | -0.28 |
| Share inside 3 kpc | 13% | a dip (not quantified here) |
| Stars under 5 Gyr in the bar's core (250 to 1000 pc) | 9% to 11% | 3% (HST) to 16 to 23% (microlensing, disputed) |

The young share in the bar's core falls from 2% to under 1%, and the intermediate from 15% to 4%, on an arm crest at 0.25 kpc (the old disk's
part of the thin disk and the bulge fill the rest).

**Proposal 2 is not met.** The B-star share at the Sun's radius between arms is 0.2% (was 0.17%) and 0.5% in an arm,
against the note's 0.04%. The derived radial factor is 0.96 at the solar circle (the Sun sits at the disk's average), so
it cannot lower the local B share. B stars are the stars formed in the last few hundred Myr, so the share follows the
*recent* rate against the lifetime average: the generator forms stars at a constant rate over 10 Gyr, and a recent rate
four times lower is a time dependence that the radial factor does not model. The next derived step is the disk's
star-formation *history* from the same gas law (the rate falls as the gas is used, unless it is replenished), which needs
an inflow rate that nothing in the shape sets; filed as a follow-up. The 0.04% target also rests on a source this note has
not opened at the time (Mamajek's memo; read since, see Citations), and is his figure from Hipparcos B types.

**Citations (re-checked against the papers, 2026-10-17, 12:40 AM EDT).** [R2321]: 61 normal spirals and 36 starbursts, index 1.4 +-0.15, confirmed.
[R1241] (Milky Way star formation 2.0 +-0.7 solar masses per year; a peak at the Central Molecular Zone and another near 5 kpc; beyond 5 kpc
log(rate per area) falls with a = -0.28 +-0.01 per kpc, an exponential scale length of 1.55 kpc, shorter than the stellar disk's 2.6 kpc; 84% inside the solar circle and 1% beyond 13.5 kpc; half of the star formation in the Molecular Ring): all confirmed, and the 1.55 kpc against 2.6 kpc is the same
sign as the derived 2.6 / 1.4 = 1.86 kpc. The remaining papers are listed in the chain above with what each says. The B-star target is [Mamajek's memo](https://www.pas.rochester.edu/~emamajek/memo_star_dens.html), read: B stars 3.9e-4 of the stars (1 in 2600; calculated from the Hipparcos (HIP2) B types), A 6.0e-3, O about 5e-7, M dwarfs
72.5% (RECONS 10 pc, N = 283); the 0.04% is confirmed as his figure, and its completeness limits are his. The numbers marked recalled in `tuning.py` stay so, except the gap edge (3.5 kpc, read).

**Out of scope here.** The Central Molecular Zone (GEN.282), metallicity (GEN.158, done), and the globular clusters.

**The majorant.** The object-first band bound (`generation.object_first._disk_share_bound`) widens the young and
intermediate weights by the factor's true smallest and largest value over the band's cylindrical radii
(`recent_factor_range`). The factor is not monotone (it dips a little from the centre, rises to the gap's edge and falls
beyond), so the extremes are at an end or at one of its two turning points, found by bisection.


## Built: GEN.469 (2026-10-17)

The disk-wide shares of the young, intermediate and old populations, which GEN.134's radial factor multiplies, now come
from the disk's own gas history (`galaxy/disk_history.py`, [disk-star-formation-history.md](disk-star-formation-history.md))
instead of the age windows' lengths over ten Gy. For a galaxy of 5e10 solar masses the local B-star share between the
arms falls from 0.19% to 0.12%; the observed 0.04% stays a test target and a strict expected failure.
