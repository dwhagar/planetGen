# Star types against distance from the galactic core (GEN.133)

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
| M dwarfs, solar neighborhood | 74-75% | 69% (10 pc census) to 75-78% | [RECONS](https://arxiv.org/pdf/1206.1022), [Mamajek](https://www.pas.rochester.edu/~emamajek/memo_star_dens.html) |
| B stars, solar neighborhood | 0.13% between arms, 0.37% in an arm | about 0.04% | [Mamajek](https://www.pas.rochester.edu/~emamajek/memo_star_dens.html) |
| A stars | 0.8-1.3% | about 0.6% | same |
| O stars | 7e-5% to 9e-4% | about 5e-5% | same |
| Bulge stars younger than 5 Gyr | 0% (bulge ages 8-12 Gyr) | 3% (HST) to 16-23% (microlensing, disputed) | [Valle et al. 2015, quoting Clarkson 2011](https://arxiv.org/pdf/1503.04570), [Bensby 2017](https://arxiv.org/pdf/1702.02971) |
| Star formation across the disk | young share flat from 2 to 15 kpc | rate peaks near 5 kpc, falls about 0.28 dex per kpc beyond; 84% inside the solar circle, 1% beyond 13.5 kpc | [Herschel Hi-GAL study](https://www.arxiv.org/pdf/2211.05573) |
| Galactic center star formation | none (bulge is old) | Central Molecular Zone holds ~5% of the Galaxy's rate | same, [Galactic center review](https://arxiv.org/pdf/1905.01309) |
| Iron gradient | none | -0.06 dex/kpc for young stars, flatter (-0.03) for 6-10 Gyr | [Anders et al. 2017](https://arxiv.org/pdf/1608.04951.pdf) |

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
4. Add a metallicity gradient only if planet occurrence is later tied to it.
