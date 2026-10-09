# Neutron stars and black holes by region and type: what is observed

Research of 2026-10-09 (GEN.100 rate retune; rates are in `interstellar-object-rates.md`). Found by web search; numbers are as the sources state them. Where a source gave no figure I say so instead of filling one in.

## What the sources say

**Distance from the core (radial)**
- Radio pulsars are not distributed like stars. Faucher-Giguere & Kaspi (2006) put the pulsar surface density peaking near 3 kpc from the center, born in spiral arms ([arXiv astro-ph/0512585](https://arxiv.org/abs/astro-ph/0512585)). A 2024 re-fit of 1,301 Parkes Multibeam pulsars peaks near 4 kpc ([arXiv 2402.11428](https://arxiv.org/pdf/2402.11428)). Lorimer et al. (2006) warn the inner profile depends heavily on the assumed electron model, so the inner-Galaxy dip is partly an observing effect ([arXiv astro-ph/0607640](https://arxiv.org/pdf/astro-ph/0607640)).
- Population synthesis of all compact remnants (Sweeney et al.) finds the black-hole-to-neutron-star ratio rises toward the Galactic center, because black holes get smaller kicks ([arXiv 2210.04241](https://arxiv.org/pdf/2210.04241)).
- The central parsecs hold a black hole cusp: Hailey et al. 2018 infer about 10,000 stellar-mass black holes within a few light years of Sgr A*, from 12 quiescent X-ray binaries. Indirect, large uncertainty ([Columbia summary](https://news.columbia.edu/node/149)).
- The Galactic center gamma-ray excess is attributed by many to thousands of millisecond pulsars in the bulge; others dispute it ([arXiv 2106.00222](https://arxiv.org/pdf/2106.00222)).

**Height above the plane (vertical)**
- Compact remnants have a scale height about 3x that of visible stars, 1260 pc, mostly from birth kicks. About 40% of neutron stars and 2% of black holes get enough speed to leave the Galaxy ([arXiv 2210.04241](https://arxiv.org/pdf/2210.04241)).
- A second simulation gets 786 pc for black holes against 306 pc for visible stars ([arXiv 2607.22814](https://arxiv.org/pdf/2607.22814)).
- Neutron star X-ray binaries sit higher than black hole ones ([Repetto et al.](https://arxiv.org/pdf/1701.01347)).
- Radio pulsars: exponential scale height about 330 pc ([Lorimer et al. 2006](https://arxiv.org/pdf/astro-ph/0607640)).
- Magnetars: only 20 to 31 pc, far thinner than radio pulsars, consistent with the most massive young progenitors ([McGill Magnetar Catalog](https://arxiv.org/pdf/1309.4167)).

**By type**
- Millisecond pulsars: 340 known in 45 globular clusters; cluster MSPs are about 10x over-abundant per star against the field. Estimates of the hidden cluster population run 1,000 to 4,700. No field total was found ([arXiv 2412.05220](https://arxiv.org/pdf/2412.05220), [arXiv 2111.08153](https://arxiv.org/pdf/2111.08153)).
- Magnetars: no total or fraction of all neutron stars was found in these searches. Only a local density of recently born ones, about 1.3e-2 per kpc^3.

## Not found

No source gave a clean "percent above or below the mean" by radius for neutron stars or black holes by type. What exists is the shape (where pulsars peak, how thick the layers are, which kinds move farthest). Turning that into multipliers needs a model choice, not a lookup.

## What the generator does now

Neutron stars and black holes follow the stellar density map exactly: same radial profile, same thickness, same arm contrast as ordinary stars, one flat rate per star.

## Proposal (not built)

1. Thicken the layer: neutron stars about 3x the stellar scale height, black holes about 2.5x, magnetars thin (about 30 pc, arms only).
2. Black hole share rising toward the core, plus a central cusp.
3. Active radio pulsars weighted to arms, peaking at 3 to 4 kpc; millisecond pulsars weighted toward the bulge.
4. Net neutron star count about 40% lower than formed, black holes 2% lower (kick losses).

## How the generator uses this (GEN.132)

The scatter no longer places neutron stars and black holes in simple
proportion to the stellar density. `galaxy/remnant_distribution.py` multiplies
each bin's density by a factor for the kind, at the bin's centerline point:

- **Height.** Neutron stars are spread over 3x the thin disk's scale height
  and black holes over 2.5x (`tuning.REMNANT_SCALE_HEIGHT_RATIO`), as a ratio
  of sech^2 profiles, capped at 4x. In the plane that is about 0.35x the
  star-proportional rate; a few scale heights up it is 3x or more.
- **Radius.** Black holes get `1 + 1.0 * exp(-R / 2 kpc)` times the stars'
  density: twice the rate at the core, 1.02x at the Sun.
- **Pulsars.** The share of neutron stars that are active pulsars is scaled
  by the Lorimer et al. (2006) radial profile over the stellar exponential,
  1 at the Sun, about 0.1 at 1 kpc and 0.5 at 3 kpc.

These sizes are **estimates** drawn from the trends above; the sources give
shapes, not percentages, so the constants in `tuning.py` are the knobs. Only
these two kinds have regional research; planetary nebulae and supernova
remnants stay proportional to the stars, and the star-type analysis in
`star-types-by-galactic-radius.md` (GEN.133) lists what could follow.
`layer_expected` (the progress bar) uses the layer-centre height factor only.
