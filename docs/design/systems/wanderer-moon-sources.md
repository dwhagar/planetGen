# GEN.560 (moon parts): ejection survival of moons and circumplanetary disks of isolated planetary-mass objects

**Status:** research, no behaviour change. It closes GEN.560 with the black hole part ([black-hole-formation-sources.md](black-hole-formation-sources.md)) and the interstellar size law ([small-body-size-laws.md](../composition/small-body-size-laws.md) section 1) already written. It feeds GEN.558 (moons of wanderers: ejection survivors and isolated-disk moons as table rows) and the moon rule of [unified-star-like-objects.md](../architecture/unified-star-like-objects.md) section 3.6.

Method: each paper was downloaded as the arXiv PDF and read in the text. "Quoted" is copied from the paper; "own reading" is a plain restatement. Source states as in the small-body note: READ is a number found in the paper this session; 2ND-HAND is a number a read paper quotes from another; GAP is a number no paper read here gives. Raw papers stay out of the repository. The plan's two recalled sources (Debes and Sigurdsson 2007, "Hong et al. 2018") are checked in section 3: the first exists as recalled, the second is not the paper that has the numbers.

## 1. Moons that survive the ejection of their planet

| Case | Number | Source | State |
|---|---|---|---|
| Earth-mass planet with a lunar-mass moon at 30 Earth radii, scattered by a Jovian planet; 2700 runs of 10^6 yr | 123 runs (4.6%) end with the Earth-Moon pair ejected still bound; 90 (3.3%) eject the Earth alone | Debes and Sigurdsson 2007, ApJ (arXiv:0709.0945), section 3 | READ |
| Same | median eccentricity 0.13 and semi-major axis 30 Earth radii of the ejected pairs; tidal heat up to 4.2e22 erg/s (100 times radiogenic) for a rocky planet | same | READ |
| Same | free-floating binary planets: 7e8 in the Galaxy for 1e10 stars, a space density of 9e-4 per pc^3 (their assumptions) | same | READ |
| Three Jupiter-mass planets, innermost at 3 AU, 100 test-particle moons per ejected planet | about 85% of moons at the Galilean distances survive; moons survive out to beyond 200 planetary radii (0.1 AU in that setup); just over 25% of runs keep all or most of their moons | Rabago and Steffen 2018, MNRAS (arXiv:1809.05639), section 3 | READ |
| Same | survival fraction against the moon's initial distance: constant at small distances, then a power law; the least-squares index is -0.607 +- 0.017 (an index of -1/2 fixed gives a worse fit) | same | READ |
| Planet-planet scattering of two or three planets (circular restricted three-body plus N-body) | the probability that the minimum planet-planet distance exceeds 0.1 R_Hill is 48% and exceeds 0.01 R_Hill is 82% (fiducial runs, m1/M about 1e-2, m1/m2 = 10); regular Solar System moons orbit at about 1e-2 R_Hill, so a significant fraction survive; a moon inside 1e-2 R_Hill survives with at least 70% for any initial separation | Huang and Lai 2026, ApJ (arXiv:2508.18239v2), section 4 and conclusions | READ |
| Same | survival decreases somewhat as the planet mass and the mass ratio m1/m2 fall; the finite size of the planets does not change the statistic while the planet radius is far below the minimum distance | same | READ |

Own reading for the draw: the survival of a moon is a function of its distance to its planet in Hill radii at ejection, flat at about 0.7 to 0.85 inside 1e-2 R_Hill and falling as a power law of index about -0.6 outside it (Rabago and Steffen), reaching about 0.48 at 0.1 R_Hill (Huang and Lai). The two papers agree inside 1e-2 R_Hill (0.82 and 0.85). They differ in scale on the tail because one fixes the planet's distance from its star and the other the planet separation in Hill units; the wider range wins (the standing rule for two sources in question). Debes and Sigurdsson's 4.6% is a different quantity: the share of scatterings of an Earth-mass pair that eject the pair intact, for one set of start distances, not a survival fraction of moons that start bound.

## 2. Circumplanetary disks of isolated planetary-mass objects

| Case | Number | Source | State |
|---|---|---|---|
| Disk fraction of free-floating planetary-mass objects (M9 or later) in IC 348, Spitzer 3.6 and 4.5 um stacked | six with disks plus three uncertain of the sample: 46 +- 13/12 percent; "comparable to more massive brown dwarfs"; the objects keep disks for several million years | Seo and Scholz 2025, MNRAS (arXiv:2501.15930), abstract | READ |
| Same in NGC 1333 | 5 of 12 planetary-mass objects show excess: 42%, large uncertainty; the disk fraction does not decline substantially with falling mass in the substellar range | Scholz et al. 2023 (arXiv:2303.12451), abstract | READ |
| OTS 44, about 12 Jupiter masses, ALMA band 6 | unresolved detection; dust mass 0.07 to 0.63 Earth masses by the empirical correlations; one scaling of dust mass with central mass holds from brown dwarfs to planetary masses | Bayo et al. 2017, ApJ Letters (arXiv:1705.06378), abstract | READ |
| 2MASS J04414489+2301513 (about 0.05 Msun) | disk radius 15 to 30 AU, first resolved substellar disk | Ricci et al. 2013, quoted by Bayo et al. 2017 | 2ND-HAND |
| Number of free-floating planets known with a disk (to 2017) | four, listed with their masses as Cha 1109-7734 (8 Mjup), OTS 44 (12), J02265658-5327032 (13) and Proplyd 133-353 (at most 13, Fang et al. 2016) | Bayo et al. 2017, introduction | READ |

Own reading for the draw: a planetary-mass object formed on its own keeps a disk with probability near 0.4 to 0.5 for a few million years, so for the moon rule the draw can use "disk present" at about 45% for the youngest bins and none after about 10 Myr; the disk dust mass is about 0.1 to 0.6 Earth masses at 12 Jupiter masses, a Moon- to Mars-mass reservoir, so the moons it can make are few and small.

## 3. Checks on the recalled sources and gaps

- **Debes and Sigurdsson 2007:** exists as recalled (arXiv:0709.0945) and gives the numbers above.
- **"Hong et al. 2018":** the arXiv search for ejection and moon survival returns Rabago and Steffen 2018 (arXiv:1809.05639), which has the survival numbers, and Huang and Lai (arXiv:2508.18239, 2025/2026), not a Hong et al. paper. The planned reference is replaced by these two; no number is attributed to Hong.
- **GAP, the disk size limit of an isolated object:** the plan text says the zone comes from the object's rotation and its tidal radius. No paper read here gives a radius for a disk around a free-floating planet except the second-hand 15 to 30 AU of a 0.05 solar mass brown dwarf and the unresolved OTS 44. For a bound planet the usual truncation is a fraction of the Hill radius (Martin and Lubow 2011 are the sources to read); not read, so it is not used. Until it is read, the disk radius of an isolated object stays a tunable stand-in in the gap plan registry, marked STAND-IN.
- **GAP, tidal ejection speed:** Huang and Lai give the ejection velocity scaling with planet mass and initial semi-major axis (their sections 3 and 4); the figures were not transcribed here because GEN.557 has not asked for them.
