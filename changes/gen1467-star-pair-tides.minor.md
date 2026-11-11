**Headline:** The stars of a close pair now raise tides on each other: a star near enough to its partner is spun to the pair's orbital rate and its spin axis is turned toward the orbit normal, its activity follows its new spin, and the pair's apsidal precession (relativity, both tidal bulges, both rotational bulges) is noted in the effect ledger. A lone star and a wide pair are unchanged.

### Added

- `physics/star_pair_tides.py`: the synchronisation time of a star under its partner's tide (Zahn's convective equilibrium tide in the form of Hurley, Tout and Pols 2002 for a star with a convective envelope, Zahn's dynamical tide for a radiative one), the pseudo-synchronous rate of an eccentric orbit (Hut 1981), the spin and spin-axis alignment left after the star's age, and the apsidal precession of a pair's orbit by cause. White dwarfs, neutron stars and black holes keep the spin they were drawn with.
- `StarSystem._tidal_spin_of_stars` applies it to every lit star of the tree before the planets and the close pair's merged star are made: a changed star's rotation period and spin axis are set, its activity (X-ray, XUV, flares, dipole) is made again from a stream of the star's own so no draw of the run moves, and the two stars of a pair node get effect-6 (precession) ledger rows.
- Stand-in rows KT1 to KT8 (GEN.725, section 6.263).

### Changed

- The rotation period, spin axis and activity of a star within a few weeks of its partner differ from the independently drawn ones. A lone star, a wide pair and any compact star are as before; the 7-day circular cut of the pair's orbit stays.

generation-output: changed
generation-affects: fill
