# The stellar mass function (GEN.927 to GEN.931)

**Status (2026-10-30):** built on `bug1/star-epoch` for the long-term star epoch (GEN.932, parent GEN.922). Boss approved the
full, no-patch fix at 5:27 PM EDT 10-29 ("let's build a robust fix so we match the curve segments for all mass levels and
star-like-object-types across the mass range").

Every star the generator makes, and every remnant, comes from one family of curves. This note says what the curves are,
what is set and what is derived, and what each consumer reads.

## The curves

| Curve | Name in code | What it counts | Where it comes from |
|---|---|---|---|
| Single stars | `imf.SINGLE` (`xi`) | Every star born, by mass | Set: [R0705] (PASP 115:763, Table 2, disk, individual stars) |
| Primaries | `imf.PRIMARY` (`P`) | The star a system is built around, by its birth mass | Derived: `xi` less the companions of heavier primaries (`generation/primary_mass_function.py`) |

**Single stars.** A log-normal in log mass below 1 solar mass (centre 0.079, width 0.69 dex) and the power law
`dN/dM ~ M^-2.3` above, joined continuously at 1 solar mass. The function is kept as a continuous piecewise power law
(`tuning.IMF_BREAKS_SOL`, `IMF_SLOPES`): the log-normal is cut at every 0.05 dex and each pair of nodes is joined by the
chord of `ln dN/dM` against `ln M`, within 0.07% of the log-normal. The sampler, the integrals and the number counts of
`physics/imf.py` stay exact on it. The birth slope is 2.3; the steeper 2.7 of the field [R2499] counts
stars still alive, and the disk history reproduces it from this slope and the ages. The old two-segment function (1.3 below
0.5 solar masses, 2.3 above) did not have the log-normal's turnover, which is what the published counts of low-mass stars
show.

**Primaries.** A system's companions are stars too, and each is lighter than its primary. So at each mass the stars born are
the primaries of that mass plus the companions of the heavier primaries that fall there:

    xi(m) = P(m) + sum over primaries M >= m of P(M) n_c(M) (share of M's companions that land at m)

`n_c(M)` is the mean count less one of the `star-count` table for the class of the star (`multistar.companion_class`), and a
companion's mass is its primary's times a ratio uniform between `max(0.1, 0.08 / M)` and 1 (`multistar.companion_ratio_range`).
Companions are never heavier than the primary, so the equation is solved from the top mass down in one pass, with no
iteration. The result peaks higher than `xi`, as the published systems' functions do (Chabrier's own for systems has its centre
at 0.22). Systems per star born is 0.697 and systems per solar mass of stars born 1.128 (`systems_per_star`,
`systems_per_solar_mass`).

What the derivation leaves out, and the reason: the companions past the seventh star of a system (a chance of about 1e-5 of
a system); the brown dwarf companion (below the hydrogen-burning limit, outside `xi`); the mergers and stripped stars of
`multistar`, which move a star between bins after its birth; and a twin excess in the mass ratio [R3108],
which the draw does not make either.

## Who reads which

| Consumer | Curve | Why |
|---|---|---|
| The star a placed system is built around (`stellar_evolution.sample_living_star`, `star_population._mass_grid`, `system_mass._placed_bins`) | `PRIMARY` | A system has one primary; its other stars come from its own tree |
| Systems per solar mass, the galaxy's system count (GEN.929) | `PRIMARY` over `SINGLE` | Both are on one scale |
| Per-star rates: H II regions, black hole and neutron star rates, mass return, nebula seeds, remnant progenitors | `SINGLE` | They are rates of events per star born, whatever system the star is in |

## The remnants of the same curve

* **White dwarfs** carry one initial-final mass relation, `tuning.WD_IFMR_PIECES` ([R0961]:21, MIST
  fit, equations 4 to 6, read): `0.080 M + 0.489` for 0.83 to 2.85 solar masses, `0.187 M + 0.184` for 2.85 to 3.6, `0.107 M
  + 0.471` for 3.6 to 7.2, the last piece carried to 8 (a 1.33 solar mass white dwarf at the top). Read by the star model, the
  member draw, the disk history and `mass_loss`, through `physics/remnant_mass.py`. The [R2230] relation the
  generator used before ran 0.04 to 0.06 under it, and another relation, the disk history's, 0.1 to 0.16 under at 5 to 7
  solar masses. The disk history's return fraction is 0.384 (it was 0.371).
* **Neutron stars** are drawn from the Gaussian of the slow pulsars ([R3396]: 1.49 +- 0.19), held to 1.1 to 2.2
  solar masses; `remnant_mass.share_above` and `moments` are the same Gaussian's for the scatter cut and the placement terms.
* **Black holes** stay uniform over 5 to 20 solar masses: the widest of the published ranges stands.

## What moves, and why it is not a defect

The massive share (8 to 150 solar masses) of the stars born goes from 0.637% to 0.713%. The bound neutron stars are fixed by
the census (0.004 per star), so the share of the massive stars left to fail and make black holes is larger: the middle draw of
the solar failed share goes from 13.6% to 17.5% (`black_hole_formation.share_cap`). A galaxy with the Milky Way's 5e10 solar
masses of stars holds 1.0e8 black holes, in line with the published estimate, and a metal-rich busy shell holds 1.6 neutron
stars to a black hole (it was 2.1).

## Source tiers

The four numbers of the log-normal are Confident Unverified, level L0 (recalled, not read in the paper, and used as
published); the Cummings pieces and the Ozel and Freire Gaussian are Verified. The primary function has no set number.
Registry section 6.50 of `research/gap-filling-plan.md`.

## Open

* GEN.931: the arm contrast by population and the A stars.
* GEN.936 and GEN.948 (built on their own branch): the main sequence from 0.07 to 0.1 solar masses and the brown dwarf cooling track are read from the [R0254] table, and the star/brown-dwarf boundary is one split at 0.08 with a class that depends on the place's [Fe/H].
* GEN.951: the remaining items of the chain.
