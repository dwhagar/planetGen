# Where ejected rogue planets drift, and where each one is stored

Plan (research written 5:20 PM EDT 2026-10-30 as the lane reported; nothing built; two questions for Boss are open in section 10). Tags: [REPO] read in the repo, [C] computed here (script `/mnt/project-files/data/rogue-drift/hill_drift.py`), [S] from a source, [SI] stand-in to register in GEN.725.

## 1. Short answers

1. **Drift is large.** A planet ejected with a kick of about 3 km/s relative to its star stays within about 0.7 kpc of the star's galactic radius (99 percent, [C]) but spreads along the orbit by roughly 1 kpc per 0.5 Gyr and by tens of kpc after a few Gyr. With ejection times uniform over 0.05 to 10 Gyr ago, the median rogue is about 7 kpc (along the orbit) from where it was ejected, and only 1.5 percent are within 100 pc [C]. A sector is 4 pc [REPO `DEFAULT_SECTOR_EDGE_PC`], so essentially no rogue is anywhere near a sector that holds its donor star.
2. **Boss's normal-vector idea, taken literally, gives the wrong drift.** Total velocity is about 220 km/s along the orbit, so a kick perpendicular to it has no along-orbit part. That is the one component that causes the long drift (section 3). Result in the sim: rogues stay within about 100 to 300 pc of the donor forever (table in 3.2). The fix is to apply the kick in the star's own frame (isotropic, or isotropic with a mild in-plane bias) and then add the system's motion afterward. The rest of the idea (pick an ejection time, run time forward from that point) is exactly right and has an exact closed form.
3. **Recommendation: option D, "virtual donors"** (section 5). The sector builds its own rogues from donors it does not need to own: each rogue draws an ejection event, runs the closed-form drift, and is placed in this sector by construction. Deterministic, sector-independent, cost per sector roughly the cost of drawing 37 ejection events (section 6). Real donors are never enumerated (that costs 10^7 to 10^9 systems per sector at the median age, section 6).
4. **The ID needs no change.** It is built from the rogue's present galactic position [REPO object-ids.md], which is the position this design computes.

## 2. How ejection happens, and the speeds

Three channels, all give the same shape: a speed at infinity of order the escape speed at the scattering orbit or lower, random direction.

* **Planet-planet scattering** (the architecture's GEN.728 path). Ejection is energy-limited by the perturber: v_inf^2 is about v_esc(a)^2 times a factor below 1, so tens of km/s at the very inner orbits and 1 to few km/s at 5 to 30 au. [R2031] give a Jacobi-energy bound on the ejection velocity and a minimum ejectable mass; the abstract gives no km/s numbers [S, abstract only]. [R2470] run star plus giant plus light planet N-body: the ejected body's mass has little effect over four decades, the giant's mass sets scale and timing (t_ejec scales as 1/M_J), eccentricity stretches only the tail: ceiling about 23 km/s, rare kicks near 80 km/s [S, abstract only].
* **Galactic population kicks**: [R3937] use a phenomenological kick with median 0.7 to 10.7 km/s across five mass bins, and integrate 10^6 particles for 0.1 Gyr: median displacement from the no-kick control 0.7 down to 0.054 kpc; a few percent leave the disc [S, abstract only]. This matches the sim here at 0.1 Gyr (median 270 pc for 3 km/s).
* **Passing stars and cluster encounters** eject the outer planets of a system, with speeds of order the encounter speed ratio times the orbital speed ([R2056, R2837] for clusters), usually 1 km/s or less at 50 au and up. Not read this session [listed, not verified].
* **Later speed-up of floaters by scattering off field stars**: [R2575]: equilibrium velocity about sigma times sqrt(2 ln(m/m_p)) [S, abstract]. This is secular heating and is ignored here (flagged, section 8).

**Replaced (GEN.1113, 2026-10-30):** the kick speed is now `v_c X` from the scatterer's mass and distance ([R2031: Eq. 26], checked against [R2470]), mean 1.15 km/s for a Jupiter at 5.2 AU, ejected-mass independent; see `docs/design/research/rogue-ejection-kick.md`. The earlier stand-in follows for the record: kick speed log-normal, median 3 km/s, sigma 0.8 in ln (so 16 to 84 percent is 1.3 to 6.7 km/s), tail to about 25 km/s. Tier: Marginal educated guess from physics, tweak L2 (tuned to the published range 0.7 to 10.7). R1 replaces it with the speed from each ejection's scattering orbit (escape speed at the scattering orbit times the energy factor).

**Direction.** Scattering is near-isotropic in the star's frame with a lean toward the orbital plane (a planet ejected from a flat system leaves within the system's inclination spread). Direction relative to the *galactic* motion has no physical reason to be perpendicular; the galactic motion does not enter the scattering. [reasoning; no source gives an anisotropy figure, so isotropic with a flag]. Stand-in [SI]: isotropic, tier Confident educated guess from physics, tweak L0.

## 3. Drift in the galaxy

### 3.1 Closed form (Hill's equations)
For two bodies on nearly the same galactic orbit, the separation obeys the linear epicyclic equations (Binney and Tremaine, Galactic Dynamics, 2nd ed., section 3.2; the Hill equations for tidal streams):

* x'' = 2 Om y' + 4 Om A x, y'' = -2 Om x', z'' = -nu^2 z, with Om = A - B, kappa^2 = -4B(A - B).
* Oort constants A = 15.3 and B = -11.9 km/s/kpc ([R0497] gives 15.3, -11.9), so Om = 27.2 and kappa = 36.0 km/s/kpc; vertical frequency nu about 76 km/s/kpc (Gaia-era midplane density 0.1 Msun/pc^3, nu = sqrt(4 pi G rho)) [S for A, B; nu is a standard estimate, [SI] tweak L1]. In the repo these should come from docs/design/galaxy/galactic-potential.md; the numbers here are only to size the effect.
* A kick (u, v, w) relative to the star at t = 0 gives
  * x = (2 Om v / kappa^2)(1 - cos kappa t) + (u/kappa) sin kappa t
  * y = -(4 Om A / kappa^2) v t + (4 Om^2 v / kappa^3) sin kappa t - (2 Om u / kappa^2)(1 - cos kappa t)
  * z = (w/nu) sin nu t

Meaning: u (radial kick) and w (vertical) only make bounded oscillations (radial amplitude u/kappa, vertical amplitude w/nu: 80 pc and 40 pc for 3 km/s). **v, the along-orbit kick, changes the guiding radius and the orbital period, so the rogue slides along the orbit at 1.285 v**: 3.9 pc/Myr for 3 km/s.

### 3.2 Simulated displacement from the donor [C] (200,000 draws, kick median 3 km/s, sigma 0.8)

| Time since ejection | isotropic: r 16/50/84 pc | in-plane (z x 0.2): median pc | Boss's normal-to-motion: r 16/50/84 pc |
|---|---|---|---|
| 10 Myr | 13 / 30 / 66 | 31 | 13 / 30 / 65 |
| 100 Myr | 98 / 274 / 716 | 350 | 53 / 144 / 365 |
| 1 Gyr | 477 / 1,746 / 4,867 | 2,285 | 23 / 55 / 131 |
| 5 Gyr | 2,263 / 8,331 / 23,219 | 10,898 | 39 / 105 / 266 |

Median vertical offset stays 12 to 17 pc (isotropic), 3 pc (in-plane). Radial excursion: |x| below 0.72 kpc for 99 percent.

**Limits of the closed form.** It is linear, valid to about 1 to 2 kpc of separation. Beyond that the along-orbit offset is an arc of the orbit (the galactic year at the Sun is about 230 Myr), and at 5 Gyr the 84th percentile (23 kpc) is more than half the orbit circumference (about 50 kpc), so the azimuth is fully mixed. Rule: compute the azimuth shift as an angle, `dphi = y / R_guide`, wrapped modulo 2 pi; the guiding radius shift is `2 Om v / kappa^2` from the line above; keep the closed form for x and z. Both pieces are bounded, so no integrator is needed. If a real potential (bar, arms) is wanted later, use the orbit integrator behind galactic-potential.md instead; same interface.

### 3.3 What this means for the rogue population
* Rogues older than about 1 Gyr have lost the donor's azimuth, so their density depends on galactocentric radius and height only: a smooth band (radial scale 0.1 to 0.7 kpc, vertical about 20 to 50 pc) around the star density, **with spiral arms, clusters and the bar smeared out**. Only rogues younger than about 100 Myr (a few percent) still follow the donor's structure.
* Stars and rogues are therefore not co-located at the sector scale, which is the physical answer to "why is there no rogue next to its star".
* Rogue density is nearly the azimuthal average of the star density; that is cheap to compute (section 5).

## 4. Options for where a rogue lives

| | A: donor's sector, drawn inside it | B: true drifted position | C: B with bounded enumeration | D: virtual donors (recommended) |
|---|---|---|---|---|
| Position correct | no: pinned to the donor, inconsistent with a multi-Gyr history | yes | only inside the cap radius | yes, by construction |
| Sector independence | yes | needs the donors of the whole galaxy | donors within the cap | yes, no neighbour sectors read |
| Cost per sector | about 1 x its donors' ejections | 10^7 to 10^9 donors | up to 10^4 to 10^7 donors | about 37 ejection draws |
| History sentence | false ("ejected 3 Gyr ago" yet sits at the donor) | true | true within cap | true |
| Donor system exists on the map | yes | yes | yes within cap | no (record only; regenerated on demand) |

Why A fails beyond looks: it puts rogue count exactly proportional to the sector's own star count, which is wrong by the clustering of the star field (a sector with 12 stars would get twice the rogues of one with 6, when real rogue counts hardly depend on which star is next door). Fraction of rogues whose donor shares or touches their sector: 4.5e-5 within 4 pc, 2.8e-4 within 10 pc [C], so A gets the true answer for essentially none.

Why B and C fail on cost. Donors within radius R at 0.1 systems/pc^3 [SI, solar-neighbourhood order; the repo's density function replaces it]: R = 30 pc 1.1e4; 100 pc 4.2e5; 300 pc 1.1e7; 1 kpc 4.2e8; 2 kpc 3.4e9 [C]. The median rogue needs the 7 kpc radius. Even testing donors only for a 100 pc cap catches 1.5 percent of rogues for 4.2e5 star rolls per sector.

## 5. Option D: how it works

Principle: the closed form is *translation invariant*. Where the rogue lands relative to its donor depends only on (ejection time, kick, galactic position). So instead of asking "which donors land here", the sector draws ejections and **defines** the donor's position as the sector position minus the displacement.

Per sector S (seed from the sector, as all draws, via `draw.py`):
1. **Count.** Poisson mean lambda_S = yield x V_S x n_star_smooth(S). `yield` is rogues per star from R1 (stand-in 5.8, [REPO `ROGUE_PLANET_RATE_PER_STAR`]); `n_star_smooth` is the star density smoothed along the orbit: for old ejections the azimuthal average at the sector's galactocentric radius R and height z, plus a radial smoothing of order 0.5 kpc; for young ones the density at the sector. Cost: evaluate the existing density function at about 16 points on the ring (azimuth steps), once per sector, then weight by the age mix. Typical lambda_S is 5.8 x 64 pc^3 x 0.1 = 37 rogues for a 4 pc sector at solar density [C arithmetic].
2. **Draw each ejection event** j from the yield table (R1): donor class (star mass, metallicity, age drawn from the population model at the sector's position), ejection time t_j (after formation, so t_j can not exceed donor age), rogue mass, kept moons, kick speed v_j and direction.
3. **Run the closed form** (3.1/3.2) for t_j and the kick: this gives the displacement Delta_j, the rogue's present velocity (host peculiar velocity + orbital motion at ejection + kick, evolved: the same solution differentiated) and the donor's ejection-time position x_S - Delta_j.
4. **Thin** by n_star(x_S - Delta_j) / n_star_ref so that an ejection that comes from a poor neighbourhood is less likely. Cheap, and removes the bias of drawing donors at x_S.
5. **Place** the rogue uniformly inside the sector. Its position (hence its ID, `objectId.format_id`) is exact, nothing is displaced.
6. **Store the donor as a record**, not a generated system: donor class, donor seed, donor distance at ejection (|Delta_j|), time since ejection, ejection speed, ejection orbit slot. When someone opens the rogue, the donor system is *regenerated from the seed* by the same architecture code (it is virtual, never stored), and the rogue's composition, moons and heat history come from that regeneration.

Properties: deterministic (every number is a function of the sector seed and event index), independent (reads only the density function, never a neighbour sector), no real star is claimed to be the donor. The history sentence: "Ejected from a K dwarf 2.4 Gyr ago at 3.2 km/s; that star lies about 4 kpc along the Galactic orbit away and is not on the map."

Caveats registered, not hidden:
* the donor properties are conditioned on the event rather than on the star population around it. The thinning in step 4 handles density only; a metallicity gradient across the displacement (0.07 dex per kpc, [SI] order) is ignored unless R1 shows it matters.
* the smoothed density assumes a stationary disc (ejection-time density equals present density). Radial migration and the ageing of the disc break this at the 10 to 30 percent level [SI].
* real nearby donors (the under-1-percent within one sector) never get a link; a later optional item can match recent ejections (t under about 2 Myr) to actual neighbour systems.

## 6. Cost

* Option D per sector: 16 density evaluations plus about 37 event draws plus 37 closed-form evaluations (a few trig calls each), then the same body generation the rogue gets today. Added cost is the donor regeneration only when a rogue page is opened or when the rogue's heat history needs the donor (first rogue body generation). Measure in R1: time of a full donor architecture roll against one rogue body.
* Option A (the cheap reference): same count of rogues, but each from a real donor already built; no extra donors needed.
* Option B/C: see the table in section 4; infeasible.
* Storage: one donor record (class, seed, three numbers) per rogue, about 40 bytes.

## 7. Tie to the map and ID

* Position and ID unchanged in form. Rogues are drawn in their own sector like today, only the count and the histories change. Collision number per position (`_claim_object_ids`) unaffected.
* Zoom levels: because rogue density is smooth on kpc scales for old rogues, the low zoom levels can use `n_rogue(R, z)` directly (no per-sector draw), and the same function is the Poisson mean at high zoom. One function, two uses, so the levels agree. Needs the Z0-Z5 plan check by the map lane.
* The donor never appears on the map; a "parent" link on the rogue's page shows the record (class, age, distance, speed), not a map jump. If Boss wants a jump, point to the donor's regenerated system page by seed instead of a position.

## 8. Open points and stand-ins to register in GEN.725

| Row | Value | Tier | Tweak |
|---|---|---|---|
| Kick speed: `v_c X` (rows RK1, RK3, RK4 of 6.87), moons RK5 to RK8 | read (GEN.1113) | Verified (RK1, RK3, RK5, RK8), marginal (RK4) | L0 to L3 |
| Kick direction isotropic | [SI] | Confident educated guess from physics | L0 |
| Oort A, B and vertical frequency nu | A 15.3, B -11.9, nu 76 | Verified (A, B), Confident Unverified (nu) | L0 / L1; replace with galactic-potential.md values |
| Donor density 0.1 systems/pc^3 for the cost table | [SI] | Marginal Unverified | L1, only sizes the cost |
| Ejection time after formation: R1 gives it; test uses uniform over donor age | [SI] | Unverified | L3 |
| Secular heating of floaters (Lau and Lai) ignored | omitted | gap | document only |

Read in full since (GEN.1113): 2508.18239, 2608.02945, 2609.21027 and [R4533]; to read still, [R0497, R2056, R2837] for the Oort constants.

## 9. Items

| Item | What | Owner (best guess, can move) |
|---|---|---|
| GEN.1113 | D1 kick research | Planets 1 |
| GEN.1114 | D2 rogue_drift.py | Aux 7 |
| GEN.1115 | D3 smoothed rogue density | Sources 2 |
| GEN.1111 | D4 virtual donors (held for Boss) | Planets 1, train 19 |
| GEN.1116 | D5 regenerate donor on demand | Planets 1 |
| GEN.1115 (note) | D6 low zoom uses the smooth density | Core 2 |
| GEN.1117 | D7 register the stand-ins | with D1 and D2 |

## 10. Questions for Boss
1. Rogues placed by option D (virtual donors, not on the map): yes (recommended), or keep real donors only (option A, wrong positions)?
2. Kick direction: isotropic in the star's frame (recommended), or still perpendicular to the galactic motion (keeps rogues within a few hundred pc of the donor)?
