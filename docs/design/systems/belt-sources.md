# Belt sources: the numbers the belt items need (GEN.1060)

Status: read 2026-10-30 (Sources 2). This is the research read of `belt-evolution.md` section 8 and `belt-origins.md`. Every number the belt items use is here with the paper it came from and a row in the stand-in registry (`gap-filling-plan.md` section 6.95, rows BR1 to BR18). The Solar System belt, the Kuiper belt, Ceres, Pluto and Eris are test targets, never inputs.

Evidence tags: READ = the paper (arXiv text) was opened and the number is in it; 2ND-HAND = the number comes from a paper that cites the original, which was not reachable; NOT REACHED = no copy found, nothing taken from it.

## 1. What was read, and what was not

| Source | Used for | Status |
|---|---|---|
| Vokrouhlicky et al. 2015, Asteroids IV (arXiv:1502.01249) | Yarkovsky drift rate and its 1/D law (stands for Bottke 2006) | READ; Bottke 2006 itself NOT REACHED |
| Mustill and Wyatt 2012 (arXiv:1110.1282) | chaotic zone half-width, Wisdom 1980 law and its numerical check (stands for Wisdom 1980) | READ; Wisdom 1980 and 1982 NOT REACHED (2ND-HAND) |
| Minton and Malhotra 2010 (arXiv:0909.3875) | loss history of asteroids in unstable resonances (stands for Moons and Morbidelli 1995) | READ; Moons and Morbidelli 1995 NOT REACHED |
| Walsh et al. 2011 (arXiv:1201.5177) | implantation and depletion of the belt in a Grand Tack | READ |
| Morbidelli et al. 2015, Asteroids IV (arXiv:1501.06204) | belt mass now and at the start, Grand Tack depletion, loss afterwards | READ |
| Raymond and Izidoro 2017 (arXiv:1707.01234) | planetesimals implanted into the belt by giant growth | READ; no single survival fraction is printed (figure 3 axes only) |
| Ida and Guillot 2016 (arXiv:1610.09643) | the snow line as a place where rock-rich planetesimals form | READ; qualitative, no number is used |
| Stern and Levison 2002 | the dynamical discriminant | 2ND-HAND from Soter 2006 (Lambda = k M^2 / P, a gap of five orders of magnitude) |
| Soter 2006 (arXiv:astro-ph/0608359) | mu = M / m, planet if above 100 | READ |
| Margot 2015 (arXiv:1507.06300) | the orbit-clearing mass and the discriminant Pi | READ |
| Lineweaver and Norman 2010 (arXiv:1004.1091) | the radius at which bodies turn round | READ |
| Kenyon and Bromley 2008 (arXiv:0807.1134) | the largest body a young icy disc grows | READ |
| Debes and Sigurdsson 2002 (arXiv:astro-ph/0202273) | planets that survive the giant branch, instability after mass loss | READ |
| Koester, Gaensicke and Farihi 2014 (arXiv:1404.2617) | white dwarfs polluted by debris, accretion rates (stands for Zuckerman 2010) | READ; Zuckerman 2010 NOT REACHED (Koester cites 20 to 30 percent from it) |
| Ivezic et al. 2001 (astro-ph/0105511) | main belt size law, 0.4 to 40 km; its indices 2.3 and 4 are differential, cumulative 1.3 and 3.0 (GEN.1061) | READ |
| Fraser et al. 2014 (arXiv:1401.2157) | Kuiper belt size law, hot and cold | READ |
| Bernstein et al. 2004 (arXiv:astro-ph/0308467) | Kuiper belt mass and its largest bodies | READ |
| Pitjeva and Pitjev 2018 (arXiv:1811.05191) | the masses of the main belt and the Kuiper belt | READ |
| Levison et al. 2008 (arXiv:0712.0553) | Kuiper belt delivered during a giant instability | READ; the disc cut at 30 to 35 AU is the Solar System's own, a test target |
| Sibthorpe et al. 2018 (arXiv:1803.00072), Marshall et al. 2014 (arXiv:1403.6186) | debris disc incidence round Sun-like stars | READ |
| Gladman and Volk 2021, Wisdom 1982, Moons and Morbidelli 1995, Zuckerman 2010, Bottke 2006 | | NOT REACHED: closed or not on arXiv; the numbers come from the papers named above |

## 2. The numbers

**Yarkovsky drift (BR1).** `da/dt` scales as 1/D for a large body (Vokrouhlicky 2015, section 2.1). The best measured case is Bennu: D = 493 m, da/dt = -18.95 +/- 0.10 e-4 au/My at a mean distance of 1.10 au; Golevka (D = 280 m) -4.52 e-4. So `da/dt x D = 9.3e-4 au km/My` at 1.1 au. Scaled to a belt with the flux (1/a^2, a stand-in that ignores the thermal-parameter dependence), a 1 km body at 2.5 au drifts 1.8e-4 au/My, the "order 1e-4" of the design note. A 40 km body drifts 0.05e-4 au/My, 0.2 au in 4.5 Gy: bodies above about 40 km do not move, as the design says.

**Chaotic zone (BR2).** The zone of overlapping first-order resonances round a planet of mass ratio `mu` is `delta a / a = 1.3 mu^(2/7)` (Wisdom 1980, read here in Mustill and Wyatt 2012, who confirm 1.35 mu^0.278 by numerical integration at zero eccentricity). Mean-motion resonances stay `a_res = a_giant (q/p)^(2/3)`.

**Loss in unstable resonances (BR3).** Dynamical chaos is the main loss for bodies of 10 km and more in the present belt. The loss history follows a logarithmic law, established by 1 My after the planets reached their orbits and holding to 4 Gy: about 50 percent of the large asteroids lost since then, equal amounts (about 30 percent of the present number) in 10 to 200 My and in 0.2 to 4 Gy (Minton and Malhotra 2010). The design's "lost in 1 to 10 My" is the resonance clearing at the start, not the later tail.

**Early depletion (BR4).** The primordial belt held about 1 Earth mass (within a factor of a few, from the disc that built the terrestrial planets), against 4.0e-4 now (Pitjeva), three orders less (Morbidelli 2015). The Grand Tack depletes it by "a factor of several hundred"; of the 1.3e-3 Earth masses of implanted S-type bodies at gas removal (efficiency 0.07 percent), 50 to 90 percent is lost in the instability and a further factor 2 to 3 by chaotic diffusion (Walsh 2011). The belt left by the first event is about 4 times the present one; half of that goes in the later instability and half by resonances (Morbidelli 2015). Raymond and Izidoro print no single fraction. For the generator: a Grand Tack type event leaves 2e-3 to 1e-2 of the belt (BR4); an instability alone removes 50 to 90 percent (BR5).

**Planet or not (BR6, BR7).** Soter: `mu = M / m` (m the mass of the other bodies in the zone), a body with mu above 100 is a planet; the Solar System has a gap, Earth 1.7e6 times Earth mass over the near-Earth objects, Ceres 1/3. Margot: the mass that clears an orbit in the main-sequence lifetime is `M_clear = 1.9e-4 C^(3/2) (a / 1 au)^(9/8) (M* / Msun)^(5/2)` Earth masses, with `C = 2 sqrt 3`, and `Pi = M / M_clear` above 1 is a planet (Jupiter 4.0e4, Mars 54; the dwarf planets are under 1). Stern and Levison's `Lambda = k M^2 / P` has a gap of five orders of magnitude between Mars and the largest asteroids and Kuiper belt objects. Margot's constant in the design note, 807, is `1 / (1.9e-4 C^(3/2))` = 816 with the printed digits.

**Round (BR8, BR9).** The potato radius, where a solid body turns spherical, is about 300 km for rock and about 200 km for ice (Lineweaver and Norman 2010): 3.4e20 kg for rock at 3 000 kg/m3 and 3.3e19 kg for ice at 1 000. The design's band for a fully relaxed hydrostatic body (ice 1e20 to 1e21 kg, rock 1e21 to 1e22) is from Bierson and Nimmo and Thomas, which were not reached for this table; it stays a stand-in band (BR9) and the potato radius is the lower, "roughly round" bound.

**The top of the size law (BR10, BR11).** A young icy disc grows planets of at most about 1 750 km radius, nearly independent of disc mass, stellar mass and age, holding at most 3 to 4 percent of the initial solids (Kenyon and Bromley 2008): a test of the top end of the order statistic. The Kuiper belt: cold classical bright-end slope alpha1 = 1.5 (+0.4 -0.2), break H_B = 6.9, hot alpha1 = 0.87 (+0.07 -0.2), H_B = 7.7 (Fraser 2014), both with a faint-end slope near 0.2; masses 0.01 (hot) and 3e-4 (cold) Earth masses. Bernstein 2004: the classical belt holds about 0.010 Earth masses, mostly in bodies near 100 km, its largest about 2 percent of its mass; the excited population's largest holds tens of percent and is about Pluto's mass. Pitjeva 2018: the main belt (4.008 +/- 0.029) e-4 Earth masses; the Kuiper belt with its 31 largest objects (1.97 +/- 0.03) e-2, the 31 largest being about 40 percent of it.

**White dwarfs (BR12 to BR14).** Planets inside about 5 au do not survive the giant branch; any beyond 5 au survive and widen as the star loses mass (a grows as M_initial / M_final). When mass is lost, two-planet systems that were marginal become unstable: of the unstable pairs about 8 percent collide, 92 percent do not, and of those about 40 percent settle with a planet on a closer, more eccentric orbit (Debes and Sigurdsson 2002); about 2 percent of young white dwarfs should show a disrupted-planet signature. Koester 2014 (cooling ages 20 to 200 My): 56 percent of 85 DA white dwarfs show heavy elements, at least 27 percent are accreting now (up to 50 percent), earlier surveys 20 to 30 percent (Zuckerman, 2ND-HAND); accretion rates a few 1e5 to a few 1e8 g/s.

**Debris discs round Sun-like stars (BR15).** 17.1 (+2.6 -2.3) percent of 275 F to K stars have an excess (Sibthorpe 2018), 20.2 +/- 2.0 in DUNES and 29 +/- 9 percent among exoplanet hosts (Marshall 2014). Most discs have a fractional luminosity near 1e-5 at blackbody radii of 7 to 40 au (about 40 au for real grains); the present Kuiper belt sits near the median, "half of stars have either depleted their Kuiper belts or had a lower planetesimal formation efficiency".

**Where belts form (BR16).** Ida and Guillot 2016: when the pebble flux exceeds a critical value (low viscosity, alpha under 1e-3), silicate dust piles up inside the snow line and rock-rich planetesimals form there. Qualitative; the snow line stays the generator's own `snow_line_au`.

**Size laws of the ledger (BR19 to BR21).** `physics/belt_ledger.py` (GEN.1061) holds two piecewise cumulative laws: the hot Kuiper belt's (q 1.0 and 4.35 about a break of 121 km, Fraser 2014, BR19) and the rocky belt's (q 1.3 below 5 km, 3.0 to 200 km, 4 above: Ivezic 2001 for the first two, a stand-in for the top, BR20). Ivezic's 2.3 and 4 are differential indices, so the cumulative exponents are 1.3 and 3.0 (corrected in `small-body-size-laws.md` 1.4). A stretch of q = 3 holds equal mass in every doubling, so the rocky law must steepen above the measured range; the break at 200 km makes 10 to 20 bodies above 1e19 kg at the main belt's mass (the Solar System has about 20), a count chosen to land there and so tier T9, tweak L3.

**Cut rule (BR17, BR18).** The cut values of Boss's rule (1e19 kg and 1e15 kg, 90 percent of the belt mass, 300 per belt) are stand-ins from the needed output, registered at tier T9.

## 3. The design's [R] tags

Every [R] in `belt-evolution.md` and `belt-origins.md` is covered: Ida and Guillot (BR16); resonance loss 1 to 10 My, Wisdom 1982, Moons and Morbidelli (BR2, BR3); Walsh, Raymond and Izidoro, Morbidelli 2010 depletion by 99 percent (BR4, BR5); Yarkovsky (BR1); YORP spin-up times (not used by the belt mass; no row); giant-branch ablation (no number used, the white dwarf rows BR12 to BR14 cover the survivors); Kenyon and Bromley (BR10); the 30 bodies above 1e19 kg and the KBO count (BR17 and the test targets above); round masses (BR8, BR9); the survival lognormal (BR4, BR5); the belts' masses (BR11, test targets).
