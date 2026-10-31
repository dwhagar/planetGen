# Where belts come from, a scatter phase, and how many Class D bodies a system needs

Status: plan (Boss asked 4:19 PM EDT 2026-10-30 and confirmed the three defaults of section 9 at 4:37 PM EDT; nothing built). No generator code changed. Asked for by Boss: (1) plan the items that make belts evolve and record on each belt where we think it came from (we know, because we made it); (2) where asteroid belts come from and what makes them; (3) whether systems need more Class D bodies given how many Sol has; (4) whether system generation needs an asteroid / dwarf planet scatter phase. Builds on `belt-evolution-and-body-ledger-2026-10-30.md` (the "ledger draft") and `reports/belt-root-cause-2026-10-30.md`.

Evidence tags: [S] read this session (arXiv abstract opened, cited); [REPO] read in the repo or measured from it; [R] recalled, not opened: the building lane reads the source before using the number, and registers a stand-in in GEN.725 (tier and tweak level) if it still cannot. The Solar System is a test target, never an input.

## 1. Short answers

1. **Origin is a stored fact.** Every belt (and outer disc) gets an `origin` cause and a `formation` record written by the generator when it is made, shown on the belt page as one sentence ("the leftover of a planet that never grew here, stirred by the giant at 5.2 AU"). The cause list is section 3.
2. **What makes a belt** is one thing in physics: a planetesimal population that did not finish growing into a planet, or was ground or scattered after. Five causes cover it (section 3). Today the generator has none of them; a belt is a 10 percent roll on a light slot.
3. **Class D:** the generator makes almost no dwarf-planet-size bodies, and it cannot make the outer ones at all. This is a real gap, but the cure is **not more Class D planets on the planet list**. It is the ledger's L1 tier of each belt and outer disc, where Sol's hundreds of 100 km to 1,000 km bodies live (section 5).
4. **A scatter phase belongs in generation**, as a step after the giants are placed and before the belts are finalised (section 6). It is the same step that tells a belt which part survives, which part becomes a scattered disc, and where the Trojans, Centaurs and resonant swarms sit.

## 2. What Class D is (my reading, stated because the docs disagree)

* **Today's definition [REPO]:** `physics/data/planet_classes.json`, row D "Small icy body": surface pressure under 10 Pa **and** water or ice at least 1 percent of the mass. Size is not in the rule. The class is derived from the finished body (epoch 121), so a Class D body is whatever small airless icy world the generator builds.
* **Older docs [REPO]:** `mass-object-class.md` and `small-bodies-composition.md` say Class D is "radius 50 to 500 km, so 100 to 1,000 km across" and "frozen hydrocarbons and ice". `tuning.PLANET_CLASSES["D"]` still carries that wording. Treat the radius as the **intended scale** and the conditions as the **rule**.
* **What can be Class D in generation [REPO]:** a planet must be at least `PLANET_MIN_MASS_EARTH = 0.1` Earth masses (Mercury is 0.055), and beyond the snow line at least `PLANET_MIN_MASS_OUTER_EARTH = 0.6016`. Lighter oligarchs are "the debris of a belt". The whole 50 to 500 km range is 1e-6 to 1e-3 Earth masses. **So no planet-list entry can be a 50 to 500 km body.** Class D appears only on (a) moons, and (b) planets of 0.1 Earth masses and up that are airless and at least 1 percent water, which is rare and is not what the class name means.
* **My reading for the redesign:** keep the conditions rule, and add that the class applies to any ledger body (belt member, dwarf planet, moon) with the same conditions. A 1,000 km icy dwarf planet and a 100 km icy moon are both D. Bodies below the dwarf-planet mass line (about 1e19 kg, section 5.2) are "minor bodies" with no class letter. Boss to confirm (section 9, question 1).

## 3. What makes a belt: the origin causes

Each cause is a rule the generator can apply and a sentence it can print. A belt can have several; the stored `origin` is the dominant one, `history` lists the rest.

| Code | Cause | Physics | Evidence |
|---|---|---|---|
| `failed-growth` | **A planet that did not finish.** Solids stall because the next oligarch would exceed the isolation mass, or a giant's growth starved the zone. | The inner disc forms planetesimals; where growth is stopped, the mass stays as a population. The planet architecture (GEN.728) already holds the chain; a gap or unfinished tail in it is a belt. | Planetesimals formed big, around 100 km and up, then the belt's size distribution has a bump near 100 km [S: Morbidelli et al., "Asteroids Were Born Big", arXiv:0907.2512; Weidenschilling, arXiv:1004.0270 gives the same 100 km peak]. |
| `giant-stirring` | **A giant excites it.** Mean-motion resonances and the chaotic zone around a giant raise eccentricities until bodies collide with each other (grinding) or hit a planet or the star. Growth never completes, and the zone is cleared of what drifts into resonance. | Chaotic zone half-width about 1.3 a (m/M*)^(2/7) [R: Wisdom 1980]; resonance positions a (q/p)^(2/3) [S physics]; ledger draft section 2.3. | Main belt gaps at 3:1, 5:2, 2:1 [R]. |
| `instability-depletion` | **A giant instability or migration swept it.** Jupiter's Grand Tack or a Nice-type scatter removes most of a belt and implants bodies from elsewhere. | Survival of a few percent; implants of outer-disc bodies into the inner belt. | Walsh et al. 2011 [R], Levison and Morbidelli, who deliver KBOs during the instability, arXiv:0712.0553 [S]. |
| `outer-disc` | **The cold leftover beyond the last giant.** Solids too sparse and orbital times too long to build a planet. A Kuiper-type belt. Self-stirring by its own largest bodies (Pluto-size) sets the grinding. | Kenyon and Bromley: KBOs form in a massive disc, then substantial collisional grinding and dynamical ejection [S: arXiv:0704.0259]. Mass today in Sol: 0.0197 +/- 0.003 Earth masses [S: arXiv:1810.09771]; 0.010 Earth masses for the classical part [S: Bernstein et al., arXiv:astro-ph/0308467]. | Herschel DEBRIS: 17.1 percent of F-K stars show cold discs, most at 7 to 40 AU blackbody radius; the data fit "all stars born with a narrow belt that then grinds", with Sol's own belt near the population median, so half of stars have depleted theirs or formed less [S: arXiv:1803.00072]. |
| `collisional-debris` | **The rubble of a destroyed planet** or a giant collision, one composition family, differentiated shards. Also transient hot dust. | A rare event, not a standing population. | DEBRIS: hot dust under 4 AU "likely transient events" [S: arXiv:1803.00072]. |
| `captured` | A belt around a pulsar or white dwarf from a disrupted disc or scattered survivors (already used for pulsar planets). | Mass-loss destabilisation, ledger draft 2.5. | Existing repo rows. |

What makes it **go away or change** over time, already in the ledger draft: collisional grinding (Wyatt 2007, `belt_mass.py`, built but not called), Yarkovsky drift, star events (engulfment keeps the surviving outer part, mass-loss widening, companion truncation).

**Link to exoplanet evidence.** Debris around known planet hosts: 10 cool discs around 37 hosts, 29 +/- 9 percent for Sun-like hosts, consistent with the 20.2 +/- 2.0 percent DUNES rate [S: Herschel exoplanet-host sample, arXiv:1403.6186]. That paper also finds more dust with lower planet masses and lower metallicity, which fits `giant-stirring` and `instability-depletion`: heavy planets sweep their belts, light ones leave them. Another paper finds no significant correlation between radial-velocity planet presence and disc properties [S, title seen only: arXiv:2005.03573], so the link is weak and stays a stand-in: **use the architecture's outcome as the driver, do not impose a correlation**. The correlation then falls out of the model, and the test is whether it appears.

### 3.1 How the generator picks the cause (rule, not a roll)

Run after the planet architecture (GEN.728) has the final chain and giant history:
1. **Inner belt candidates:** every gap in the oligarch chain where the next mass step is below the isolation mass, and the zone between the snow line and the first giant. Each gets `failed-growth`. If a giant lies within the chaotic-zone ratio of the zone, add `giant-stirring`: mass multiplied by a survival factor that falls with giant mass.
2. **Instability flag:** the architecture already moves and scatters giants. If it moved a giant through the zone, apply `instability-depletion` (survival stand-in about 0.01 to 0.3, widest range until read).
3. **Outer disc:** past the last giant, always built from the disc solids out to a disc edge (stand-in), then ground by Wyatt for the star's age. Occurrence of a **detectable** cold disc is then an output, to be compared with Sibthorpe 2018 and DEBRIS (17 percent F-K) as a check, not an input. Where the generator already stores the GEN.787 descriptor, that descriptor becomes the outer disc's summary and the two stop being separate things.
4. **Rare cause:** `collisional-debris` at a small stand-in rate per system age tied to the system's dynamical history.
5. Write `origin`, `formation_snow_line`, and a one-line history for the belt page.

Stand-ins this section adds: survival factors, the disc edge, the stirring factor (all T2 or T3, L0 to L1 until read).

## 4. Counting Sol, so we know what to compare to

[R] unless marked. Orders of magnitude; the lane reads the surveys.

| Population | Bodies of 100 km and up | Bodies near 400 km and up ("round" scale) | Largest |
|---|---|---|---|
| Main belt | about 200 to 250 | 4 or 5 (Ceres, Vesta, Pallas, Hygiea; Interamnia nearly) | Ceres, 940 km |
| Jupiter Trojans | about 25 to 30 | none | Hektor, about 220 km |
| Kuiper belt and scattered disc | order 1e5 [R: Bernstein 2004 fits, Fraser 2014; the Bernstein paper puts the classical mass "predominately in the form of ~100 km bodies", S] | roughly 15 to 20 candidates (Pluto, Eris, Haumea, Makemake, Gonggong, Quaoar, Sedna, Orcus, Salacia, Varda, Huya about 410 km [S: arXiv:2205.12882]) | Eris and Pluto, about 2,300 km |
| Moons of 50 to 500 km radius | about 25 to 30 (Saturn: Mimas, Enceladus, Hyperion, Phoebe, Janus; Uranus: Miranda, Puck; Neptune: Proteus, Nereid, Larissa; others) | only a few at the top | Tethys, 530 km |
| IAU dwarf planets | 5 (Ceres, Pluto, Eris, Haumea, Makemake) | | |

The cold classical belt has **no objects of 400 km or more**; its size distribution has an exponential cutoff there [S: OSSOS, arXiv:2107.06120]. The big bodies are in the excited (scattered) population, which is the point of section 6. 

## 5. Do systems need more Class D? Yes, but as ledger members

### 5.1 What the generator does today [REPO]
* No planet can be lighter than 0.1 Earth masses (0.6 beyond the snow line). Dwarf planets as planets: **zero**.
* `DWARF_PLANET` exists only as a moon-parent kind and a composition draw (`moon_composition.py`, `small_body_composition.draw_icy_body`). Nothing makes a dwarf planet in a belt or an outer disc.
* There is no outer disc at all, so the part of Sol with the most dwarf planets does not exist in any system.
* Belts carry no members, so a belt with Ceres in it is a belt with nothing in it.
* I could not run the census this session: the stellar data files (SC1 isochrones) are not on this machine, so the per-star counts of Class D from a sample run are not measured. The Class D count per system from the code above is not zero only through moons and the rare airless icy planet. The building lane takes the number from `scripts/test.py` census once the data are fetched; section 7 item 1.

### 5.2 Population model (recommended)
Do not add Class D to the planet slot loop. Draw them from the belt and disc size laws:

* `N(>D) = N100 (D / 100 km)^(-q)` per population, with the slope `q` and the cutoff from the size-law note (`small-body-size-laws.md`: main belt cumulative q 3.0 from 5 to 40 km (Ivezic's differential 4) and a steeper stand-in above, see `belt_ledger.ROCKY_LAW`; Kuiper bright slope, hot and cold from Fraser 2014 [REPO table]).
* `N100` from the belt mass and size law, not typed in; the top body is the order statistic `N(>D) = 1` (ledger draft 3.1). Test target: main belt with 0.0005 Earth masses gives a Ceres-size top body and 200 to 250 bodies of 100 km and up; Kuiper-type with 0.02 Earth masses gives order 1e5 bodies of 100 km and up and the 15 to 20 round candidates.
* **Dwarf planet line:** a body is a dwarf planet when its mass is above the hydrostatic mass for its composition (ice about 1e20 to 1e21 kg, rock above that, ledger draft 3.2) and the orbit-clearing number is under the planet line. The count follows from the size law; **do not put a count in as a target**.
* **Class D label:** applied to every L1 body that is airless and at least 1 percent ice (conditions rule). The rocky ones inside the snow line are Class C or A by the rules the planet classes already hold, or have no letter if they are under the dwarf line.
* **What the user sees:** a rich outer disc lists tens of individual dwarf planets; a main-belt analogue lists Ceres-class bodies; a thin belt lists none. Mass is conserved.

Expected effect, order of magnitude for a Sun-like star with a Sol-like outer disc: from zero dwarf planets per system to about 5 to 20, and L1 members (the 1e19 kg and up cut recommended in the ledger draft section 11.1, rule E, cap 300) from none to tens or hundreds per rich disc. For an M dwarf the disc is lighter and the count falls with it; for a system where `instability-depletion` happened the inner belt has a few and the outer more.

### 5.3 Moons
Sol has about 25 to 30 moons in the 50 to 500 km radius range, mostly irregular captured moons and small regular ones. The moon draw already has a `captured` row. The scatter phase (section 6) is the natural source for captured irregular moons (giants capture scattered bodies); if the moon counts are low after the architecture lands (GEN.921 already tracks moons), tie capture to the same event.

## 6. The scatter phase

**Yes, system generation needs one.** Without it the outer disc would be one flat ring and there would be no reason for a dwarf planet to be on an eccentric orbit, for Trojans to exist, or for a belt to be partly emptied. Boss's idea fits the physics.

What it does, once per system, after the architecture, before belts are finalised and before validation:
1. **Inputs:** final giants (mass, a, e), their migration history, the outer disc (mass, edge) and each inner belt.
2. **Scatter fraction:** the fraction of each disc that was scattered is a function of the last giant's mass and migration extent (stand-in; Levison et al. [S: arXiv:0712.0553] needs Neptune scattering and a disc truncated at about 30 to 35 AU for the observed structure). A system whose giants never moved keeps a cold disc and no scattered population; a system with a Nice-type instability makes a hot population.
3. **Outcome classes** (each a tag on a ledger member and a share of the disc mass): `cold-classical` (stays, low eccentricity and inclination, no large bodies, as OSSOS), `hot-classical`, `resonant` (3:2 and 2:1 with the outer giant), `scattered-disc` (high eccentricity, perihelion near the giant), `detached`, `centaur` (crossing the giants), `trojan` (1:1 with a giant, L4 and L5), `ejected` (leaves the system; contributes to rogue and interstellar bodies, which the generator already has).
4. **Implant:** a small share of outer bodies is implanted into the inner belt (the C and D taxon mix in Sol's outer main belt) and a share is captured as irregular moons.
5. **Write** the mass shares, the tags on the L1 draw, and one history line ("when the outer giant moved from 18 to 30 AU it scattered 82 percent of the disc").

The step is a rule over the architecture's own events; it needs no N-body run. A sparse fast check (does a drawn scattered body cross a giant's Hill sphere within the system age) can be a test, not a runtime cost.

**Limits to name:** no direct observation of how many other systems did this. Exoplanet evidence is indirect (debris discs at 7 to 40 AU, Fomalhaut-type narrow rings). We declare the scatter fraction a stand-in (tier educated guess from physics, tweak L1 to L3) and keep it a single registered row.

## 7. TODO items

| Item | What | Owner and train (best guess, can move) |
|---|---|---|
| GEN.1087 | Census of Class D and dwarf-size bodies | Aux 3, now |
| GEN.1060 | Research read (extended with the sources here) | Sources 2, now |
| GEN.1088 | belt_origin.py | Planets 1 |
| GEN.1089 | scatter.py | Core 4 |
| GEN.1090 | dwarf_population.py | Sources 2 |
| GEN.1091 | Class D wording (docs only) | Aux 2, done: the wording is a conditions rule, stated in mass-object-class.md, small-bodies-composition.md, D-small-bodies-moons.md, unified-system-barycentre.md and the comment on `tuning.PLANET_CLASSES["D"]` |
| GEN.1068, GEN.1069 | Belt birth and storage carry the origin (lock) | Planets 1 and Sources 2, train 15 |
| GEN.1092 | Wire scatter and origin (lock) | Sources 2 with Planets 1, train 16 |
| GEN.1082 | Captured moons from the scatter step wait for the census | with GEN.921 |

The belt chain itself is in [belt-evolution.md](belt-evolution.md).

## 8. Stand-ins and tier

| Value | Tier | Tweak | Replace by |
|---|---|---|---|
| Survival factors for `giant-stirring` and `instability-depletion` | educated guess from physics | L1 | Walsh 2011, Raymond and Izidoro 2017 |
| Disc edge and scatter fraction | educated guess from physics | L1 to L3 | Levison 2008 simulations; DEBRIS fit |
| `collisional-debris` rate | educated guess from needed output | L3 | transient-dust statistics |
| Cut for L1 members | design | L0 | Boss decision in ledger draft 11.1 |
| Dwarf-planet mass line | confident unverified | L0 | Bierson and Nimmo |

## 9. Questions for Boss (defaults taken; change on his word)

1. **Class D scope.** Default: the class stays "small, airless, at least 1 percent ice", and applies to belt members and dwarf planets as well as moons. The alternative is a new class letter for dwarf planets of belts, leaving D for moons only.
2. **Where dwarf planets live in lists.** Default: they appear on the system page as members of their belt or outer disc, not in the planet list.
3. **Scatter phase size.** Default: a rule over the architecture's events with stored mass shares, no N-body run.

## Sources
* Mass of the Kuiper Belt, arXiv:1810.09771 [S]
* OSSOS exponential cutoff, arXiv:2107.06120 [S]
* Asteroids Were Born Big, arXiv:0907.2512 [S]
* Towards initial mass functions for asteroids and KBOs, arXiv:1004.0270 [S]
* Bernstein et al., size distribution of trans-Neptunian bodies, arXiv:astro-ph/0308467 [S]
* Kenyon and Bromley chapter, arXiv:0704.0259 [S]
* Levison et al., Kuiper belt structure during the instability, arXiv:0712.0553 [S]
* Sibthorpe et al., DEBRIS Sun-like sample, arXiv:1803.00072 [S]
* Herschel exoplanet-host sample, arXiv:1403.6186 [S]
* Huya occultation, arXiv:2205.12882 [S]
* Repo: `physics/data/planet_classes.json`, `tuning.py` (PLANET_MIN_MASS_*), `composition/small-bodies-composition.md`, `architecture/mass-object-class.md`, ledger draft, belt root-cause report.
