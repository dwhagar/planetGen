# The atmosphere mixture of a body: from regimes to mole fractions (GEN.451, GEN.479, GEN.454, GEN.499)

Status: every family implemented and tested (`physics/mixture_draw.py`, `draw_mixture`) and called by the generator since GEN.411 (the fixed point
of `physics/surface_air.py` replays it at every step). Lane I stores the result and owns the pressure. The trace and minor gases of table `m` are drawn after the life stage by
`physics/mixture_trace.py` (`draw_trace_gases`, which adds N2O, PH3, HCl, HF, NO2 and CH3Cl at stand-in amounts where the rows draw them
present); the ozone rows are lane F's.

## The split

| Step | Owner | Where |
|---|---|---|
| Which regime the air is in (redox, delivery, outgassed class, retention class, envelope) | GEN.478 | `physics/atmosphere_recipes.py`, 104 rows |
| Free energies of the gases and condensates | this | `physics/species_thermo.py`, from the NASA Glenn database (99 species) |
| The equilibrium gas of an element budget at T and P, optionally over a mantle of a given fO2 | this | `physics/gas_equilibrium.py` |
| Regimes to mole fractions of the gases of the table | this | `physics/mixture_draw.py`, `draw_mixture` |
| The subject of a generated body: the quantities the tables and the chemistry read, built from its stored fields and its star | GEN.756 | `physics/mixture_subject.py`, `subject_of_body` |
| What an ocean takes out of the drawn air: the sulfur uptake (the water vapour cap is the air solve's saturation) | GEN.757 | `physics/mixture_sinks.py`, `sulfur_remaining_fraction` |
| Total pressure, freeze-out of what the surface condenses, sum of partials = pressure, storage, accessor | lane I | GEN.411, GEN.422: `physics/surface_air.py`, `physics/atmosphere.py` (built 2026-10-24) |
| Clouds, hazes, ozone | lanes G and F | GEN.414, GEN.413 |

A class letter is never an input to the chemistry. The recipe tables say which regime applies; the composition is what the element
budget, the temperature, the pressure and the oxygen fugacity make of it. The research's per-regime speciations (part 05 section 3.4)
are the test of the solver, not a table the generator reads.

## Interface

    draw = mixture_draw.draw_mixture(Subject("A", {...}))         # any class
    draw.mole_fractions      # gas-table formulas, sums to 1
    draw.drawn               # the recipe rows: {"a": "a6", "b": "b2", "c": "c3", "c-after": "c5"}
    draw.redox_diw           # log10 fO2 above the iron-wuestite buffer [R1951]
    draw.inventory           # atoms per 100 g of body: H, C, N, S, Ar
    draw.constants           # every design constant this draw took, for the record
    mixture_draw.required_inputs("A")   # the quantities the subject must carry

`required_inputs` is computed from the recipe tables, so the generator's wiring (GEN.479) cannot drift from what the tables read. A
missing quantity raises `MixtureError` naming it. The draws go through `util/draw.py` after the recipe draw, in a fixed order: a change
of the order is a generator-epoch change.

## The subject of a body and the ocean's sinks (GEN.756, GEN.757)

Two small pieces of GEN.411 that need nothing of its climate loop; the loop passes in what it owns.

`mixture_subject.subject_of_body(body, ts_k=, pressure_bar=, ocean_present=, long_lived_geology=, water_lost=, star_distance_au=)`
returns the `Subject` of a planet or moon, holding exactly `required_inputs(class)`. It derives the mass, radius and escape speed, the
age of the star (the reference age of a starless body, as `interior_heat` takes it), the insolation, the equilibrium temperature at the
body's stored Bond albedo, the XUV exposure index (`activity.exposure_index`), the position against the star's snow line, the
runaway-greenhouse ratio against the Kopparapu inner edge, the zone against the star's habitable zone, the host type, and for a giant
the temperature of its 1 bar layer (`giant_draw.ONE_BAR_OVER_EFFECTIVE` times its effective temperature), its dayside and its
effective temperature. The surface temperature and pressure, the ocean and the geology are the loop's. A moon passes the distance of
its parent from the star. A quantity that cannot be built is refused by name. Limits: the dayside temperature is the no-redistribution
bound (`DAYSIDE_OVER_EQUILIBRIUM`, `(8/3)^(1/4)`: the generator has no redistribution model), the layer pressure of a giant is 0.1 bar
(`PHOTOSPHERE_BAR`), a hot A, B or O star takes the F rows and a brown dwarf the M rows, and the snow ratio uses the body's present
distance because no formation distance is stored.

`mixture_sinks.sulfur_remaining_fraction(formula, ts_k, gravity, ocean_column_kg_m2)` keeps `1 / (1 + kH m g M / P0)` of each sulfur gas (SO2, H2S: Henry's law over `m` kg of ocean per square metre). The water vapour cap over an ocean is the saturation of the air solve, and the outgassed inventory and the ocean and carbonate sinks of the carbon, water and nitrogen are GEN.726's. A function that capped the whole drawn mixture in one call (`apply_ocean_sinks`) was removed by TEST.524: nothing called it. Limits: the physical solubility only (both gases dissociate in water, so a real ocean takes more) and the Henry law's 273.15 to 323.15 K with the nearest end standing in outside.

## Where the degassed sulfur goes (GEN.815)

The supply law brings sulfur to the surface; the sulfur sink decides how much of it the air keeps (`physics/sulfur_balance.py`). It
replaced the stand-in `SULFUR_SURFACE_KEPT` (1e-5 to 1e-2 of the sulfur stayed in the air with an ocean; removed when `surface_air` called this
module), under which a wet world could hold hundreds of bar of S2. Two solid sinks take sulfur out of the gas: sulfate (SO2 and SO3 with calcium as CaSO4) and sulfide (H2S, S2
and COS with iron as FeS and FeS2). Which one a mole goes to is the redox of the gas it arrived in: `oxidised_share_of` reads the share of
the sulfur that is bonded to oxygen from the equilibrium gas of the outgassing stage.

* **Capacity.** The crust the volcanism built holds the calcium and the iron: `crust_column_kg_m2(released_j, radius_m)` is the Earth
  twin's crust per processed mantle (`CRUST_PER_PROCESSED_MANTLE`) times the mantle `surface_pressure` says the interior has processed,
  over the area. One calcium binds one sulfur as sulfate, one iron one sulfur as sulfide. Sulfur beyond the capacity stays in the air whole,
  so a body whose inventory is far larger than its crust (the open question of the inventory brackets, GEN.481) shows it.
* **Rate.** The sulfur that goes to the sinks arrived at its mean rate; the air holds `min(1, tau / age)` of it (the steady state of
  outgassing against removal, as the carbonate balance is). `tau` is the lifetime against rain-out and oxidation over liquid water
  (`WET_LIFETIME_DAYS`, 0.5 to 7 days) and against sulfation of the surface rock on a dry surface (`DRY_LIFETIME_YEARS`, 0.5 to 5 million
  years, [R1343] for Venus).
* **Output.** `balance(...)` returns the sulfur in the air, in the sulfate sink and in the sulfide sink, mol per square metre; they
  sum to the degassed sulfur. The air's sulfur is what the caller puts into the closed system for the equilibrium, and the ocean's
  Henry's-law uptake (`mixture_sinks`) acts on it afterwards.

Limits (stated): sulfate and carbonate share the calcium and the balance does not take the carbonate's share (so the sulfate capacity is
an upper bound where a carbonate sink is full); the mean-rate steady state gives a dry world almost no sulfur in the air, while Venus holds
about 14 mbar of SO2, which only recent resupply explains (Fegley and Prinn say so too); the lifetimes are brackets, not rate laws of
temperature; the share of the sulfur the melt keeps is `surface_pressure.DEGAS_EFFICIENCY`. Five stand-ins, registry rows AM45 to AM49.

## Where the carbon of a world with no land goes (GEN.814)

Silicate weathering (`surface_pressure.carbonate_pco2_bar`) needs land, and an ocean world has none, so its carbon dioxide had no sink: the
findings of `reports/air-fixed-point-findings-2026-10-24.md` (finding 3) have 31 of 40 water-rich worlds hot (median 1470 K). The model is
`physics/waterworld_carbon.py`, af[R2381, R2489] (Nat. Commun. 8, 15423).
`sink_kind(sealed, land_fraction)` is the dispatch the air solve takes: sealed seafloor, then ocean chemistry; else land, then the
continental balance; else seafloor weathering.

* **The seal.** An ocean heavy enough that the seafloor pressure `g m` reaches `SEAL_PRESSURE_PA` (1 to 3 GPa) stops the mantle melting enough
  to exchange carbon with the water, and high-pressure ice (VI and VII) sits on the floor. The system of air, ocean and shallow crust is then
  closed for carbon `SEAL_TIME_GY` (0.1 Gy) after the last giant impact; the caller gives it the carbon degassed up to then.
* **A sealed world's air is ocean chemistry** (`ocean_pco2`). The carbon is shared between the air, dissolved dioxide, bicarbonate and
  carbonate, and calcite, by the carbonate equilibria (`CARBONATE_FITS`, fresh water) and the charge balance of the cations the early ocean
  leached from `ALTERED_CRUST_KM` (0.1 to 50 km) of basalt (`leached_cations_mol_m2`: sodium 2.8 wt% and calcium 8.5 wt% of the rock). More
  cations than carbon holds the carbon as carbonate and leaves the air nearly empty (the floor `MIN_PCO2_BAR`); more carbon than cations
  leaves the rest to the air at tens of bar. The air solve (`surface_air`) takes this pCO2 at the temperature of its climate loop; a stand-alone solver of the pCO2 and the temperature together (`sealed_state`, the coolest crossing of a temperature scan) was removed by TEST.524: nothing called it.
* **An unsealed world with no land** has working volcanoes and a seafloor that dissolves and precipitates carbonate (`seafloor_pco2_bar`):
  the rate is Arrhenius in the pore-water temperature (`SEAFLOOR_ACTIVATION_KJ_MOL`, 40 to 110), a power of the pCO2
  (`SEAFLOOR_PH_EXPONENT / 1.34`) and of the heat flow, and it is 0.45 of 6.6 Tmol C per year (`SEAFLOOR_SHARE_OF_OUTGASSING`) at the modern
  state. The temperature does the work: at the modern heat flow a 300 K world balances at about half a bar, a 288 K one needs hundreds of
  bar, a 310 K one a few millibar.

What the sweep of Kite and Ford's prior shows (tested): at 300 K, carbon 1e-6 to 1e-1 of the ocean's mass and a rock layer of 0.1 to 50 km, a
sixth of the worlds sit at 0.2 to 20 bar (their "frequent" range), nearly three quarters hold their carbon whole in the ocean (under a
micro-bar), and 7 percent are above 20 bar. The water-world temperatures of the findings report are not recomputed here: that needs Lane I's solve (4c), and
the greenhouse module's water-vapour feedback is too strong (finding 1), so the cold end of the stand-alone solve (removed by TEST.524) sat at the freezing point
(274 K for 240 W/m2) and not under it. With the real greenhouse a carbon-rich ocean (CO2 a tenth of the water's mass) had no balance below water's
critical point, as for Kite and Ford, who find the high-carbon end hot.

Limits (stated): the ocean is an ideal solution with fresh-water constants (seawater's first dissociation constant is about three times the
fresh water's, so a salty ocean holds more bicarbonate at a given pCO2) held to 273 to 373 K; only calcite precipitates (no magnesium or iron
carbonate, no dolomite, no hydrous carbonates, no clathrate); the cations are one leach of a basalt of a drawn thickness, not a reaction
network, and sulfate competes for the same calcium (`sulfur_balance`); the surface ocean's temperature is used throughout; the seal is a
pressure threshold, not a solidus calculation; the seafloor law is one pore-space rate law and not the inverse model of Krissansen-Totton and
Catling (no Michaelis-Menten variant, no pelagic fraction); the class O inventory brackets (GEN.481) reach a CO2 to H2O mass ratio of 0.2,
above the 1e-1 top of Kite and Ford's carbon. Eleven stand-ins, registry rows AM52 to AM62.

### The sinks in the air solve (GEN.830, GEN.823, GEN.837)

`surface_air._shares` takes the carbon that the surface holds at each step's temperature from the sink that `sink_kind` picks for the body
(`water_world_of`: sealed when the ocean's column `g m / area` reaches the drawn seal pressure): the continental balance as before for a
body with land, `seafloor_pco2_bar` for an unsealed ocean world, and `ocean_pco2` for a sealed one, whose supply is also the carbon
degassed by `SEAL_TIME_GY` (`interior_heat.released_by`) and not by the age of the star. Before this, an ocean world had no sink for the
carbon at all, and a class O planet in zone e came out of the solve at 9e9 Pa.

*The hydrogen cap (GEN.837).* The delivered water of recipes b2 and b3 is 0.1 to 0.2 percent of a planet's mass, and the body's own water
(`bulk_water_fraction`, 0.01 to 0.03 percent for a class M) is 3 to 50 times less: the draw gave a dry planet an ocean's worth of steam
(55 to 870 bar at 620 to 843 K). The supply of hydrogen is now limited to `tuning.MANTLE_OUTGAS_SHARE` (0.3, the DESIGN stand-in of GEN.639,
the share of a mantle's water that volcanism brings up) of the water the mantle holds (`Setting.mantle_water_kg`: the bulk water less the
surface's); the hydrogen of CH4, NH3 and H2S is limited with it, since they are made of the same atoms. The recipe rows themselves (b2, b3) are
Aux 17's re-fit against the bodies' bulk water.

Measured on 12 planets of each class at zone e (G2V host, the same seeds before and after):

| class | hot (above 500 K) before | after | median temperature before | after |
|---|---|---|---|---|
| M | 3 | 1 | 286 K | 271 K |
| N | 2 | 0 | 237 K | 224 K |
| O | 11 | 3 | 1557 K | 283 K |
| P | 3 | 0 | 267 K | 234 K |

What is left: four hot worlds of 48 (three class O at 4,300 to 33,000 bar, one class M at 578 bar). They are sealed or carbon-rich worlds whose greenhouse
takes the surface past the critical point of water (647 K, as Venus, whose 92 bar of carbon dioxide holds 735 K on 157 W/m2): the ocean can
no longer stand as a liquid, and the sinks that need one switch off. All the surface water is then the air, a steam envelope that is the inventory of the planet and
not a runaway of the solve: the weight of the air, `g m / area`, is the mass of the air, and `test_planets` bounds it by the planet's own water plus
the volatile share of a CI chondrite (9 percent of the mass) and not by a flat pressure. The greenhouse module's water feedback is too strong
(GEN.812, Aux 4) and a sealed world's CO2 to H2O bracket is above Kite and Ford's (GEN.824); when either moves, these four move toward cool
ocean worlds. The class of such a world is a matter of GEN.646 (class from conditions), which sees a steam envelope.

## Rocky outgassing

1. **Redox.** `a1` (magma-ocean stage: [R4154]) for a body without long-lived geology; for a tectonic body it moves to `a2` (the
   later upper mantle) as `1 - exp(-age / tau)` with `tau` seeded in 0.1 to 1.5 Gyr. A white-dwarf host takes `a4`, a small body `a3`.
   The recipe's class row `a5`/`a6`/`a7` describes the magma-ocean stage only; the chemistry reads the blended number.
2. **Inventory.** Delivered water and CO2 from the `b` class ([R3378: Table 1]: CO2/H2O = 0.2 for the wet classes, dry class
   CO2 0.018 to 0.1 wt%, GEN.813), N and S per carbon in design brackets (BSE and CI chondrite values inside them), 40Ar from the potassium of the rock
   and the age.
3. **The volcano's gas**: equilibrium with the mantle at the magmatic temperature (1100 to 1500 K) and the pressure of the layer, with
   the oxygen fugacity held.
4. **Cooling** in a closed system to the surface temperature, or to the quench temperature (400 to 900 K, [R2686, R4863]) when the surface is colder. Water, sulfur and the other condensables a cool surface takes out are the pressure
   bookkeeping's (GEN.422): this is the gas before that.

Holding fO2 at the surface temperature instead of at the vent was tried first and is wrong: below about 900 K the buffer pins fO2 at
1e-30 bar or less, a gas no longer exchanges oxygen with the mantle, and the result was carbon disulfide worlds.

## The families

`draw_mixture` draws the recipe tables a to l once (`RECIPES.draw_all`) and picks the family from the rows that came out, not from the
class letter (`mixture_draw` module text has the table). Each family is one module over the same solver:

| family | module | chemistry |
|---|---|---|
| `rocky` | `mixture_rocky` | the volcano's gas over the mantle, cooled in a closed system (above) |
| `giant`, `brown-dwarf`, `hot-jupiter`, `envelope` | `mixture_giant` | [R0168] budget at the drawn metallicity and C/O, Visscher condensation (clouds, NH4SH, ice caps) removing what the layer's temperature takes out, adiabat to the layer, Zahnle and Marley quench |
| `venus`, `mars` (`j1`, `j2`) | `mixture_evolved` | the outgassed gas minus its water (`WATER_RESIDUAL`) and most of its sulfur, cooled to the surface temperature (not quenched) over a crust that holds the body's redox: CO2 with a few percent N2 (Mars-like 0.978 CO2, 0.016 N2; the TEST values are 0.96 and 0.019) |
| `abiotic-oxygen` (`j5`) | `mixture_evolved` | the outgassed gas rescaled to a carbon of 0.01 to 0.2 per molecule of an ocean of one to ten Earth oceans, the hydrogen of 30 to 100 percent of it lost and the crust taking up to 90 percent of the oxygen left: nothing holds the oxygen, so the gas has O2 |
| `titan` (`j3`) | `mixture_evolved` | N2 balance (the nitrogen of the rock plus the ammonia ice's: GEN.821), CH4 at 30 to 100 percent of its saturation pressure over the pressure of the air, H2 a trace |
| `frost` (`j4`) | `mixture_evolved` | Raoult over a surface of N2, CH4 and CO ice: `p_i = chi_i p_sat,i(T)`, the nitrogen setting the total |
| `coma` (`k`) | `mixture_coma` | the ice's composition: CO2/H2O, CO/H2O, CH3OH, NH3, HCN and the 67P species as drawn ratios to water, the class row ordering CO2 and CO |
| `lava` (`i` at `tp_k >= 1400 K`, `l1`) | `mixture_lava` | the vapour over an ideal-solution BSE melt: element potentials from the melt's oxides and the held oxygen fugacity, every gas of those elements at its own partial pressure, no iteration |
| `airless` (`l2`, `l1` below 1400 K) | | no gas: `mole_fractions` is empty |

After the family's chemistry, `mixture_escape.apply` (GEN.693) removes the gases a small warm body cannot keep, for the families that
outgas (`rocky`, `venus`, `mars`, `abiotic-oxygen`). The Jeans parameter `lambda = m v_esc^2 / (2 k T_exo)` of each gas against `lambda =
54` (the escape speed at six times the rms speed, kept for 4.5 Gyr, moved by `ln(age/4.5)`) gives the share that stays, a logistic of
width 4; a gas of which less than 0.1 percent stays is gone, the rest is renormalised, and an air of which less than 1 percent stays
is airless. The exobase temperature is a stand-in (`EXOBASE_RATIO`: 4 times the equilibrium temperature, 1.3 for an air that is half
carbon dioxide, times the XUV flux to the 0.25, clamped to 0.5 to 3) fitted to Earth, Mars and Venus and replaced when a thermosphere
model exists. Earth gives hydrogen 15 against 54 (gone) and helium at 0.2 percent kept (a few ppm, as measured); a 3.5 km/s world at
310 K keeps none of its outgassed hydrogen, carbon monoxide or carbonyl sulfide and is airless. The giants, the envelopes, the lava
vapour, the ice-supported air and the Titan-like air (whose hydrogen is a photochemical steady state) are not filtered. The total
pressure is not changed here.

The oxygen cap (`mixture_oxygen`, Boss 10-24): no air is pure oxygen. A lifeless world that lost its water (the abiotic-oxygen family) is
capped at 95 percent oxygen with at least 0.5 percent nitrogen plus argon; a living world is capped at 35 percent (`BIOTIC_O2_MAX`, from the
Phanerozoic peak and the fire limit), for the life stage that sets the oxygen of a living world to apply with `cap_oxygen`. The analysis is
in `/mnt/project-files/reports/pure-oxygen-atmosphere-2026-10-24.md`.

The lava vapour reproduces [R2380]'s vapour pressure over a BSE melt (1.8e-4 bar at 2000 K, 1.1e-2 bar at 2500 K against 2e-4
and 8e-3) and the sequence of [R3059]: Na and K first, then SiO, Mg and Fe from 2000 to 2500 K, SiO above O2 beyond about 2900 K,
mean molar mass 33 g/mol at 2000 K and 42 g/mol at 3000 K. It does so only with the alkali oxides' activity coefficient of 3e-9 to 3e-8
(an ideal solution puts 1 to 6 bar of sodium and potassium over the melt at 2000 K), which is a stand-in (below).

## Stand-ins and limits

- The recipe tables had no rows for the planned rocky classes U, Y and Z in tables a to d, l, m and o, although the rows of `j` are gated
  on them (and `j3`, `j4` need a retention row `d4` or `d5`): `ROCKY` in `scripts/atmosphere_recipes_data.py` now lists them. No body has
  these classes yet, so no output changes.
- GEN.674 closed these gaps of the draw (all unwired until GEN.411, so no generator epoch): (1) a comet whose perihelion is beyond the
  coma's gate (rows `k` pass only inside 3 AU, or for an icy asteroid) drew the `none` of group k and then fell through to the rocky
  family; it is now `airless` (an asteroid still goes on to the rocky rows). (2) In a group whose fallback is `none` (j, k) the weights
  are probabilities of occurrence, not shares: a lone row of weight 0.1 is drawn one time in ten, and a frost row 100 K outside its
  window is no longer normalised to 1 (frost was drawn at 181 K); `c_j1` and `c_j5` now mean what their brackets say (the share of
  eligible worlds that are Venus-like and that keep abiotic oxygen). (3) A body with an ocean and a long-lived geology kept only a
  share of the outgassed carbon and sulfur (`INVENTORY["weathering_sink"]`, 1e-5 to 1e-2): GEN.411 removed that bracket, and the sinks are
  now the carbonate-silicate balance (GEN.726, `surface_pressure`), the sulfate and sulfide sinks of the crust (`sulfur_balance`, GEN.815, wired
  in `surface_air._shares`) and the ocean's uptake of sulfur (`mixture_sinks`, GEN.757), which `surface_air` passes to the draw as `supply`, a callable from the inventory to the share of each
  element degassed (see `surface-pressure-supply-law.md` section 6).
- Gaps that stay, each with its reason (GEN.674 items 3, 5 and 6; the filed finding numbers are in the item): rows `j3` and `j4` are
  gated on the retention rows `d4` and `d5`, and the shoreline (a fit to hot rocky exoplanets, extrapolated to an escape speed of 1.2
  km/s) gives a Pluto-like body `d1` (airless), so the Pluto-like family is never drawn for it. Boss decided (fork 7, 10-18) that
  these families are an output of the escape and retention physics; the retention physics of cold volatile bodies (the sublimation
  supply against hydrodynamic escape) is not built, and changing the gate without it would be a design flag, which he ruled out. GEN.751
  added K2, KO, K2O and Na2O to the gas table and the species table (research part 01 section 5, computed from the NASA fits, no
  phase data; the file spells aluminium `AL`), so the lava vapour lists them (KO is about 4e-7 bar over a BSE melt at 2500 K). GEN.818 then gave the melt an Al2O3 activity (its mole fraction, `mixture_lava.ALUMINA_ACTIVITY_COEFFICIENT` = 1, tier 3: [R1341] find the activity-coefficient models of the BSE melt agree with an ideal solution except for the alkalis,
  abstract only), added AlO, Al2O, AlO2 and Al2O2 and the liquid Al2O3 to the species table (99 species), and put Al in the vapour's
  elements: a BSE melt at 2500 K now has about 2e-7 of its vapour as aluminium, mostly AlO and Al (Al is refractory), 6e-5 at 3000 K.
- Design brackets this module added, each to be replaced by the research item named: `WATER_RESIDUAL`, `SULFUR_SINK`,
  `WATER_PHOTOLYSED`, `OCEAN_WATER_WT_PERCENT`, `CARBON_PER_OCEAN_WATER`, `OXYGEN_SINK`, `TITAN_HUMIDITY`, `TITAN_HYDROGEN`, `AMMONIA_OUTGASSED`, `AMMONIA_PER_WATER`, `AMMONIA_TO_N2`, `FROST_ICE`
  (`mixture_evolved`; GEN.481 and the ice composition item), `ALKALI_ACTIVITY_COEFFICIENT` and `COMPOSITION_SCATTER` (`mixture_lava`; the
  melt activity model of GEN.499 and the bulk composition of lane J), the k1 to k3 ordering rule of the coma.

- N/C, S/C, K and the 40Ar degassing fraction are design brackets (`INVENTORY`), replaced by GEN.481 (volcanic delivery) and GEN.492
  (radioisotopes) when they land. The delivered water classes are research values [R3378].
- Sulfur is S2 in the gas; the solver has no S8 or S3, and sulfur condenses below about 700 K anyway (GEN.422).
- Above about 1600 K at low pressure the gas holds atoms and radicals (H, O, OH, C, N, S) that the gas table does not list: the draw
  refuses with `MixtureError` when they exceed `UNLISTED_LIMIT`. Those bodies are lava worlds (GEN.499); the gas table additions are
  filed with the TODO thread.
- Missing from the rocky path, to come: the oxygen that life makes (life is computed after the planet; the ordering belongs to GEN.422)
  and the temperature feedback of the weathering sink (GEN.481).
