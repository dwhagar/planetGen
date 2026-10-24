# The atmosphere mixture of a body: from regimes to mole fractions (GEN.451, GEN.479, GEN.454, GEN.499)

Status: every family implemented and tested (`physics/mixture_draw.py`, `draw_mixture`), not yet called by the generator. Lane I
(GEN.411) stores the result and owns the pressure. The trace and minor gases of table `m` are drawn after the life stage by
`physics/mixture_trace.py` (`draw_trace_gases`, which adds N2O, PH3, HCl, HF, NO2 and CH3Cl at stand-in amounts where the rows draw them
present); the ozone rows are lane F's.

## The split

| Step | Owner | Where |
|---|---|---|
| Which regime the air is in (redox, delivery, outgassed class, retention class, envelope) | GEN.478 | `physics/atmosphere_recipes.py`, 104 rows |
| Free energies of the gases and condensates | this | `physics/species_thermo.py`, from the NASA Glenn database (89 species) |
| The equilibrium gas of an element budget at T and P, optionally over a mantle of a given fO2 | this | `physics/gas_equilibrium.py` |
| Regimes to mole fractions of the gases of the table | this | `physics/mixture_draw.py`, `draw_mixture` |
| Total pressure, freeze-out of what the surface condenses, sum of partials = pressure, storage, accessor | lane I | GEN.411, GEN.422 |
| Clouds, hazes, ozone | lanes G and F | GEN.414, GEN.413 |

A class letter is never an input to the chemistry. The recipe tables say which regime applies; the composition is what the element
budget, the temperature, the pressure and the oxygen fugacity make of it. The research's per-regime speciations (part 05 section 3.4)
are the test of the solver, not a table the generator reads.

## Interface

    draw = mixture_draw.draw_mixture(Subject("A", {...}))         # any class
    draw.mole_fractions      # gas-table formulas, sums to 1
    draw.drawn               # the recipe rows: {"a": "a6", "b": "b2", "c": "c3", "c-after": "c5"}
    draw.redox_diw           # log10 fO2 above the iron-wuestite buffer (Hirschmann 2021)
    draw.inventory           # atoms per 100 g of body: H, C, N, S, Ar
    draw.constants           # every design constant this draw took, for the record
    mixture_draw.required_inputs("A")   # the quantities the subject must carry

`required_inputs` is computed from the recipe tables, so the generator's wiring (GEN.479) cannot drift from what the tables read. A
missing quantity raises `MixtureError` naming it. The draws go through `util/draw.py` after the recipe draw, in a fixed order: a change
of the order is a generator-epoch change.

## Rocky outgassing

1. **Redox.** `a1` (magma-ocean stage: Sossi 2020) for a body without long-lived geology; for a tectonic body it moves to `a2` (the
   later upper mantle) as `1 - exp(-age / tau)` with `tau` seeded in 0.1 to 1.5 Gyr. A white-dwarf host takes `a4`, a small body `a3`.
   The recipe's class row `a5`/`a6`/`a7` describes the magma-ocean stage only; the chemistry reads the blended number.
2. **Inventory.** Delivered water and CO2 from the `b` class (Ortenzi 2020 Table 1: CO2/H2O = 0.2 for the wet classes, dry class
   CO2 0.6 wt%), N and S per carbon in design brackets (BSE and CI chondrite values inside them), 40Ar from the potassium of the rock
   and the age.
3. **The volcano's gas**: equilibrium with the mantle at the magmatic temperature (1100 to 1500 K) and the pressure of the layer, with
   the oxygen fugacity held.
4. **Cooling** in a closed system to the surface temperature, or to the quench temperature (400 to 900 K, Zahnle and Marley 2014,
   Liggins 2022) when the surface is colder. Water, sulfur and the other condensables a cool surface takes out are the pressure
   bookkeeping's (GEN.422): this is the gas before that.

Holding fO2 at the surface temperature instead of at the vent was tried first and is wrong: below about 900 K the buffer pins fO2 at
1e-30 bar or less, a gas no longer exchanges oxygen with the mantle, and the result was carbon disulfide worlds.

## The families

`draw_mixture` draws the recipe tables a to l once (`RECIPES.draw_all`) and picks the family from the rows that came out, not from the
class letter (`mixture_draw` module text has the table). Each family is one module over the same solver:

| family | module | chemistry |
|---|---|---|
| `rocky` | `mixture_rocky` | the volcano's gas over the mantle, cooled in a closed system (above) |
| `giant`, `brown-dwarf`, `hot-jupiter`, `envelope` | `mixture_giant` | Asplund 2021 budget at the drawn metallicity and C/O, Visscher condensation (clouds, NH4SH, ice caps) removing what the layer's temperature takes out, adiabat to the layer, Zahnle and Marley quench |
| `venus`, `mars` (`j1`, `j2`) | `mixture_evolved` | the outgassed gas minus its water (`WATER_RESIDUAL`) and most of its sulfur, cooled to the surface temperature (not quenched) over a crust that holds the body's redox: CO2 with a few percent N2 (Mars-like 0.978 CO2, 0.016 N2; the TEST values are 0.96 and 0.019) |
| `abiotic-oxygen` (`j5`) | `mixture_evolved` | the outgassed gas rescaled to a carbon of 0.01 to 0.2 per molecule of an ocean of one to ten Earth oceans, the hydrogen of 30 to 100 percent of it lost and the crust taking up to 90 percent of the oxygen left: nothing holds the oxygen, so the gas has O2 |
| `titan` (`j3`) | `mixture_evolved` | N2 balance, CH4 at 30 to 100 percent of its saturation pressure over the pressure of the air, H2 a trace |
| `frost` (`j4`) | `mixture_evolved` | Raoult over a surface of N2, CH4 and CO ice: `p_i = chi_i p_sat,i(T)`, the nitrogen setting the total |
| `coma` (`k`) | `mixture_coma` | the ice's composition: CO2/H2O, CO/H2O, CH3OH, NH3, HCN and the 67P species as drawn ratios to water, the class row ordering CO2 and CO |
| `lava` (`i` at `tp_k >= 1400 K`, `l1`) | `mixture_lava` | the vapour over an ideal-solution BSE melt: element potentials from the melt's oxides and the held oxygen fugacity, every gas of those elements at its own partial pressure, no iteration |
| `airless` (`l2`, `l1` below 1400 K) | | no gas: `mole_fractions` is empty |

The lava vapour reproduces Kite et al. 2016's vapour pressure over a BSE melt (1.8e-4 bar at 2000 K, 1.1e-2 bar at 2500 K against 2e-4
and 8e-3) and the sequence of Miguel et al. 2011: Na and K first, then SiO, Mg and Fe from 2000 to 2500 K, SiO above O2 beyond about 2900 K,
mean molar mass 33 g/mol at 2000 K and 42 g/mol at 3000 K. It does so only with the alkali oxides' activity coefficient of 3e-9 to 3e-8
(an ideal solution puts 1 to 6 bar of sodium and potassium over the melt at 2000 K), which is a stand-in (below).

## Stand-ins and limits

- The recipe tables had no rows for the planned rocky classes U, Y and Z in tables a to d, l, m and o, although the rows of `j` are gated on
  them (and `j3`, `j4` need a retention row `d4` or `d5`): `ROCKY` in `scripts/atmosphere_recipes_data.py` now lists them. No body has these
  classes yet, so no output changes. Even so a Pluto-like body (escape speed 1.2 km/s) gets the retention row `d1` (airless) from the
  shoreline, so `j4` is never drawn for it; filed with TODO.
- The recipe rows `k` carry no gate on the perihelion: a comet or icy asteroid is drawn a coma whatever its distance; filed with TODO.
- Design brackets this module added, each to be replaced by the research item named: `WATER_RESIDUAL`, `SULFUR_SINK`,
  `WATER_PHOTOLYSED`, `OCEAN_WATER_WT_PERCENT`, `CARBON_PER_OCEAN_WATER`, `OXYGEN_SINK`, `TITAN_HUMIDITY`, `TITAN_HYDROGEN`, `FROST_ICE`
  (`mixture_evolved`; GEN.481 and the ice composition item), `ALKALI_ACTIVITY_COEFFICIENT` and `COMPOSITION_SCATTER` (`mixture_lava`; the
  melt activity model of GEN.499 and the bulk composition of lane J), the k1 to k3 ordering rule of the coma.

- N/C, S/C, K and the 40Ar degassing fraction are design brackets (`INVENTORY`), replaced by GEN.481 (volcanic delivery) and GEN.492
  (radioisotopes) when they land. The delivered water classes are research values (Ortenzi 2020).
- Sulfur is S2 in the gas; the solver has no S8 or S3, and sulfur condenses below about 700 K anyway (GEN.422).
- Above about 1600 K at low pressure the gas holds atoms and radicals (H, O, OH, C, N, S) that the gas table does not list: the draw
  refuses with `MixtureError` when they exceed `UNLISTED_LIMIT`. Those bodies are lava worlds (GEN.499); the gas table additions are
  filed with the TODO thread.
- Missing from the rocky path, to come: the carbonate-silicate CO2 sink of a body with an ocean and plate tectonics, and the oxygen
  that life makes (life is computed after the planet; the ordering belongs to GEN.422).
