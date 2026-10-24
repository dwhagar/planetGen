# The atmosphere mixture of a body: from regimes to mole fractions (GEN.451, GEN.479, GEN.454, GEN.499)

Status: rocky outgassing implemented and tested (`physics/mixture_draw.py`), not yet called by the generator. The giant, brown dwarf, hot
Jupiter, envelope and lava-world parts follow (GEN.454, GEN.499). Lane I (GEN.411) stores the result and owns the pressure.

## The split

| Step | Owner | Where |
|---|---|---|
| Which regime the air is in (redox, delivery, outgassed class, retention class, envelope) | GEN.478 | `physics/atmosphere_recipes.py`, 104 rows |
| Free energies of the gases and condensates | this | `physics/species_thermo.py`, from the NASA Glenn database (87 species) |
| The equilibrium gas of an element budget at T and P, optionally over a mantle of a given fO2 | this | `physics/gas_equilibrium.py` |
| Regimes to mole fractions of the gases of the table | this | `physics/mixture_draw.py`, `draw_mixture` |
| Total pressure, freeze-out of what the surface condenses, sum of partials = pressure, storage, accessor | lane I | GEN.411, GEN.422 |
| Clouds, hazes, ozone | lanes G and F | GEN.414, GEN.413 |

A class letter is never an input to the chemistry. The recipe tables say which regime applies; the composition is what the element
budget, the temperature, the pressure and the oxygen fugacity make of it. The research's per-regime speciations (part 05 section 3.4)
are the test of the solver, not a table the generator reads.

## Interface

    draw = mixture_draw.outgassed_mixture(Subject("A", {...}))     # rocky; the other families join with their PRs
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

## Stand-ins and limits

- N/C, S/C, K and the 40Ar degassing fraction are design brackets (`INVENTORY`), replaced by GEN.481 (volcanic delivery) and GEN.492
  (radioisotopes) when they land. The delivered water classes are research values (Ortenzi 2020).
- Sulfur is S2 in the gas; the solver has no S8 or S3, and sulfur condenses below about 700 K anyway (GEN.422).
- Above about 1600 K at low pressure the gas holds atoms and radicals (H, O, OH, C, N, S) that the gas table does not list: the draw
  refuses with `MixtureError` when they exceed `UNLISTED_LIMIT`. Those bodies are lava worlds (GEN.499); the gas table additions are
  filed with the TODO thread.
- Missing from the rocky path, to come: the carbonate-silicate CO2 sink of a body with an ocean and plate tectonics, and the oxygen
  that life makes (life is computed after the planet; the ordering belongs to GEN.422).
