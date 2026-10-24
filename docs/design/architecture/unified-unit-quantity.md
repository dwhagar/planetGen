# One unit and quantity class for everything the project measures

**Status:** plan (decided). Written after 11:51 PM EDT on 2026-10-17 from `main` at bc97f75; approved by Boss on the six cards Q1 to Q6 of section 9, taps 11:56 to 11:57 PM EDT on 2026-10-17 (as the coordinator relayed), each the recommended option. Nothing is built. It adds a section to the unified mass-object plan ([unified-star-like-objects.md](unified-star-like-objects.md), [mass-object-class.md](mass-object-class.md)). The provisional labels UN1 to UN12 of section 8 are filed as: UN1 = GEN.568, UN2 = GEN.569, UN3 = GEN.570, UN4 = DB.76, UN5 = GEN.571, UN6 = UX.139, UN7 = UX.140, UN8 = GEN.573, UN9 = GEN.574, UN10 = DB.77, UN11 = DOC.70, UN12 = GEN.575; the exact-factor step of UN5 (Q6) is its own item, GEN.572; the relation items of section 11 are filed as UN13 = GEN.585, UN14 = GEN.586, UN15 = GEN.587, UN16 = GEN.588 (cards R1 to R4, 12:41 to 12:43 AM EDT on 2026-10-18). Two changes from section 8's text: the removal of the duplicate constant names is the last step of GEN.571 (after its sweeps merge), not of GEN.568, so UN1 has no wait on UN5; and Q6 is its own item because the exact factors change a scatter output (see GEN.572). Boss tapped 'Pull into 8.1' on the card for GEN.572 at 12:04 AM EDT on 2026-10-18 (as the coordinator relayed): it is built in 8.1, before the one wiped rescatter, and is a prerequisite of GEN.301.

**Boss, 11:51 PM EDT:** "I want part of this plan to be the creation of the unified 'unit' code, it will store any kind of unit we use, it stores and updates itself if you have a Watt unit then when you put in Watt it convers it to every other format we use immediately so that when the time comes it can just be referenced. This is what we did with point in space, which every mass object must have a point in space."

Evidence tags: **[M]** measured in the repository (count, grep, AST walk), **[E]** my estimate from those counts, **[C]** computed here by arithmetic. Nothing here is a citation of a paper, so there are no [R] items.

## Decided (Boss, 11:56 to 11:57 PM EDT on 2026-10-17, card taps, all on the recommended option)

- **Q1** Other units are worked out lazily on first read and always in step (the point-in-space way).
- **Q2** A field holds a plain float in its declared native unit; a `Quantity` is built at the edges.
- **Q3** A column's stored unit is declared per column in the registry, chosen for precision.
- **Q4** The unit suffix in names stays, and a test checks it equals the declared unit.
- **Q5** Rolled out inside the host chain items, file group by file group, each sweep bit-identical on the golden set; the registry and the class first.
- **Q6** The two rounded factors (63241.077 AU per ly, 3.26156 ly per pc) become the exact ones in the first sweep, as one epoch-changing change that rides the pending epoch.

## Decided (Boss, 12:41 to 12:43 AM EDT on 2026-10-18, card taps R1 to R4, all on the recommended option)

- **R1** Dimensions in the registry plus named relations with a dimension test; plain floats in the kernels; no operator arithmetic on `Quantity` in the kernels (UN16, optional, stays unscheduled).
- **R2** Relations declared as data (name, input kinds, output kind, formula), with the float kernel and the `Quantity` wrapper generated from the one declaration.
- **R3** SI inside the kernels, other unit systems (AU, year, solar mass) generated from the registry where a kernel wants them, so 4 pi squared comes from G's dimension.
- **R4** The relations replace the hand-written laws with the UN5 sweep, file group by file group, only where the golden set is unchanged.

Order and lane check (written 11:58 PM EDT): see "Lane check" at the end of section 8. The cards in section 9 below are kept as asked.

## 1. Short answer

Yes, and the project is closer to it than it looks: it already has four partial unit tables that do not know about each other (section 3). The proposal is **one registry of units, one `Quantity` value class that holds a number as given and answers in every unit the project uses, and one declaration of each field's quantity** (class field, database column, browser). It copies the point-in-space class (`SpatialPosition3D`) in six ways (section 2) and differs from it in one: the point class derives other frames lazily, and I recommend the same here (card Q1), because converting eagerly to ten forms for a hundred million objects costs memory and time for forms nobody reads.

Size: about 12 items (section 8). Direct line saving is modest, **[E] 300 to 700 lines net**, because the ad-hoc conversions are small each; the value is that a unit can no longer be wrong silently (35 of the 127 attributes of the mass-object classes carry a unit only in their name, and 422 of 926 numeric columns do **[M]**), and that the hybrid `MassObject` serializer (E1, approved 11:48 PM EDT) gets its field table for free. It is not additive to the 800 to 1,800 lines of the class note: the same field table serves both.

## 2. What the point-in-space class does that this should copy

`physics/position.py` (714 lines) and `db/point_columns.py` (265 lines) **[M]**:

1. **One class, every body has one.** "Every object in the galaxy is one `SpatialPosition3D`" (point_columns docstring). A test (`test_point_columns.py`) lists every table that holds an object and fails when one has no point, mass, vector and epoch, so a new kind cannot ship without them. The unit class should have the same gate: a test that fails when a numeric column or field has no declared quantity.
2. **The value is kept as given, in the form given.** "The position that was last set is kept as given, in the frame it was given in, and the other frames are derived from it" (position.py docstring). A moon placed by its system offset keeps its metres exactly. For units that means the number entered in solar masses stays in solar masses, and kg is derived; never a round trip through a converted value.
3. **Every form is always in step.** Set any coordinate in any of three frames and three forms and the rest follow (`set_galactic_spherical`, `set_system_cartesian`, nine setters). For units: `Quantity(3.828e26, "W")` answers `.W`, `.kW`, `.Lsun`, `.erg_s` and the display ladder's rungs, all consistent.
4. **Derived forms are worked out on first read, not on every write.** `_sync` only marks the forms stale; `_fresh` computes them when read (PERF.46: "a body is moved several times as it is placed, and only the last position is ever read"). This is what keeps the cost down, and it is why I recommend lazy forms (Q1).
5. **A paired quantity follows its partner.** `set_mass(kg)` sets `mu = G x mass`; `set_mu` sets the mass back. The unit class generalises this to declared relations (mass and mu, radius and diameter, period and frequency).
6. **Refuses nonsense at the edge.** Not finite, negative radius, speed at or above c raise `ValueError`. The unit class checks finiteness, sign (a mass cannot be negative) and an optional range, per quantity.
7. **The storage layer knows the spellings in one place.** `point_columns.py` is "the one place that knows the spellings" (`position_x_mpc`, `center_x_pc`, `position_x_km`; mass in kg or solar); `PointLayout` names each table's columns with frame and unit, and the SQL takes its column lists from it. The unit class extends this from the four point parts to every numeric column.

What it does **not** do, and the unit class must: it covers one family (length, plus velocity, mass and epoch beside it). Nothing guards, for instance, that a luminosity column is watts, and the display ladders and the storage maps are separate tables.

## 3. What exists in the code

Four partial unit tables plus ad-hoc factors **[M]**:

| Piece | File | Lines | What it knows |
|---|---|---|---|
| Constants | `physics/constants.py` | 472 | Exact factors, with duplicates under two names: `AU_M` and `AU_TO_M`, `LIGHTYEAR_M` and `LY_TO_M`, `KM_M` and `KM_TO_M_FACTOR`, `SPEED_OF_LIGHT_M_S` and `SPEED_OF_LIGHT_KMS`; used on 517 lines in 70 files (385 more in tests). `AU_TO_KM` alone 157 lines, `SOLAR_MASS_TO_KG` 65. |
| Distance helpers | `physics/units.py` | 143 | 8 functions (ly, pc, mpc, au conversions); 27 files import it, 16 call sites. Docstring says the generator "keeps its native light-year units ... and never calls this". |
| Display ladders | `util/ladders.py` (415) + `util/format.py` (606) + `html/static/unitladder.js` (188) + `unitladders.js` (165, generated) | 1,374 | 11 quantities (length, speed, time, temperature, pressure, gravity, mass, power, flux, number density, magnetic field), each rung with symbol, factor, name, offset, preference system. 21 `format_*` wrapper functions in `format.py` (lines 219 to 447, about 230 lines). This is the best unit table in the project and is display-only. |
| Storage maps | `db/point_columns.py` | 265 | `LENGTH_M` (m, km, pc, mpc, au) and `MASS_KG` (kg, solar): two more copies of factors. |
| Position class | `physics/position.py` | 714 | Its own `SPEED_OF_LIGHT_MS`, a bare `1000.0` in `set_velocity_kms`, `length_unit_m` for the whole position. |
| Describe | `describe/model.py`, `quantity.py` | | `qty(kind, si, **extra)`: a quantity is already a dict `{"kind", "si"}`, 4 uses; `quantity.py` redeclares `_KM_M = 1000.0`. |

Ad-hoc conversions outside those **[M]**:

- Factors re-declared outside `constants.py`: `_AU_PER_LY = 63241.077` (`generation/phenomena_plausibility.py:79`), `AU_PER_LY = 63241.077` (`phenomena/asteroid_field.py:35`), `LY_PER_PC = 3.26156` (`physics/hii.py:40`), `KM_PER_PC` (`galaxy/straight_line.py:24`), `_LY_TO_AU_FALLBACK` (`web/maps/phenomenonmap.py:43`, with a try-import for a deployment that cannot import the package), `MPC_PER_PC = 1000` (`bright_stars.py:99`), `BAR_TO_PA` (`rogue_surface.py:67`), `ERG_PER_J` (`activity.py:45`), `MINUTES_PER_YEAR` (`navigation.py:269`), `_KM_PER_S_TO_AU_PER_S` (`system_tree.py:33`). The AU per ly (63241.077) and the ly per pc (3.26156) are rounded, not the exact ones: they differ from `LIGHTYEAR_M` and `PARSEC_M` by 1.3e-9 and 1.2e-6 in relative size **[C]**.
- A bare `* 1000` or `/ 1000.0` beside a unit name: 96 lines in 37 files (167 lines with a bare 1000 in 58 files, of which the 96 sit next to km, m, mpc, kms, kg, pc or au) **[M]**.
- Browser: `AU_KM = 149597870.7` declared in `orbitpositions.js`, `systemview3d.js` and `galaxysystem.js`; `KM_PER_PC` in `galaxysystem.js` **[M]**.
- Schema: 926 numeric column declarations, 422 with a unit in the name: length 169 (mpc 10, pc 60, km 82, au 10, ly 7), mass 31 (kg 14, solar 15, earth 2), time 88 (Gy 21, years 25, unix 22, hours 6, yr 11, ...), speed 68, angle 41 (deg 37, rad 4), temperature 13, power 7 **[M]**.
- The same quantity under several units in the mass-object classes: mass as `mass_kg`, `mass_solar`, `initial_mass_sol`, `primary_mass_kg`, `primary_mass_solar`; angles as `_deg` (8 of 35); time as `_years`, `_gy`, `_hours`, `_ms`; length as `_km`, `_au`, `_ly`. 35 of 127 attributes in 8 class files carry a unit suffix **[M]**.

So the unit of a value is a **naming convention**, checked by nothing. A wrong suffix, or a value in AU passed where km is read, is found only if a test happens to cover it.

## 4. The design

### 4.1 Vocabulary

- **Quantity kind** (what is measured): length, mass, time (duration), instant (a point in time, the Unix epoch), speed, temperature, angle, power, pressure, density, flux, ... The ladders already say `quantity="length"`.
- **Unit**: a named factor for one kind, with an optional offset (°C, °F), a display symbol, an ASCII key and a name.
- **Quantity** (the value class): a number in one unit, as given.

The word "quantity" is already the project's (`Ladder.quantity`, `qty()`), so I use it for the value and keep "unit" for the table row.

### 4.2 The registry

One table in one module, built from `constants.py` (which stays the home of the exact numbers):

```python
UNITS = {
  "length": Unit("m", base=True), Unit("km", 1e3), Unit("au", AU_M), Unit("ly", LIGHTYEAR_M), Unit("pc", PARSEC_M),
             Unit("mpc", ...), Unit("kpc", ...), Unit("Rsun", SOLAR_RADIUS_M), ...,
  "mass":   Unit("kg", base=True), Unit("Msun", SOLAR_MASS_TO_KG), Unit("Mearth", ...), Unit("Mjup", ...),
  "power":  Unit("W", base=True), Unit("kW", 1e3), ..., Unit("Lsun", SOLAR_LUMINOSITY), Unit("erg_s", 1e-7),
  "temperature": Unit("K", base=True), Unit("C", offset=273.15), Unit("F", factor=5/9, offset=255.372...),
  ...
}
```

Each quantity kind also records a dimension (length, mass, time, temperature exponents) so the tests can check that a compound unit's factor agrees with its parts (km/s = km / s), and the allowed range (finite; mass at least 0; speed below c). The display ladders keep their rung order, thresholds, asides and preference systems as display data, but their `symbol` and `factor` come from the registry, so a factor is stated once. The generated JS table (`build_unit_ladders.py` precedent) carries the registry to the browser.

### 4.3 The Quantity class

```python
power = Quantity("power", 3.828e26, "W")   # truth: the number and unit as given
power.W, power.kW, power.Lsun, power.erg_s # every unit the project uses, always in step
power.as_(("MW"))                          # explicit form; raises for a unit of another kind
power.set(1.0, "Lsun")                     # change the truth; every form follows
power.value, power.unit                    # the stored number, never converted
power.to_dict() -> {"kind": "power", "value": 3.828e26, "unit": "W"}   # the describe/API shape
```

- **Truth as given** (copy 2). Other forms are `value * factor_given / factor_wanted` from the truth, never from another derived form, so no drift builds up. A conversion between two units whose ratio is exactly a power of ten (km and m) stays exact in the common range.
- **Forms** (copy 3 and 4): computed on first read from the truth and cached on the instance, marked stale by `set`. `__slots__` keeps the object to about 5 words.
- **Relations** (copy 5): declared pairs, e.g. `mass <-> mu` with G, `radius <-> diameter`, `period <-> frequency`; the point class's `mass_kg`/`mu` becomes the first declared one.
- **Checks** (copy 6): finite, sign and range from the registry; wrong-kind unit refused with a clear message.
- **Offsets**: temperature in °C is absolute; a temperature *difference* is its own kind (`delta_temperature`) with no offset, because 10 K warmer is 10 °C warmer but not "283 °C". The ladders already carry an offset for °F.

### 4.4 Hot paths

Sector fill builds millions of objects. The earlier units research measured a plain divide at 0.17 to 0.5 microseconds, an `astropy` unit conversion at 4.5 to 5.6, and a `Quantity` round trip at 22 to 37 (`docs/design/architecture/units-and-number-formatting.md`, section 2 **[C]** there). So the rule: **generation maths stays on plain floats in the field's declared native unit, and the registry is consulted at the edges** (when a value is stored, loaded, shown or handed to the API). The class carries the unit and the forms; the numeric kernels (draws, physics) do not wrap every float. A field declaration gives both worlds:

```python
class StellarObject(MassObject):
    mass       = Field("mass", native="kg")        # a plain float on the instance; the unit lives in the declaration
    luminosity = Field("power", native="W")
```

`obj.mass` reads the float (fast, in kg); `obj.q.mass.Msun` or `obj.quantity("mass").Msun` builds the `Quantity` when a form other than the native one is wanted. This also gives the hybrid `MassObject` its serializer table: every field's name, kind, native unit and column come from the one declaration (E1 option 3).

## 5. How it fits the hybrid MassObject

- Every unit-valued field on the core (`mass`, `radius`, `age`, `spin`, the parent-relative state vector, galactic motion) and on each profile (stellar `luminosity`, `temperature`; planetary `surface_gravity`, `pressure`; small body `diameter`) is a `Field(kind, native)`. The nine classes' 313 slots **[M]** become one table per class family.
- **A point in space stays what it is** (every mass object has one), but it reads its length unit from the registry (`length_unit_m` becomes `length_unit="au"`), loses its own `SPEED_OF_LIGHT_MS` and bare `1000.0`, and its `mass_kg`/`mu` pair becomes a declared relation. Position truth keeps its per-frame forms; nothing in section 2 changes.
- **Ages and instants**: age is a duration (Gyr as given), `epoch_unix` is an instant. They are separate kinds so an age is never added to an epoch.
- **Dimensionless values** (fractions, ratios, counts, mixing ratios) are not quantities; the registry has a `ratio` kind only for those that need a symbol or percent display.

## 6. Storage

- **One column per value, in one declared native unit**, as in the code. Not one column per unit and not all forms: the forms are computed.
- **The native unit is chosen per column for precision, not forced to SI.** Reason **[C]**: a double has a spacing (ulp) that grows with the number; at 8 kpc from the galactic center (2.5e20 m) it is 32,768 m, so a position stored in metres from the center cannot place a 10 m body, and stored in kpc the same. Hence the existing frames and units (galaxy pc, sector mpc, system km or au) are kept, and the registry declares them.

  | Distance from the reference | Spacing of a double, metres |
  |---|---|
  | 1 AU | 3.1e-5 |
  | 100 AU | 2.0e-3 |
  | 10,000 AU | 0.25 |
  | 100,000 AU | 2.0 |
  | 1 pc | 4.0 |
  | 1 kpc | 4,096 |
  | 8 kpc | 32,768 |

  Consequence for the 10 m floor: a 10 m body can be placed to a third of its own size only within about 10,000 AU of its parent, which is also where the stability limits of a bound orbit sit (multistar `FRAGILE_AXIS_AU` 10,000). Anything farther is parent-relative in a wider frame or is a wanderer with its own galactic frame at a lower precision. This is a limit of the floating-point type, not of the unit class, and the registry makes it a number a test can read: each length column declares its reference frame and the largest distance for which a 10 m body keeps a metre.
- **Mass range**: 10 m rock about 1.05e6 kg up to about 2e40 kg is 34 orders of magnitude; a double holds it; the unit class stores the number as given, so `mass_solar` for a 10 m rock is 5e-25 in solar masses with full relative precision.
- **Column registry** (extends `point_columns.py`): `column -> (table, kind, native unit)` for every numeric column, generated from one list, replacing the unit suffix as the source of truth. The suffix stays as a readable name, and a test checks `suffix == unit key` (card Q4). `PointLayout` becomes a thin client of it.
- **No data migration is needed** for the first steps: values and units stay as stored; only declarations are added. A column that is later re-unified (two mass columns in kg and solar for one concept) is a normal migration with the generator-epoch note below.
- **Gate**: a test like `test_point_columns.py` that fails when any numeric column has no declaration, so a new table cannot ship without a unit.
- **Built (DB.76, UN4).** `db/column_units.py` holds the list: `TABLES` (table -> unit -> columns, 974 numeric columns of 63 tables), `LOCAL_UNITS` (the few units `physics/unit_registry.py` does not hold yet) and `SUFFIXES` (the unit a name's suffix stands for). `PointLayout` reads its position and mass units from it. `tests/test_column_units.py` is the gate of this section: every numeric column declared once, every declaration a column, and the suffix equal to the declared unit (one exception, `nebula_field.cell_k`, an index). Four quantities whose code states no unit are declared `unspecified`. The load and store of the generic insert (DB.74) is the next client; the dimensions and relations of section 4.2 (UN1) are not built.

### 6.1 Built (DB.77, UN10): the precision table

`db/precision.py` gives each of the 182 length columns a role (`position` 63, `separation` 42, `size` 61, `derived` 12, `relative` 4), and each position or separation a frame and a reach: `bound-orbit` (10,000 AU, the note's figure), `stability-ceiling` (2.2 million AU, the Hill sphere of a 150 solar-mass star at 8 kpc, the farthest a comet's closed orbit goes), `sector` (the half-diagonal of a 4 pc cell) and `galaxy-edge` (262,140 pc, the largest Edge the dials allow). A size keeps its relative precision at any size and declares no reach. From the unit's spacing at the reach the table gives the smallest body each table places to a tenth of its size (`floor_body_m`):

| Tables | Smallest body placed well | Why |
|---|---|---|
| planets, moons, asteroid belts, orbit nodes, scatter nodes | 2.4 m | km column at 10,000 AU (0.24 m) |
| stars, scatter members | 2.7 m | AU column at 10,000 AU |
| comets | 620 m | km column at the stability ceiling (62 m); a comet inside 10,000 AU is placed as well as a planet |
| bright stars, the phenomenon scatter | 9,200 km | mpc column in a galaxy frame at the largest Edge |
| star systems, free bodies (black holes, neutron stars, rogue planets, interstellar comets, asteroid fields, nebulae, remnants, quasars, facilities, vessels), nebula block centroids, sector path knots, container centres | 9,000 km | pc column in a galaxy frame at the largest Edge (the spacing at 8 kpc is 33 m, at 1 kpc 4 m) |

So the 10 m floor body holds where it is parent-relative inside 10,000 AU (a test, `REQUIRED_FLOOR_M`, fails when a planet or moon column cannot place one), and a body stored free in a galaxy frame is placed to a tenth of its size only above the table's figure: a rogue 10 m rock is a point to the kilometre at 1 kpc, to tens of kilometres at the Edge. The way out is the one the decision names, a parent-relative state in a wider frame for a body that joins a system (the 8.4 revision) and a lower declared precision for a wanderer; the table makes the cost a number instead of an assumption. `tests/test_precision.py` pins the figures, so moving a column to another unit or changing a reach shows up as a failing test and an edit of this table.

## 7. Risks

- **Bit-identical output.** Replacing `x * 1000.0` with a registry factor is exact, but replacing `ly_to_au(x)` with a path through km or a rounded `63241.077` by the exact one changes the last digits of generated values and so the golden files and the generator epoch. Rule: the first sweep (UN5) must produce bit-identical output on the golden set; where the exact factor differs from a rounded literal it ships as its own change with an epoch entry (one pending epoch rides until a release stamps it; epochs 12 to 18 stamped, 19 pending at 9:38 PM EDT) and is grouped with another epoch-changing item.
- **Performance** (section 4.4): kept off the hot path by design; a PERF check on the sector fill benchmark is part of UN3.
- **Offset units and differences** (section 4.3).
- **Breadth**: 517 lines in 70 files plus 385 in tests use unit constants. Sweep by file group with the owning lane, inside the host chain items, not as one change (card Q5).
- **JS and Python agree**: the browser table is generated, and a test fails when the committed file is stale, as for the ladders (`test_unit_ladders.py`).
- **Name**: `Quantity` collides with `astropy.units.Quantity` in anyone's head; astropy is no longer a runtime dependency (GEN.211); only tests (`test_astropy_constants.py`, `test_unit_ladders.py`) use it. Keep the registry free of an astropy import.

## 8. Items, for the TODO thread to number

Provisional labels UN1 to UN12, in this order. Lanes as in the unified plan (Core 4 host chain and the MassObject core, Core 1 storage, Core 5 small sweeps, Core 2 browser), but the TODO thread decides.

- **UN1 Registry.** One module with every unit the project uses, built from `constants.py`, with kind, dimension, offset, range and the declared relations; duplicate constant names removed (`AU_TO_M`, `LY_TO_M`, `KM_TO_M_FACTOR`, the copy of c in `position.py`). Tests: factors against the exact definitions; compound units agree with their parts; round trips. Core 5 or Core 4. No dependencies.
- **UN2 `Quantity` class.** Truth as given, lazy cached forms, relations, checks, `to_dict`. Depends on UN1.
- **UN3 `Field` declaration and the MassObject table.** Native unit per field; the table-driven serializer of U4a reads it. Depends on UN1 and UN2 and on U4a; Core 4. Includes the sector-fill benchmark check.
- **UN4 Column registry.** Generalise `point_columns.py` to every numeric column; load and store through it; the "no numeric column without a unit" test. Depends on UN1; Core 1; fits U11a.
- **UN5 Replace the ad-hoc conversions.** The 10 re-declared factors, the 96 bare-1000 sites, the 16 `units.py` calls, the `describe` constant; bit-identical on the golden set; in file groups. Depends on UN1; Core 5.
- **UN6 Display reads the registry.** Ladder rungs take symbol and factor from the registry; the 21 `format_*` wrappers (about 230 lines) become one `format_quantity(Quantity)` entry plus thin names where the call sites need them; UX.107 call sites move with it. Depends on UN2.
- **UN7 Browser table.** The registry generated to JS beside `unitladders.js`; the three `AU_KM` and the `KM_PER_PC` declarations removed. Core 2.
- **UN8 Position class on the registry.** `length_unit` by key, the bare 1000 and the copy of c removed, `mass`/`mu` a declared relation. Depends on UN1 and UN2.
- **UN9 Describe and API.** `qty()` builds a `Quantity`; the `{"kind","si"}` shape becomes `{"kind","value","unit"}` with `si` kept only where an API consumer needs it (no backward shims otherwise). Ties to UX.107 and UX.119.
- **UN10 Precision table.** Each length column declares its reference frame and the largest distance at which a 10 m body keeps a metre; a test reads it. Depends on UN4.
- **UN11 Docs.** `units-and-number-formatting.md` and `database-schema.md` updated; the unit rules go in the file header of the registry.
- **UN12 Sweep the mass-object classes.** The 35 unit-suffixed attributes in the 8 class files, and the same quantity under two units (mass kg and solar, angle deg and rad), unified by the hybrid, one declared native unit and the other form derived. Part of U4a; Core 4.

Estimated size **[E]**: UN1 to UN2 about 400 to 600 new lines with tests; UN5 and UN6 remove 300 to 500 lines; UN3, UN4 and UN12 fold into the savings already estimated for the MassObject core (not additive). Net direct saving 300 to 700 lines.

**Lane check.** The order above holds, with three points for the lanes:

- **UN1, UN2 have no lane conflict.** They are new modules and touch no existing file except removing the duplicate constant names in `constants.py` (UN1), which all 70 files import from; do that removal last in UN1, after UN5's sweeps, or keep the old names for one PR only if the sweep is not yet merged (no shim past that).
- **UN3, UN12 wait on U4a** (the MassObject core), which Core 4 owns after GEN.315/264 (Core 3). Until U4a exists they cannot start; UN1, UN2, UN4, UN5, UN7 and UN8 do not wait on it.
- **UN5 touches the most files** (517 lines in 70 files), so it will meet every lane's open branches. Sweep a file group only when no lane has an open PR on it; the TODO thread can sequence by file group. UN6 overlaps UX.107, which Core 4 already owns (description layer), so UX.107 and UN6 should be one owner. UN4 touches `point_columns.py` and the tables Core 1 changes in its identity and repair chains (DB.43 to 17); Core 1 should own UN4 and sequence it after its current storage item.
- None of UN1 to UN12 blocks 8.1 or GEN.301; they ride the host chain items as decided.

## 9. Cards

**Q1. When are the other units worked out?** (Boss said "immediately".)
1. Lazily on first read, and always in step (the point-in-space way). Same answer to the caller, no cost for forms nobody reads. **Recommended.**
2. Eagerly at entry, all forms stored. Literally "immediately"; about 8 stored numbers per length value and 5 per mass, and PERF.46 found that cost on positions.
3. Eagerly for the three most used forms of each kind, lazily for the rest. A middle path, with a rule to maintain.

**Q2. What does a field hold?**
1. A plain float in the field's declared native unit, with a `Quantity` built at the edges. Fast in generation; the unit lives in the declaration. **Recommended.**
2. A `Quantity` object in every field, always. Safest against a wrong unit, 20 to 100 times slower per operation in the kernels **[C]** (section 4.4).
3. Plain floats everywhere with the registry only used for display and storage. Smallest change, but the unit stays a convention inside the code.

**Q3. What is a column's stored unit?**
1. Declared per column in the registry, chosen for precision (galaxy pc, sector mpc, system km or au), as stored. **Recommended.**
2. SI for every column. Simplest rule; cannot place a 10 m body away from its parent (section 6 table).
3. SI for everything except positions, which keep their frames. A smaller rule list, but the same exception reasoning as option 1.

**Q4. What happens to the unit suffix in names (`mass_kg`, `radius_km`)?**
1. Keep the suffix and add a test that it equals the declared unit. No renames, no migrations. **Recommended.**
2. Drop the suffix from column and field names (`mass`, `radius`) and let the registry say the unit. Cleaner; a rename migration for 422 columns and the call sites.
3. Keep suffixes on columns only, none on class fields. Half of option 2's work.

**Q5. How is it rolled out?**
1. Inside the host chain items, file group by file group, each sweep bit-identical on the golden set; the registry and class first. **Recommended.**
2. One dedicated unit refactor before the host chain. Cleaner start; touches the 70 files and the same classes the MassObject core will rewrite, so some work is done twice.
3. Only for new code, existing code swept later. Least risk; the duplicates and the two rounded constants stay.

**Q6 (only if Q5 option 1): do the two rounded factors (63241.077 AU per ly, 3.26156 ly per pc) become the exact ones in the first sweep?**
1. Yes, in the first sweep, as one epoch-changing change with the next pending epoch. **Recommended.**
2. Yes, but only when a lane already changes stored values for another reason.
3. No, keep the rounded values as tuning constants (they are used for display and plausibility, not for positions).

## 10. What this does not do

It does not change what the site shows (ladder thresholds, asides and preference systems stay), does not add a dependency, does not convert generation maths to wrapped objects, and does not choose which physical constants are used (that stays `constants.py`).

## 11. Converting between kinds: derived units and named relations (added 12:15 AM EDT on 2026-10-18; decided 12:43 AM EDT)

**Boss, 12:09 AM EDT:** "Do we need any unit converters? Not between the same type of unit (length to length) but say from one type of unit to another by putting 2 or more other units together in the right order."

### 11.1 Short answer

Yes, but it is not one more converter. It is two small things the registry of section 4.2 needs: **a dimension on every unit** (so the registry knows that a watt is mass x length squared per time cubed, and that surface gravity is length per time squared), and **a short list of named relations** that put quantities together in the right order (mass and radius give surface gravity, luminosity and distance give flux). The code does this by hand about 50 times, in three different unit systems, and each hand version carries its own unit juggling. Full operator arithmetic on `Quantity` (multiply two quantities and get a third) is optional and I would keep it out of the generation kernels.

### 11.2 What the code writes by hand **[M]**

The same law appears in several files, each with its own unit handling:

| Relation (kinds involved) | Sites in `src/planetgen` | Where |
|---|---|---|
| Schwarzschild radius, 2GM/c^2 (mass -> length) | 7 lines in 6 files | `compact_remnant.py:324`, `quasar.py:89`, `multistar.py:425`, `phenomena_plausibility.py:393,519`, `stability.py:34`, `mathcheck.py:288` |
| Surface gravity, GM/r^2, then divided by Earth's (mass, length -> acceleration) | 5 lines in 4 files | `planets.py:591`, `plausibility.py:156,203`, `rogue_surface.py:135`, `db/query.py:2047` |
| Density and sphere volume (mass, length <-> density) | 11 lines in 6 files | `admin/edits.py:350-351`, `rogue.py:254`, `bulk_mass.py:20,26`, `planets.py:177,178,556,1136,1260`, `plausibility.py:156` |
| Luminosity, radius, temperature, 4 pi R^2 sigma T^4 (length, temperature <-> power) | 9 lines in 8 files | `compact_remnant.py:574`, `star.py:1407`, `phenomena_plausibility.py:419`, `stellar_evolution.py:593`, `brown_dwarf.py:74`, `pre_main_sequence.py:69`, `radiation.py:202`, `planets.py:664,714` (equilibrium temperature) |
| Kepler's third law and mean motion (length, mass -> time) | 13 lines in 8 files, in **three unit systems** | AU, years and solar masses (`kepler.py:101`, `planets.py:340`, `stability.py:49`, `companion_classes.py:55`, `multistar.py:589-594`); SI (`hydrosphere.py:215`, `formation.py:124`, `position.py:150`); AU^3/yr^2 and km/s (`state_vectors.py:180`, `system_tree.py:33`, `body_positions.py`) |
| Escape and break-up speed, vis-viva (mass, length -> speed) | 3 lines in 3 files | `star.py:765`, `spin.py:195`, `kepler.py:343` |
| Anything using G | 45 lines in 26 files | |
| 4 pi, or a unit juggle `radius_km * 1000.0` beside a radius | 38 lines in 25 files; 24 lines of the juggle | |
| Earth's gravity as a divisor | 15 lines outside `constants.py` | |

So about 50 sites write one of six laws, and the units are handled by hand at each (a `/ 1000.0` here, `constants.KM_TO_M_FACTOR` there, a division by `EARTH_GRAVITY` in a third place). Two things in the code already point the way:

- **Natural-unit shortcuts.** Orbital code works in AU, years and solar masses so that G becomes 4 pi^2 (`kepler.AU3_PER_YR2_PER_SOLAR_MASS`, `T = sqrt(a^3/M)`). That is a derived-unit conversion done by hand once, in the head of the author, and then repeated as a bare constant. A dimension-aware registry computes it: G has dimension length^3 mass^-1 time^-2, so its number in any chosen base units (SI, or AU-year-solar-mass) follows without anyone remembering 4 pi^2.
- **`physics/mathcheck.py`** (1,012 lines, 58 checks) already runs before everything else and has an "invariant" group that checks "unit conversions round-trip, the constants agree with each other" (TEST.65). A dimension check of every relation belongs in that gate.

### 11.3 What I propose

1. **Dimensions are required, not optional, in the registry (UN1).** Each unit has an exponent vector over mass, length, time, temperature (angle and counts stay dimensionless kinds). A watt is `[M1 L2 T-3]`; a compound unit's factor is checked against its parts (km/s = km / s), and the derived kinds the code uses (density, acceleration, energy, force, mass flow, flux, number density, angular momentum) are listed with their dimensions and display units. Cost: about 100 lines and a table, no runtime cost.
2. **Unit systems.** The registry can state a coherent system by naming its base units (SI; AU, year, solar mass). The numeric value of any constant (G, c, sigma) in a system is computed from its dimension and the base units, so the 4 pi^2 constant and similar shortcuts are generated and tested, not typed.
3. **Named relations in one module (`physics/relations.py`).** Each relation is declared once as data: its name, its input kinds, its output kind, and a formula. Two entry points come from the declaration: a **float kernel** in SI (fast, for generation: `surface_gravity(mass_kg, radius_m)`) and a **Quantity wrapper** for the edges (`surface_gravity(Quantity mass, Quantity radius)` returns a Quantity of kind acceleration, every unit and the Earth multiple answered). About 15 relations cover the table above plus Hill radius, Roche limit, equilibrium temperature, inverse-square flux and mass flow.
4. **A dimension test for every relation.** The test runs the formula once on probe objects that carry a dimension and checks the result has the declared output dimension; a formula written with the ordinary operators passes through unchanged, and a mistake (km where m was meant, a missing square) fails. This is the "a mismatch is refused" part of Boss's question, done in the test gate and not at run time. It needs the formulas to use `+ - * / **` (a `sqrt` goes through a `relations.sqrt` that accepts probes).
5. **Optional: operator arithmetic on `Quantity`** (`mass / volume` gives a Quantity whose kind is found by dimension, density if it matches a listed kind, else an anonymous derived unit; adding mismatched dimensions raises). Worth it at the edges (describe, API, one-off values, tests); not in the kernels, because the unit research measured a Quantity round trip at 22 to 37 microseconds against 0.17 to 0.5 for a float operation **[C]** (section 4.4).

### 11.4 Cost and saving **[E]**

- New: dimensions and systems about 100 to 150 lines; the relations module and its tests about 250 to 350 lines; operator arithmetic, if chosen, about 150 lines.
- Removed: the roughly 50 hand-written sites collapse to about 15 relation functions and calls; 60 to 120 lines and the unit juggle at each site; the three Kepler unit systems become one law with three views; the Schwarzschild, gravity, density and luminosity forms stop disagreeing (they are written once).
- Net lines: about zero to plus 150. **The saving is not lines.** It is that a law has one definition, a unit slip is found by a test, and Boss's own example (put in a watt, get every power format) extends to put in a mass and a radius and get every density and gravity format.
- Risk: moving a formula can change its last digits (an operation order changes), so every relation is introduced with the same bit-identical rule as the unit sweep (card Q5): a relation replaces a site only where the golden set is unchanged, and a site that differs in the last digit goes with an epoch-changing change.

### 11.5 Items to add to section 8 (decided: R1 to R4 tapped on the recommended options, 12:41 to 12:43 AM EDT on 2026-10-18; TODO files UN13 to UN16)

- **UN13 Dimensions and unit systems in the registry.** Extends UN1 (dimension required; the derived kinds and display units; unit systems and constant values computed from dimension). Core 5 or Core 4.
- **UN14 `physics/relations.py`.** The declared relations and the two entry points, with the dimension test in the `mathcheck` gate. Depends on UN13, UN2. Core 4.
- **UN15 Replace the hand-written relation sites** in file groups (the table of 11.2), bit-identical on the golden set. Part of the UN5 sweep, Core 5, in files not open in another lane.
- **UN16 (optional) operator arithmetic on `Quantity`** for the edges. Depends on UN14. Core 5.

### 11.6 Cards (all four tapped on the recommended option, 12:41 to 12:43 AM EDT on 2026-10-18, as the coordinator relayed) (answered, recommended option on all four; see the Decided block at the top)

**R1. How far does dimension awareness go?**
1. Dimensions in the registry and named relations with a dimension test; plain floats in the kernels; no operator arithmetic. **Recommended.** Catches the slips where they are written, costs nothing at run time.
2. Also operator arithmetic on `Quantity` at the edges (describe, API, one-off values). Adds UN16 for a nicer edge code; two ways to build a derived value.
3. Full dimension-aware arithmetic in all code, kernels included. Safest against slips; 20 to 100 times slower in the kernels **[C]** (section 4.4), so it forces a rewrite of the generation maths.
4. No dimensions: only same-kind conversion as in section 4. Smallest change; the 50 hand-written sites and the three Kepler unit systems stay.

**R2. How are relations declared?**
1. As data (name, input kinds, output kind, formula), with the float kernel and the Quantity wrapper generated from it. **Recommended.** One definition per law.
2. As hand-written SI kernels in one module, with the dimension test beside them and no wrapper. Simplest; edges convert by hand.
3. As Quantity-only functions. Clean to read; every call pays the Quantity cost.

**R3. Which unit system do the kernels use?**
1. SI everywhere inside a kernel, converted at the edges, with the AU-year-solar-mass shortcut generated from the registry where a kernel runs in it for speed or readability. **Recommended.**
2. Keep AU-year-solar-mass for orbital code and SI for the rest, as in the code, and only register both. No rewrite; the three Kepler views remain by hand.
3. One system for everything (AU-year-solar-mass). Constants are neat; km-scale and kg-scale code gets awkward numbers.

**R4. When?**
1. With the UN5 sweep and Q5's rule: a relation replaces a site only where the golden set is unchanged; the registry and the relations module first. **Recommended.**
2. A dedicated relations sweep first. Cleaner start; touches 26 files at once and meets every lane's open branches.
3. Only for new code. Least risk; the duplicates stay.

## Built: the Python text output takes a Quantity (GEN.348)

The describer holds each value in the unit it is stored in and writes the Document's `si` from the `Quantity`
(`describe/model.measured`, `describe/build._held`), so `describe/build.py`, `describe/phenomena.py` and
`describe/quantity.py` multiply by no unit factor; the registry gained `ld` (light-days) for a quasar's broad-line region. A page
measurement is `fmt.quantity_html(Quantity)` or `fmt.measured_html(unit, value)` (a dash for a value its kind refuses, never a
traceback), which the system pages, the object lists, the profile age labels, the population ages and the body radius
ratio use. Text is unchanged: 380 generated Documents rendered in four formats were identical before and after, and the
Document's numbers differ by at most 3e-16 relative (the AU to km step and the light-day product round differently).
Left for later: `dose` has no registry kind, so `describe/quantity.py` still writes mSv by a factor of 1000; the thin
`format_distance_km`-style wrappers stay (they hold no factor); the class pages' ranges already hold their numbers in the
unit they print, so they convert nothing.

The column registry (`db/column_units.py`) no longer has its own `light-day`: `broad_line_region_light_days` is declared in the registry's `ld`, so a unit lives in one place.

## Built: the registry in the browser (UX.140)

`scripts/build_unit_registry.py` writes the registry to `html/static/unitregistry.js` (every unit with its kind, symbol, names,
factor and offset, every kind with its base unit, dimension and bounds, and `unit`, `toBase`, `fromBase`, `scale` and `convert`,
the same arithmetic as `unit_registry.py`). `tests/test_unit_registry_js.py` fails when the file is stale, checks every
unit's numbers and every conversion between the units of a kind (all pairs, 14 values each) against Python through Node for exact
equality, and fails when a browser script declares a unit factor. `orbitpositions.js`, `systemview3d.js`, `galaxysystem.js`,
`orbitclock.js` and `systemreach.js` read `AU_KM`, `KM_PER_PC`, `SECONDS_PER_YEAR` and the 50 AU of the nominal system from it
(the last one was the rounded 206264.806 AU per parsec; it is the exact factor now, a 1e-9 relative change in a pixel-size
estimate). The browser sweep of the other conversions and the JavaScript `Quantity` are UX.128.

## Built: field declarations, stellar family (GEN.570, first slice)

`MassObject.FIELD_UNITS` is the one table of a family: each persisted attribute with the unit its instance holds it in, a key of the unit registry for a quantity or a word of `mass_object.PLAIN` (`text`, `ratio`, `count`, `index`, `flag`, `dex`, `unit-vector`, `id`) for a field that holds none. `SERIALIZABLE_FIELDS` is the table's names (`list(FIELD_UNITS)`), so the persisted list and the units cannot disagree, and the shared groups (galactic motion, spin, star activity, the disk descriptor) are declared once in their modules (`GALACTIC_MOTION_UNITS`, `spin.SPIN_UNITS`, `activity.STAR_ACTIVITY_UNITS`, `pre_main_sequence.DISK_UNITS`) and spread into each family's table. An instance still holds a plain float; `obj.quantity("mass")` builds the `Quantity` in the declared native unit (a pair such as the habitable zone as a tuple), so `star.quantity("mass").Msun` works and nothing in a constructor changed: the sector-fill path does no more work than before, the class-level dict is read only at the edges.

Done for `StellarObject` (`Star`), `BlackHole`, `NeutronStar` and `BinaryStarProxy`, and for the small bodies and the free planet (`Comet`, `InterstellarComet`, `AsteroidField`, `RoguePlanet`); the rogue surface fields are declared in `rogue_surface.ROGUE_SURFACE_UNITS`. `tests/test_field_units.py` fails when a persisted field has no declaration, when a declaration is neither a registry unit nor a plain word, when a name's unit suffix and its declaration disagree (card Q4), and when a declared unit measures something other than the column the attribute is stored in (the `Field` lists of `store.py` against `db/column_units.py`); the family still to declare (`Planet`, with the atmosphere, magnetism, hydrosphere, radiation and habitability groups, several in units the registry does not hold yet) is listed in the test and must leave the list when it declares. Declaring found the quantities held under two units that GEN.575 unifies: a remnant's `mass_solar` beside the inherited `mass` in kg, its `event_horizon_radius_km` beside `radius`. The registry gained the rate units `1/s`, `1/yr` and `1/Gy` (frequency), which the column registry no longer holds locally.
