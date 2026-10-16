# Which later phenomena are scattered, backfilled or fill-only

Boss (2026-10-16, 12:30 AM EDT), to Core 3: say how each new phenomenon of goal issue #1677 applies to the scatter and back-scatter phases, and what stays sector-fill-only (like comets and rogue planets) and what gets scattered. This note is the answer (GEN.277). It is a recommendation: nothing here is built, and nothing is built from it until Boss approves.

Informs: GEN.113, GEN.149, GEN.151, GEN.153 to GEN.155, GEN.189, GEN.244 to GEN.246, GEN.254 to GEN.257, and the rule set of GEN.235. Read with [phenomena-catalogue-and-fill.md](phenomena-catalogue-and-fill.md) (the kind list and fill stages), [multi-star-systems.md](multi-star-systems.md) (companions at scatter time), [intergalactic-and-halo-contents.md](intergalactic-and-halo-contents.md) (the sector profile) and [anomalies.md](anomalies.md) (GEN.113).

Evidence tags as in the catalogue: **[S]** seen in the repository or a search summary, **[C]** computed with the repository's code, **[R]** recalled and unconfirmed, **[I]** inferred.

## Summary

- **Four ways an object can come about**, and every item gets exactly one: a **scatter row** (a point, a seed and a kind placed after `plan`, with stellar companions stored at scatter time and no planets), the **backfill** (the sector-by-sector top-up around filled space: the mass rings of GEN.187 and the below-cut remnant draws of GEN.168), **sector fill only** (made when a sector is filled, from a roll on its own seed, like comets and rogue planets), and **derived** (no object: a value or an overlay computed from other things).
- **Result for the fourteen items** (table below): two items stay scatter rows and change only what the row stores (GEN.149, GEN.151); one object becomes a new scatter class (magnetars, GEN.153); two are components that feed the existing passes (GEN.244 and, by derivation, GEN.254); one is a state of stars the passes already place (GEN.246); one is an attribute drawn at build, fill-only (GEN.245); five are derived values or overlays (GEN.154, GEN.189, GEN.255, GEN.256, GEN.257); one is a documented decision (GEN.155) and one is the analysis behind the rest (GEN.113).
- **No new column on `phenomenon_scatter` or `bright_stars`.** Magnetars and remnant classes ride on the existing `subtype`; the planetary-nebula row stores a mass and a luminosity in the existing `mass_solar` and `map_luminosity_sol`. The only new columns are the disk descriptor on `stars` (shared with GEN.242) and the age and profile columns of the component table (GEN.237), both Core 1's.
- **One real call for Boss** (the rest stand as defaults): whether magnetars are scatter rows (recommended) or fill-time draws only (section 3).

## 1. The phases, as the code has them

| Phase | What it does | Where [S] |
|---|---|---|
| **Scatter** | after `plan`, object first, one Poisson count per shell from a certified density bound: bright stars (9,000 L_sun and up), the mass pass (born at 14 solar masses and up), black holes at 14 and up, intermediate-mass holes, planetary nebulae, supernova remnants, hypervelocity stars, the nucleus. A row is a point, a seed, a kind, a `subtype`, a velocity, a mass and a map luminosity. Stellar companions are stored with it (GEN.260); planets, moons and belts never are | `bright_stars.py`, `phenomenon_scatter.py` (`SCATTER_CLASSES`, `PHENOMENON_SCATTER_COLUMNS`) |
| **Backfill** | around generated sectors only. The bright-star backfill adds the stars the scatter left out by initial mass in rings (1, 2, 5 and 8 solar masses, GEN.187) and the luminosity bands (GEN.44); the phenomena backfill is the below-cut draw: the neutron stars and black holes under 14 solar masses, Poisson per sector from the sector's own stream (GEN.168). Neither stores anything before the sector is filled | `bright_stars.backfill_mass_cells`, `phenomenon_scatter.below_cut_draws` |
| **Fill only** | a sector rolls it when it is filled, from its own seed, for a kind that is common, small or local: rogue planets, brown dwarfs, comets, asteroid fields, hosted nebulae, the descriptors of a star | `generate_sector_phenomena`, `add_star_hosted_nebulae` |
| **Derived** | no object: computed from a position, a mass or other rows when asked | `activity`, the pages, the map overlays |

**What decides the phase.** Four questions, in this order:
1. *Is it an object at all?* A flag, a score or an overlay is derived.
2. *Is it so rare that a fill, drawing only inside the sectors that are filled, could never make it visible to a neighbour?* If something else must see it from outside its sector (a hazard score, a map before the fill, a nebula that needs it), it is a scatter row. The count does not matter for cost: every kind in this note is under 1e5 rows against 6e7 already stored [C: multi-star section 2].
3. *Is it a state or attribute of a star that already exists?* Then it lives in the star's record and follows whichever pass placed the star.
4. *Does a field or component decide where it is, not the density?* Then it is a component or a seeded field that the passes read.

Everything else (common, local, dense) is fill only, as comets are.

## 2. The decisions, item by item

| Item | Phenomenon | Phase | Why (observed rate and locality) | What the row or build stores |
|---|---|---|---|---|
| GEN.149 | planetary-nebula central star | **scatter row, unchanged kind** | about 20,000 in the Galaxy [S]; the row exists and the star is made with the nebula | the central star as member 0 (mass 0.5 to 0.7, median 0.59), its 20% close companion at scatter time (multi-star 4.4); `mass_solar` = that mass, `map_luminosity_sol` = the star's 1e2 to 1e4 L_sun (a main-sequence law would give 0.2 L_sun); no new column |
| GEN.151 | supernova remnant size | **scatter row, unchanged kind** | 400 to 3,000 visible [S]; 3.2e3 rows today | `subtype` = the remnant class letter (R to W) drawn at scatter time, so the map can size it before the fill (the class fixes the age window and the ambient density range); the exact radius still comes from the Sedov-Taylor law at build, from the row's seed; no new column |
| GEN.153 | magnetar | **new scatter class** of `neutron-star` (recommended; section 3) | about 100 active (30 to 600) [S: anomalies.md]; 1e-7 of neutron stars, so a fill would almost never put one in a given neighbourhood | `subtype` = `magnetar`, `mass_solar` drawn, young-population weight (below); B, spin and age at build from the seed. A **pulsar** needs no row: it is an attribute of every neutron star, drawn at build from its age (built) |
| GEN.154 | Einstein radius | **derived** | a number from mass and distance | nothing; a scatter row's `mass_solar` already gives it for an unfilled black hole or neutron star |
| GEN.155 | nuclear cluster | **not built** (decision recorded) | one object; the Sgr A* sector holds a star every 0.005 pc and cannot be filled star by star | nothing |
| GEN.189 | gamma-ray and nuclear-activity ozone loss | **derived** | a rate per Gyr from galactocentric radius alone, so scatter, backfill and fill give the same value for the same place | the existing `lethal_event_rate_per_gyr` on each star; no new column |
| GEN.244 | young open clusters, OB associations | **component** (GEN.237's table), feeding the existing passes | about 3,000 catalogued, tens of thousands in all [R]; local, defined by a profile, not by density | a component row (position, age 1 to 1,000 Myr, Plummer scale, [Fe/H], mass); **no scatter row per star**: the cluster adds its density to the cell's `lambda` and 100% young to its `mix`, so the mass pass and the bright pass place its O and B stars with no change and the backfill rings (1 to 8 solar masses) and the fill draw the rest. Trapezia and runaway origins ride on the member tree (GEN.273, GEN.274) |
| GEN.245 | debris disk (A and F stars) | **fill only, an attribute** | about 20% of A and F stars [R]; A and F stars are 1.4 to 2.1 and 1.0 to 1.4 solar masses, so none is pre-placed (the mass pass is 14 and up, the bright floor 9,000 L_sun) | a descriptor drawn at build from the star's seed, the same in a fill and in a backfill ring; the `disk` columns on `stars` shared with GEN.242 (disk class, mass, outer radius) |
| GEN.246 | luminous blue variable, hypergiant | **state of scattered stars**, no new kind | tens in the Galaxy [R]; all at least 25 solar masses and 1e5 L_sun, so the mass pass and the bright pass already place every one | the star's record at its age (a short window before Wolf-Rayet or explosion); the **ejecta nebula** is a hosted nebula made at fill, host first (like a Wolf-Rayet ring, GEN.243), not a scatter row. Backfill: none (the rings stop at 8 solar masses) |
| GEN.254 | superbubble field | **derived field** (seeded, like `nebula_field`) | superbubbles fill about 12% of the disk [S]; the Local Bubble is 100 to 200 pc | nothing stored; the profile function reads the field, so the cloud pass and the remnant sizes see the cavity (section 4) |
| GEN.255 | jets | **derived at render** | each needs a spin axis (GEN.104) and a state | nothing stored. Herbig-Haro jets: protostars in class Q cores, which are fill (GEN.241). Microquasar jets: a scattered black hole or neutron star whose stored member is a Roche-lobe donor in the accreting state (GEN.269); the map can draw the jet before the fill from that state. Pulsar jets: `pulsar_type` at build |
| GEN.256 | radiation hazard score | **derived**, from neighbouring rows | hot stars, Wolf-Rayet stars, remnants and magnetars within a few hundred parsecs, and the 10 pc magnetar flag | nothing stored; reads scatter rows through the cell-key region query (stage C: DB.34, GEN.276) plus the objects of filled sectors. This is why magnetars must be rows (section 3) |
| GEN.257 | Fermi bubbles | **derived** (overlay) | two lobes of about 8 kpc, no objects | nothing |
| GEN.113 | the anomaly documents | **analysis** | `anomalies.md` is the list (Tier 1 to 3, the errors found, the 20 ly to 4 pc conversion) | the Tier 3 fiction kinds, if ever wanted, are **fill only**, one `lore-anomaly` kind with `K_lore` 0: fictional objects have no observed rate to scatter by and must not show on a map before someone asks |

## 3. The one call: magnetars as scatter rows

**Recommendation (default): scatter them as a class of `neutron-star`**, beside the intermediate-mass black holes (`SCATTER_CLASSES` gains `("neutron-star", "magnetar")` last, so earlier classes keep their draws).

- **Rate.** About 100 active in the Galaxy (birth 1 per century times 1e4 yr [S, R]; range 30 to 600). That is about 3.4e-10 per star (100 over 2.95e11 systems) [C]. About 100 rows against 6.1e7 in the scatter.
- **Placement.** The young population with arm contrast, not total density: a new gate `young stars` in the registry (the gates are `stellar density`, `nebula field`, `hot stars`, `none`), reusing the thickening of `remnant_distribution.placement_factor` with a young scale height (a magnetar is under 1e4 years old and has not drifted).
- **Why a row and not a fill draw.** The fill's neutron-star draw would give a given sector a magnetar with probability about 4e-9; the Galaxy would hold about 50 in total but no neighbour could ever see one from an unfilled sector. A hazard score within 10 pc of an active magnetar (GEN.256) and a map that shows them before the fill both need the object to exist before its sector is filled. The same argument already places supernova remnants, Wolf-Rayet stars and black holes in the scatter.
- **No double count.** The neutron-star fill draw must stop producing magnetars itself: the built age lottery of GEN.153 (a neutron star younger than 1e4 yr is a magnetar with 10% chance) then applies only to the **cores of young supernova remnants**, where the age is the remnant's, and the free draw's age floor rises to 1e4 yr (`COMPACT_REMNANT_AGE_RANGE_GY` lower bound 1e-5 Gy). Both changes are in the follow-up item below, not in GEN.153.
- **Cost.** About 100 rows; the build reads `subtype == "magnetar"` and forces the magnetar branch with the row's seed.

**Built (GEN.278).** `("neutron-star", "magnetar")` is the last of `SCATTER_CLASSES`, placed by the registry gate `young stars` (`registry.SCATTER_CLASS_GATES`): a band's angle bins are weighted by the young population's density with its arm contrast (`density.population_densities_array`, now accepting one height per band), so none sit in the halo, a cluster or the thick disk, and the remnant thickening of ordinary neutron stars does not apply. The rate is `MAGNETAR_PER_NEUTRON_STAR` (1e-7) of the neutron stars' per star, divided by the young population's 0.48% of all stars in a Milky Way-sized galaxy (`MAGNETAR_YOUNG_STAR_SHARE`), which gives about 150 rows over a full scatter against the 100 asked. A magnetar ignores the mass cut (the cut leaves the ordinary neutron stars below it to the fill, and a magnetar is never drawn by the fill: `share_above` is 1). `COMPACT_REMNANT_AGE_RANGE_GY` now starts at 1e-5 Gy, so a free-standing neutron star is never young enough to be one; the cores of young supernova remnants keep the 10% birth lottery. The row builds `NeutronStar(magnetar=True)` with the row's mass and an age in `MAGNETAR_AGE_RANGE_GY`. No migration: `subtype` holds `magnetar`. The Galaxy Map lists them with the neutron stars as "Magnetar (uncharted)".

**Alternative (fill only):** keep the age lottery alone. Cheaper by one item, correct in expectation over filled sectors, but the 10 pc flag and the map can see only magnetars in filled sectors, and a filled galaxy would show about 50 where the sky has 100. Boss's call; if no answer, the recommendation stands.

## 4. What this does to the one rule set and the sector profile (GEN.235)

The profile function of [intergalactic-and-halo-contents.md](intergalactic-and-halo-contents.md) section 5.3 stays a pure function of the shape, the component table and the address. The later phenomena add terms; none adds a rule.

| Profile field | Gains | From |
|---|---|---|
| `kind_factors` | `magnetar`: young share times arm contrast, 0 in the halo, clusters and old satellites; `supernova-remnant`: unchanged (gas fraction), now also times 1 inside a superbubble | GEN.153, GEN.254 |
| `mix` | the young-cluster components contribute 100% `young` with their age window | GEN.244 |
| `gas` | the ambient density in cm^-3 inside a superbubble cavity (0.001 to 0.01 in place of the disk's 1; the Local Bubble is about 0.005 [R]); the molecular-cloud keep rate is 0 inside it | GEN.254 |
| `lambda` | the Plummer density of each cluster inside its tidal sphere | GEN.244 |
| a descriptor rule, not a field | A and F stars roll a debris disk (about 20%, [R]); stars under 10 Myr roll a protoplanetary disk (GEN.242) | GEN.245 |

Two consequences worth stating:
- **Remnant size depends on the cavity.** The Sedov-Taylor radius goes as `n^(-1/5)` and the snowplow cap as `n^(-3/7)` (GEN.151). A remnant in a cavity of 0.005 cm^-3 reaches 45.7 ly times 0.005^(-3/7) = 440 ly (135 pc) before the cap, against 46 ly at n = 1. A scatter row's class letter fixes the density range, so the map's typical radius uses the class range; the exact radius at build uses the profile's ambient density when the position is inside a cavity.
- **A seeded field, not a derived one.** The superbubble field is seeded like `nebula_field` (cell 50 pc, radius 50 to 500 pc, about 12% filling) rather than derived from the OB-association components. Deriving it from GEN.244's clusters would be more physical (a superbubble is the work of an association's supernovae) but would make GEN.254 wait on GEN.237 and GEN.241 instead of GEN.235. Default: seeded; the Local Bubble anchors it.

**Built (GEN.254).** `galaxy/superbubble_field.py`: the plane is cut into 250 pc columns (not 50: a bubble reaches 500 pc, so a point would read 8,000 columns), each column's bubbles drawn from its own seed (`superbubble-cell`), a Poisson count whose mean puts 12% of the plane at the solar circle inside a bubble and follows the molecular gas (`nebula_field.gas_factor`) elsewhere. A bubble has a centre in its column at a height of about 60 pc, a log-uniform radius from 50 to 500 pc (the Local Bubble's 100 to 200 pc among them), and 0.001 to 0.01 atoms per cm3 of ionised gas at a million kelvin. Nothing is stored. `nebula_field.cloud_kept` loses a cloud whose centre is inside one (the draw and the object IDs still count it; `--redo-nebula-field` removes the ones already stored). `SectorProfile.gas.superbubble` names the cavity a place is in and its ambient gas replaces the disk's (`profile_at(..., galaxy_seed=)`). `GET /api/galaxy/superbubbles` lists the cavities for the Galaxy Map, which draws each as three faint circles; `GET /api/sectors/<id>` carries `superbubble`, and the sector page says so under its map. The remnant size at build (GEN.151) and the remnant factor of `kind_factors` are not yet read from it.

## 5. The backfill, case by case

| Backfill pass | Effect of the later phenomena |
|---|---|
| **Mass rings (1, 2, 5, 8 solar masses; GEN.187)** | stars in a young cluster component come with it: the cell's `lambda` and `mix` include the component, so the ring draws cluster members like any others and the scatter's own pre-placed share is subtracted once. A and F stars in the rings roll their debris disk at build. No LBV or hypergiant is possible (25 solar masses and up are all pre-placed). |
| **Luminosity bands (GEN.44)** | unchanged. The planetary-nebula central stars (1e2 to 1e4 L_sun) are scatter rows of their own, so they are not in the bright-star bands either way (no double count; the rows with 9,000 L_sun and up that are central stars stay under the planetary-nebula class). |
| **Below-cut remnants (GEN.168)** | neutron stars and black holes under 14 solar masses are drawn per sector from the sector's stream. Pulsars come out by age. Magnetars come out only as the cores of young remnants (or not at all if the fill-only alternative is chosen). Intermediate-mass holes and remnants stay as they are. |
| **Phenomenon passes at fill** | superbubble cavities zero the molecular clouds drawn in them; jets and the hazard score read what the other passes made. |

## 6. Registry rows (GEN.250) these decisions imply

One row each, the columns of [phenomena-catalogue-and-fill.md](phenomena-catalogue-and-fill.md) section 6. Kinds and subtypes not listed are unchanged.

| kind (subtype) | family | mode | scatter_pass | gate | system_model | carries | notes |
|---|---|---|---|---|---|---|---|
| `neutron-star` (`magnetar`) | compact-remnant | scatter | `band-density` | `young stars` (new gate) | none | none | `id_code` is the neutron star's; subtype only |
| `supernova-remnant` (class letter) | nebula | scatter | `band-density` | `nebula field` | none | its core, at build | the subtype becomes the class letter instead of `None` |
| `planetary-nebula` | nebula | scatter | `band-density` | `stellar density` | central-star | its white dwarf and its companion | row stores the central star's mass and luminosity |
| `young-cluster` | star | component | none | none | none | its members through the passes | GEN.237's table, kind value added |
| `superbubble` | nebula | field | none | none | none | none | seeded from the galaxy seed, like `nebula_field` |
| `debris-disk` | star | flag | none | none | none | none | attribute of A and F stars |
| `lbv`, `hypergiant` | star | scatter | `bright-star` (existing passes) | `hot stars` | star-system | its ejecta nebula, host first | stellar states, not rows of their own |
| `lore-anomaly` | small-body | fill | none | none | none | none | `K_lore` 0, off unless enabled |

The completeness test of the registry then covers these: every row has a builder, a rate or a field, and a page class.

## 7. Columns and migrations (for Core 1)

- `phenomenon_scatter`: **none.** `subtype` carries `magnetar` and the remnant class letter; `mass_solar` and `map_luminosity_sol` carry the planetary-nebula star.
- `bright_stars`: **none** (LBV and hypergiant are states of rows that exist).
- `stars` (and the members of `scatter_members`, GEN.260): the disk descriptor shared by GEN.242 and GEN.245, three columns (disk class, mass in Earth masses, outer radius in AU), nullable. Core 2 owns GEN.242; this is one set of columns, added once, in whichever of the two migrations comes first.
- Component table (GEN.237): a `kind` value for `young-cluster`, age in Myr, Plummer scale radius, and a flag `hosts_nebula` for those under 10 Myr. Core 1 owns the migration; if GEN.237's table already carries age and a profile kind, no column is new.
- `COMPACT_REMNANT_AGE_RANGE_GY` and the `SCATTER_CLASSES` tuple are settings and code, not schema.

## 8. Changes to the TODO items (for the TODO thread)

| Item | Change |
|---|---|
| GEN.149 | add: the scatter row stores the central star's mass and map luminosity; the 20% close companion goes into the member tree at scatter time (GEN.261) |
| GEN.151 | add: the scatter row's `subtype` is the remnant class letter; the build's ambient density is the profile's inside a superbubble |
| GEN.153 | unchanged as built (magnetar by age, pulsar fraction by age); **new item** below does the scatter class |
| **new** | **Magnetar scatter class** (Core 3, after GEN.250, GEN.235 and GEN.153): `("neutron-star", "magnetar")` in `SCATTER_CLASSES`, the `young stars` gate, the age floor of the free neutron-star draw, the build forcing the magnetar branch; test: about 100 rows over a full scatter, none in the halo |
| GEN.244 | add: the component feeds `lambda` and `mix` only; no per-star scatter rows; the backfill reads the component term; trapezia and runaway origins come from the member tree |
| GEN.245 | add: fill-only attribute, shares the `stars` disk columns with GEN.242 |
| GEN.246 | add: states of stars the passes place; the ejecta nebula is hosted at fill; no scatter work |
| GEN.254 | add: seeded field, read by the profile; remnant sizes and the cloud keep rate see it |
| GEN.255 | add: derived at render; microquasar jets read the stored accreting member |
| GEN.256 | add: reads scatter rows through stage C; needs the magnetar rows to see the 10 pc flag before the fill |
| GEN.257, GEN.154, GEN.189, GEN.155 | no change; derived or decided |
| GEN.113 | closed by `anomalies.md` and this note; the Tier 3 fiction kinds are fill only if ever built |

## 9. Decisions and defaults

1. **Magnetars are scatter rows** (section 3). Alternative: fill only. Needs Boss's word; the default stands without it.
2. **Superbubbles are a seeded field**, not derived from the clusters (section 4).
3. **The planetary-nebula row stores the central star's mass and luminosity** (section 2), so the map agrees with the sector page, as MAP.165 did for black holes and neutron stars.
4. **Young clusters are components, not scatter rows of their stars**, so a cluster's star is never counted twice (section 2).
5. **Fiction anomalies, if ever wanted, are fill only** with `K_lore` 0.

## Sources

The folded notes' own sources. Observed counts are as in [phenomena-catalogue-and-fill.md](phenomena-catalogue-and-fill.md) and [anomalies.md](anomalies.md); figures marked [C] are computed with the repository's code or its stated formulas. The superbubble ambient density (0.005 cm^-3) is recalled, not confirmed.
