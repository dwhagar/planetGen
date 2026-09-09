# Galaxy Coordinate System — Design Proposal

**Status:** proposal, not implemented. This document is the design pass
`docs/TODO.md` Phase 4 calls for ("needs its own design pass") before
`galaxyGen.py` is written. No generation code changes here — the schema
migration below is written out concretely so it can be reviewed, but it has
not been applied to `stellarObjects/schema.sql`.

**Scope:** the galaxy-scale coordinate system and the `sectors` schema
changes it needs. Explicitly **out of scope**: `galaxyGen.py` itself,
batch-generation performance, and any model of *where in the galaxy stars
are actually dense* (spiral arms, bulge, metallicity gradients). Where
those touch this design, they're flagged as open questions rather than
decided here — see §7.

## 0. Decisions already made (recap)

Per the brief, these are fixed and not relitigated below:

- Sectors are arranged **radially** from the galactic center, not on a flat
  grid/lattice.
- A sector's galactic position is defined by **spherical coordinates
  measured from the galactic center**, but **stored/exposed as xyz**
  (the sector's center point in a galaxy-frame Cartesian system).
- Sector-*internal* geometry stays exactly as it is today: milliparsecs,
  relative to the sector's own cubic center (`star_systems.position_x/y/z_mpc`,
  `sectors.edge_mpc`). This design only adds a *second*, larger-scale
  coordinate — it does not touch the first.

## 1. Coordinate system definition

### Origin and axes

- **Origin**: the galactic center, `(0, 0, 0)`.
- **Axes**: right-handed Cartesian. `X`/`Y` span the galactic plane; `+Z`
  points toward galactic north (out of the disk). This is a generation-space
  convention for this project only — it does not need to match any real
  astronomical reference frame (galactic longitude/latitude, J2000, etc.),
  since nothing here is meant to correspond to real sky coordinates.

### Spherical convention

Reuses the same `(theta, phi)` naming `spaceSector.py`'s own
`_random_unit_direction` already uses — **not** the ISO physics convention,
which swaps the two names. Kept for internal consistency with existing
code, since anyone extending this into `galaxyGen.py` will already be
reading that function:

- `r` — radial distance from the galactic center, in **parsecs** (see §2).
- `theta` — azimuthal angle in the `X`/`Y` plane from `+X`, `[0, 2*pi)`.
- `phi` — polar angle from `+Z`, `[0, pi]` (0 = galactic north pole, `pi/2`
  = in-plane, `pi` = galactic south pole).

Conversion to Cartesian (identical form to `_random_unit_direction`, scaled
by `r`):

```
x = r * sin(phi) * cos(theta)
y = r * sin(phi) * sin(theta)
z = r * cos(phi)
```

The inverse (needed anywhere code wants "how far out, and how far off the
plane, is this sector" from a stored xyz row):

```
r     = sqrt(x^2 + y^2 + z^2)
theta = atan2(y, x)
phi   = acos(z / r)                 # undefined at r == 0; treat as phi = 0
```

### What's stored per sector

Per the fixed decision, **xyz is the persisted form** — `r`/`theta`/`phi`
are a derivation used by the placement algorithm (§3) at generation time,
not separate columns, with one exception: `r` (as `galactic_radius_pc`) is
worth persisting anyway, for the same reason `star_systems.quadrant` is
persisted even though it's derivable from `position_x/y/z_mpc` — "sectors
within radius R of the core" is an obviously common query, and a stored,
indexed column serves it directly instead of needing `sqrt(x*x+y*y+z*z)` in
every `WHERE` clause (which SQLite can't use an index for). `theta`/`phi`
have no equivalent query use, so they're not stored. See §4 for the exact
new columns.

## 2. Unit choice: parsecs (pc)

**Recommendation: store galaxy-scale distances in parsecs (`_pc` suffix
columns), not milliparsecs and not light-years.**

### Why not milliparsecs

`schema.sql`'s own header comment gives the rule that put `mpc` on
sector-internal columns in the first place: km gets unwieldy at that scale.
The same argument now cuts the other way. The Milky Way's disk radius is
commonly cited around 15 kpc (15,000 pc). In milliparsecs that's:

```
15,000 pc * 1000 mpc/pc = 15,000,000 mpc      -- and that's the radius alone
```

Two sectors on opposite sides of the galaxy would be up to **30,000,000
mpc** apart — an 8-digit number for what should be this schema's *biggest*
distance, while `sectors.edge_mpc` (~3,526) is one of its smallest. That's
the exact "ugly number" problem the brief calls out, just moved up one
level.

### Why not light-years

Light-years are the generator's own *native* unit for sector-internal
geometry (`SpaceSector.edge_ly`, `hill_radius_ly`, everything in
`spaceSector.py` before it gets converted to `mpc` at the database
boundary). It's tempting to reuse it here too, and the numbers are at least
plausible — 15,000 pc is about 48,923 ly — but `ly` and `mpc` are
**different unit families**: converting between them needs the existing
`ly_to_milliparsecs`/`milliparsecs_to_ly` helpers, which round-trip through
AU (`LY_TO_AU`, `AU_PER_MILLIPARSEC`). Every calculation that relates a
sector's galaxy-scale position to its own internal geometry — "how many
sector-widths from the core is this shell", "does this sector's cube reach
past the next shell in" — would need that AU round-trip for no reason.

### Why parsecs

A milliparsec *is* a parsec, scaled by `1/1000` — the existing convention
already anchors sector-internal geometry to the parsec family, just at
`10^-3` scale. Parsecs continue that same family at `10^0` scale:

```
edge_pc = edge_mpc / 1000          # exact, no AU round-trip
```

Worked example, `DEFAULT_SECTOR_EDGE_LY = 11.5`:

| Quantity | Value |
|---|---|
| `edge_ly` (generator-native) | 11.5 ly |
| `edge_mpc` (stored today) | ~3,526 mpc |
| `edge_pc` (proposed) | ~3.526 pc |
| Galaxy radius, same unit family | 15,000 pc |
| Galaxy radius / sector edge | ~4,254 — a plain, dimensionless ratio, no unit-family crossing |

Compare the same ratio computed via `ly`: `48,923 ly / 11.5 ly` — same
answer (~4,254), but only because you happened to convert both sides back
to the same family first; nothing about the stored columns (`mpc` vs. `ly`)
would have let you divide them directly.

### Conversion helpers needed (`stellarObjects/utils.py`)

Following the existing `ly_to_milliparsecs`/`milliparsecs_to_ly` pattern
exactly:

```
mpc_to_pc(mpc)   -> mpc / 1000                          # exact, trivial
pc_to_mpc(pc)    -> pc * 1000                           # exact, trivial
pc_to_ly(pc)     -> pc * (AU_PER_PARSEC / LY_TO_AU)     # ~3.2616 ly/pc, for display/prose
ly_to_pc(ly)     -> ly * (LY_TO_AU / AU_PER_PARSEC)     # inverse
```

`AU_PER_PARSEC` and `LY_TO_AU` already exist in `physical_constants.py`
(lines 73 and 84) — no new physical constant is needed, only the four thin
wrapper functions above, mirroring `ly_to_milliparsecs`'s own
docstring-and-one-line-body shape. `pc_to_mpc`/`mpc_to_pc` are exact
(power-of-ten scaling, no floating-point drift beyond a single multiply);
`pc_to_ly`/`ly_to_pc` exist purely for human-readable prose (e.g. rendering
a sector's summary as "~48,923 ly from the galactic core") the same way
`table_*` columns keep `ly`/`AU` display strings alongside `km`-stored raw
values elsewhere in this schema.

## 3. Radial sector-tiling scheme

### The geometric problem, stated plainly

A cube cannot tile the surface of a sphere without gaps or overlaps — this
isn't a shortcut being avoided, it's not possible in general (the reason
there's no "cubical tiling of S²," unlike the flat plane, which cubes-as-
squares tile perfectly). Any radial, shell-based layout of cubic sectors
has to either accept small seams at shell boundaries, or abandon cubes.
This proposal keeps cubes (§0 fixes `edge_mpc` as the sector shape, and
changing it is out of scope) and **explicitly accepts small gaps/overlaps
between adjacent sector cubes** as the practical resolution — justified in
"Why this is fine" below, not swept under the rug.

### Shells

Divide radial distance into shells of thickness equal to the sector edge
length, so cubes stack outward along the radial direction without gaps
*in that one direction* — the only direction cubes can tile exactly, since
adjacent shells share a flat spherical-ish boundary the way stacked plates
do:

```
shell k spans radius [k * edge_pc, (k+1) * edge_pc)     for k = 0, 1, 2, ...
nominal shell radius:  r_k = (k + 0.5) * edge_pc
```

### How many sectors per shell

Approximate each sector's footprint on the shell's spherical surface as
`edge_pc^2` (its cross-section perpendicular to the radial direction), and
divide it into the shell's surface area:

```
N_k = round(4 * pi * r_k^2 / edge_pc^2)
    = round(4 * pi * (k + 0.5)^2)          # edge_pc cancels — count depends only on k
```

This is a clean, dimensionless formula: the count per shell doesn't depend
on the absolute sector size at all, only on how many shells out you are.
Worked values:

| Shell `k` | `r_k` (sector-edges out) | `N_k` (sectors in this shell) |
|---|---|---|
| 0 | 0.5 | 3 |
| 1 | 1.5 | 28 |
| 2 | 2.5 | 79 |
| 3 | 3.5 | 154 |
| 5 | 5.5 | 380 |
| 10 | 10.5 | 1,385 |
| 100 | 100.5 | 126,923 |
| 4,253 (outermost shell reaching 15,000 pc at `edge_pc ~= 3.526`) | 4,253.5 | ~227,000,000 |

Growth is quadratic in `k` — expected, since a shell's surface area grows
with `r^2` while its thickness (and hence sector footprint) stays fixed.
Summed across all shells out to the galaxy's edge, total addressable
sectors come out close to `galaxy volume / sector volume`:

```
(4/3) * pi * 15,000^3 pc^3  /  3.526^3 pc^3   ~=   3.2 x 10^11 sectors
```

**~320 billion addressable sectors** if the whole galactic volume were
treated as a uniform sphere at this sector size. This number matters for
§7 (eager vs. lazy generation) — it's not a count anyone should try to
generate up front.

### Placing sectors within a shell

Within shell `k`, distribute its `N_k` sector centers evenly across the
sphere of radius `r_k` using a **Fibonacci (golden-angle) sphere
sequence** — a standard, deterministic way to place N points roughly
evenly over a sphere with no pole-clustering, using the same "uniform
area, not uniform angle" concern `_random_unit_direction`'s own docstring
already flags for a single random direction, just made deterministic and
extended to N points instead of 1 random one:

```
for i in 0 .. N_k - 1:
    phi_i   = acos(1 - 2 * (i + 0.5) / N_k)      # polar angle, pole-to-pole
    theta_i = (2 * pi * i / golden_ratio) mod (2 * pi)   # golden-angle azimuth step (~137.5 deg)
    center_i = r_k * (sin(phi_i)*cos(theta_i), sin(phi_i)*sin(theta_i), cos(phi_i))
```

This is deterministic and reproducible: `(shell_index, shell_slot_index)`
is a stable address that always maps to the same galaxy position, which is
useful for `galaxyGen.py` later (a shell can be regenerated or extended
without perturbing sectors already placed in other shells) — but that
determinism is itself a design choice worth flagging; see §7.

### Why this is fine despite the gaps/overlaps

- Sectors are **sparse content cells**, not solid masonry. This project's
  own stellar density (`LOCAL_STELLAR_DENSITY_LY3`, ~1 system per 352
  cubic light-years) means a default-size sector (~1,521 ly³) holds only a
  handful of systems on average (`expected_system_count`) — most of a
  sector's cube is already empty space. A seam of a few percent of an
  edge-length between two adjacent cubes has no content or gameplay
  consequence; it isn't a wall or a floor that needs to line up.
- "Nearest sector" / adjacency should be defined by **center-to-center
  distance**, exactly the pattern this codebase already uses one level
  down (`SpaceSector.nearest_neighbors`/`distance_between` for systems
  within a sector). It never needed shared-face cube topology at the
  system level, and doesn't need it at the sector level either.
- The alternative — a sector shape that *does* tile a sphere exactly — has
  no good answer for "cubic," since no polyhedron tiles a sphere's surface
  the way a cube tiles a plane. Any exact solution would mean giving up
  the fixed cubic sector shape entirely, which §0 already puts out of
  scope.

### The disk-vs-sphere gap (real galaxies aren't spherical)

Shells as defined above populate a full sphere uniformly in `phi.` A real
disk/spiral galaxy is flattened — the Milky Way's stellar disk has a scale
height on the order of ~300 pc (thin disk) to ~900-1,000 pc (thick disk),
tiny next to its ~15,000 pc radius. Treating the whole sphere as equally
"real" would put as many addressable sectors far off the galactic plane
(large `|z|`, i.e. `phi` near 0 or `pi`) as near it, which is not how the
galaxy this package is modeling actually looks.

This design deliberately **separates addressing from population**:

- The shell/Fibonacci-sphere scheme above defines the full set of
  addresses a sector *could* occupy — this is the coordinate system, and
  it's what this document is scoped to.
- Whether `galaxyGen.py` actually generates content for a given
  `(shell_index, shell_slot_index)` address is a separate, later decision
  — e.g. gating generation on `|z| / r` against a disk-density envelope,
  or a radial density falloff for the bulge/disk/halo. That's a galactic
  population/astrophysics model, not a coordinate system, and per this
  doc's scope (§0) it's flagged as an open question (§7) rather than
  decided here.

For scale, restricting to a cylinder of the galaxy's full 15,000 pc radius
but only ~1,000 pc total thickness (a generous stand-in for "thin + thick
disk") cuts the addressable volume by about 20x versus the full sphere
(`pi * 15,000^2 * 1,000` vs. `(4/3) * pi * 15,000^3`) — from ~320 billion
down to ~16 billion sectors. Still enormous, but it's the first concrete
argument that a disk envelope isn't just astronomically motivated, it's
also a large, free reduction in how much of the address space ever needs
sectors actually generated.

### Cube orientation

Each sector's cube needs an orientation, not just a center point, since
`edge_mpc` alone doesn't say which way the cube's local axes point in the
galaxy frame. The natural default needs no extra stored state at all: align
the cube's local "up" (its own `+Z`, i.e. the direction its own octant
scheme's `phi`-sign axis points) with the **radial direction** from the
galactic center to the sector's own center (`normalize(x, y, z)`) — like a
brick pointing outward, consistent for every sector, and fully derivable
from the already-stored center point. The remaining one degree of freedom
(rotation *around* that radial axis, i.e. "roll") needs a convention too;
the simplest is to derive it from a fixed reference (e.g. project galactic
`+Z` onto the plane perpendicular to the radial direction, and align the
cube's local `+X` to that projection — breaks down only exactly on the
galactic axis where a sector's radial direction *is* `+Z`, an edge case
with no sectors in practice at galaxy-core-adjacent shells).

Because this whole scheme is derivable from the center point with a fixed
convention, **no orientation column is proposed** in §4. If per-sector
random roll (rather than a fixed convention) turns out to be wanted later,
that would need a stored angle — flagged as an explicit open question in
§7 rather than added preemptively.

## 4. Schema migration proposal (v3 -> v4)

New nullable columns on `sectors`, mirroring exactly how
`star_systems.position_x/y/z_mpc`/`quadrant` are nullable ("never placed"
is a real, valid state — every sector generated by today's
single-sector-only tooling has no galaxy position at all):

```sql
ALTER TABLE sectors ADD COLUMN center_x_pc REAL;
ALTER TABLE sectors ADD COLUMN center_y_pc REAL;
ALTER TABLE sectors ADD COLUMN center_z_pc REAL;

-- Derived from center_x/y/z_pc (sqrt(x^2+y^2+z^2)) but persisted anyway,
-- same "quadrant" precedent -- lets "sectors within radius R of the core"
-- use a plain indexed range scan instead of a non-indexable expression.
ALTER TABLE sectors ADD COLUMN galactic_radius_pc REAL;

-- Which radial shell (see Sec 3) this sector's center falls in, and its
-- placement index within that shell's Fibonacci-sphere ordering. Together
-- these are a stable, reproducible "address" a shell can be regenerated
-- or extended from without perturbing sectors already placed elsewhere.
ALTER TABLE sectors ADD COLUMN shell_index INTEGER;
ALTER TABLE sectors ADD COLUMN shell_slot_index INTEGER;

CREATE INDEX IF NOT EXISTS idx_sectors_galactic_radius_pc ON sectors(galactic_radius_pc);
CREATE INDEX IF NOT EXISTS idx_sectors_shell_index ON sectors(shell_index);
```

Because SQLite's `ALTER TABLE ADD COLUMN` can't add a multi-column `CHECK`
after the fact, the "all four galaxy-position columns are NULL together or
all set together" invariant (matching the existing
`quadrant`/`position_x_mpc` null-together convention on `star_systems`)
would need to move into `schema.sql`'s `CREATE TABLE sectors` body directly
the next time this table's DDL is rewritten wholesale, e.g.:

```sql
CHECK (
    (center_x_pc IS NULL) = (center_y_pc IS NULL) AND
    (center_y_pc IS NULL) = (center_z_pc IS NULL) AND
    (center_z_pc IS NULL) = (galactic_radius_pc IS NULL)
)
```

`shell_index`/`shell_slot_index` are left out of that CHECK deliberately —
they're bookkeeping for the placement *algorithm*, not physically implied
by a center point the way `galactic_radius_pc` is, so a sector could in
principle have a galaxy position without having been placed by this
particular shell scheme (e.g. a hand-authored one-off position).

`PRAGMA user_version` moves from `3` to `4`. Per `docs/database-schema.md`'s existing
history format:

> **v3 -> v4**: added `sectors.center_x/y/z_pc`, `galactic_radius_pc`,
> `shell_index`, `shell_slot_index` — a sector's position in a galaxy-scale,
> spherical-coordinates-from-center layout (stored as Cartesian; see
> `docs/design/galaxy-coordinate-system.md`). NULL for a sector never
> placed in a galaxy (every sector generated by today's `sectorGen.py`
> single-sector tooling, and any migrated from an older database).
> `migrate_database`'s `_migrate_v3_to_v4` would need to be added, following
> `_migrate_v2_to_v3`'s exact shape (attach the old file read-only, copy
> every table straight across except `sectors`, which needs the explicit
> pre-v4 column list instead of `INSERT ... SELECT *` since the current
> schema has more columns than a v3 database's `sectors` table does), and
> per `migrate_database`'s own docstring, **`_migrate_v1_to_v2` and
> `_migrate_v2_to_v3` would both need updating too** to leave the new
> columns NULL, since each of those functions maps straight to the
> *current* schema in one hop rather than chaining through intermediate
> versions.

No backfill step is possible for `_migrate_v3_to_v4` beyond leaving the new
columns NULL — there is no way to recover a legacy sector's intended galaxy
position after the fact, exactly the same situation `location`/`quadrant`
are in for a `star_systems` row that was never placed in a sector at all.

## 5. Interaction with existing physics — flagged, not resolved

`stellarObjects/physical_constants.py` already has:

```python
GALACTIC_CENTER_DISTANCE_LY = 25800  # Distance from Sol to the Galactic Center in light-years
```

`Star.calculate_system_perimeter` (`starData.py:341-361`) and
`BinaryStarProxy`'s equivalent (`doubleStar.py:270-282`) both use this
fixed constant for **every** generated system, regardless of which sector
it's in — i.e., today, every star in the game is implicitly assumed to sit
at Sol's own galactocentric distance (~25,800 ly, ~7.9 kpc) when computing
its Hill sphere relative to the galaxy (`system_perimeter`). That value
then drives `hill_radius_ly`/`required_separation_ly` in `spaceSector.py`
— the minimum-separation physics used to place systems within a sector.

Once sectors have a real galactic radius (`galactic_radius_pc`), this
constant becomes **wrong in a specific, physically meaningful way** for
any sector not near 7.9 kpc out: a sector near the galactic core should
have a *smaller* Hill sphere ceiling per star (stronger galactic tidal
shear, tighter Oort-cloud limits, systems packed closer together), and a
sector out in a sparse halo shell should have a *larger* one. This is a
real interaction this design surfaces but does not resolve — changing
`calculate_system_perimeter` to take a distance parameter (falling back to
`GALACTIC_CENTER_DISTANCE_LY` when no sector position is known, for
standalone/non-galaxy-placed generation) is generation-code, out of this
document's scope, and is called out explicitly in §7 as a judgment call:
is this worth doing now, alongside the coordinate system, or deferred as
its own follow-up once `galaxyGen.py` exists and it's clear how often
sectors actually land far from 7.9 kpc?

## 6. Worked end-to-end example

A sector placed at shell `k = 50`:

```
edge_pc         = 3.526                          (DEFAULT_SECTOR_EDGE_LY)
r_50            = (50 + 0.5) * 3.526             = 178.06 pc
N_50            = round(4 * pi * 50.5^2)         = 32,047 sectors in this shell
```

Picking slot `i = 12,000` of those 32,047 via the golden-angle sequence
(illustrative, not solved to full precision here):

```
phi_i   = acos(1 - 2 * 12000.5 / 32015)   ~= 1.318 rad  (~75.5 deg from galactic north)
theta_i = (2*pi*12000 / 1.618...) mod 2*pi ~= 2.564 rad (~146.9 deg)

center_x_pc ~= 178.06 * sin(1.318) * cos(2.564) ~= -144.4
center_y_pc ~= 178.06 * sin(1.318) * sin(2.564) ~=   94.2
center_z_pc ~=  178.06 * cos(1.318)             ~=   44.6

galactic_radius_pc = sqrt((-144.4)^2 + 94.2^2 + 44.6^2) ~= 178.1   (matches r_50, as it must)
```

In light-years for a display string (`pc_to_ly`): `178.06 pc * 3.2616 ~=
580.8 ly` from the galactic core — small enough that
`GALACTIC_CENTER_DISTANCE_LY`'s fixed 25,800 ly assumption (§5) would be
off by a factor of ~44x for a system generated in this sector, which is
exactly the kind of case that flag matters for.

## 7. Open questions — for the user's judgment, not decided here

1. **Eager vs. lazy shell/sector generation.** §3 shows shell counts grow
   quadratically (~227 million sectors in the single outermost shell of a
   15,000 pc galaxy) and total addressable volume runs to ~320 billion
   sectors (or ~16 billion restricted to a disk-shaped envelope). `galaxyGen.py`
   almost certainly cannot generate "the whole galaxy" eagerly — but
   whether it should generate whole shells on demand, individual sectors
   on demand (e.g. as a player/query "visits" a region), or some other
   unit of batching is a `galaxyGen.py` design question this document
   deliberately leaves open (per docs/TODO.md, that script's own design is
   gated on this one, not the reverse).
2. **Disk-density envelope.** §3 flags that a pure uniform sphere doesn't
   match a real disk galaxy's shape, and sketches the *concept* of gating
   generation by `|z|`/scale height, but does not specify a concrete
   profile (thin/thick disk split, bulge, spiral-arm overdensity,
   metallicity-vs-radius gradients feeding into star-type selection).
   That's a galactic astrophysical population model, arguably its own
   design pass, layered on top of the coordinate system this document
   defines.
3. **`GALACTIC_CENTER_DISTANCE_LY` becoming position-dependent** (§5) —
   whether `calculate_system_perimeter` should start taking the owning
   sector's actual `galactic_radius_pc` (converted to ly) instead of the
   fixed Sol-distance constant, and if so, whether that lands alongside
   this migration or as a separate follow-up once real galaxy-scale data
   exists to see how much it matters in practice.
4. **Cube roll/orientation.** §3 proposes a fixed, zero-extra-storage
   convention (radial-outward `+Z`, projected-galactic-north `+X`). If a
   per-sector random roll is wanted instead (e.g. for visual/narrative
   variety), that needs a stored orientation value (a single roll angle
   suffices, given the radial axis itself is already fixed by the center
   point) — not proposed here since it's pure additional complexity with
   no identified requirement yet. §9's `sector_vertices` table implements
   the fixed convention concretely (as the tangent-plane basis its exact
   prism vertices are derived in), but still doesn't add a roll degree of
   freedom -- that half of this question remains open.
5. **Deterministic Fibonacci-sphere placement vs. randomized placement.**
   §3's scheme is fully deterministic — the same `(shell_index,
   shell_slot_index)` always yields the same position, which makes shells
   independently regenerable/extensible, but also means sector placement
   has none of the organic irregularity `spaceSector.py`'s own
   system-within-a-sector placement deliberately has (Poisson-disk growth,
   `secrets`-backed randomness). Whether galaxy-scale sector placement
   should also carry some jitter (e.g. a small random offset within each
   Fibonacci-sphere cell, sized so it can't push a sector's cube into a
   neighboring shell) is an open aesthetic/design call, not a geometric
   necessity.
6. **Non-uniform sector size.** This entire scheme assumes one fixed
   `edge_pc` for every sector (today's actual default, `edge_mpc` is
   already a per-row column so it technically *can* vary — but §3's shell
   math assumes a single edge length shared by a whole shell's sector
   count formula). Larger, cheaper sectors farther from the core, or
   smaller/denser ones near it, would need a materially different
   version of §3's `N_k` formula (shell thickness and sector footprint
   would no longer be shell-wide constants). Not addressed here; assumed
   out of scope unless the user wants it folded in now.
7. **One galaxy per database.** This design assumes a single galactic
   center/origin per database file (there is no `galaxy_id` concept
   anywhere in the schema or in `docs/TODO.md`'s framing, which only ever says
   "the galaxy," singular). If multiple independent generated galaxies
   ever need to coexist in one database, every new column in §4 would need
   a `galaxy_id` companion — flagged in case that's a real future
   requirement rather than a hypothetical one.
8. **Sector-to-sector adjacency as a queryable relationship.** §3's
   "nearest sector by center-to-center distance" answer works fine as a
   one-off computed query (mirroring `SpaceSector.nearest_neighbors`), but
   if `galaxyGen.py` or a future web UI wants persistent "neighboring
   sectors" (e.g. for travel/lore, similar to how `star_systems.location`
   persists each system's nearest neighbors today), that would need its
   own join table and its own write-time computation — not proposed here,
   since nothing in the current brief calls for it yet.

## 8. Generation unit: sector enumeration by radius (Track C addendum)

**Status:** implemented (`src/stellarObjects/galaxyGeometry.py`). This
section is a Track C addition, written after this document's original
proposal (sections 0-7 above, including the open questions in §7 —
several since resolved by the user for this track, see the top of this
document's "Decisions already made" recap and Track C's own brief) and
after the user decided that `galaxyGen.py` needs one primitive supporting
both whole-shell batch generation and on-demand generation of the
neighborhood around an *already-placed* sector — not just "generate a
whole shell" as §7 question 1 originally framed the choice. That
primitive is specified and analyzed here before its implementation, per
this track's own requirement that the geometry sub-problem get a short
design pass of its own.

### The primitive

> Given an arbitrary point `P = (x, y, z)` in galaxy-space (parsecs, not
> necessarily the galactic origin, and not necessarily an already-placed
> sector's center) and a radius `R`, enumerate every `(shell_index,
> shell_slot_index)` address (§3's tiling scheme) whose sector center
> falls within `R` of `P`.

Both `galaxyGen.py` modes reduce to one call of this primitive:

- **Batch mode** (`--shell k`): `P` = the galactic origin, `R` = shell
  `k`'s own outer radius (or simply `shell_radius_pc(k, edge_pc)` plus a
  hair of slack) — returns (up to) that whole shell's `N_k` addresses.
- **Local-neighborhood mode** (`--center-sector <id> --radius-pc R`): `P`
  = that sector's own stored `(center_x_pc, center_y_pc, center_z_pc)`, a
  small `R` — returns just the handful-to-few-dozen addresses genuinely
  near it, regardless of which enclosing shell they fall in (a
  neighborhood can straddle a shell boundary) and without requiring that
  enclosing shell to have been generated as a whole first.

### Why P is not always the origin, and why that's the hard part

Every position in this scheme is derived from `(shell_index,
shell_slot_index)` via `sector_position_pc` (§3's Fibonacci-sphere
formula), which is O(1) per address — trivial in the batch case, since
`P` = origin means "distance to `P`" is just each shell's own fixed radius
`r_k`, so the *whole* shell either qualifies or doesn't (no per-slot
geometry needed at all). The local-neighborhood case is harder precisely
*because* `P` is generally not the origin, not on any particular shell's
"axis" of anything, and not aligned with the Fibonacci-sphere's own
golden-angle indexing in any convenient way. A shell far out (say `k` in
the thousands, `N_k` in the hundreds of millions per the §3 table) must
never be iterated in full just to answer "which of your slots are within
3 sector-widths of this one sector I already have" — that would make the
very use case this primitive exists for (generate a local starmap around
a point of interest) slower than generating the whole enclosing shell it
is trying to avoid.

### Chosen approach: two cheap, exact prunes, then brute-force the rest

No spatial index (k-d tree, octree, or similar) is built over the address
space — the Fibonacci-sphere scheme already has enough structure to prune
analytically, in closed form, without one. Two independent pruning passes
run before any per-slot position is ever computed:

1. **Which shells can possibly qualify** (`_candidate_shell_range`).
   Every slot in shell `k` sits at the *exact same* radius `r_k` (§3's
   placement is a fixed-radius sphere per shell, not a radius range — the
   `[k*edge_pc, (k+1)*edge_pc)` interval in §3 is the shell's conceptual
   thickness, not where individual slots actually land). So the question
   "could shell `k` hold a point within `R` of `P`" has an exact answer,
   not an approximation: yes iff `|r_k - |P|| <= R` (the minimum possible
   distance from `P` to *any* point on a sphere of radius `r_k`, achieved
   when that point lies on the ray through `P`, is exactly `|r_k - |P||`;
   the triangle inequality guarantees no closer point exists). This bounds
   the candidate shell range to `k` in roughly
   `[(|P|-R)/edge_pc - 0.5, (|P|+R)/edge_pc - 0.5]` — a width of about
   `2R/edge_pc` shells, independent of how far out `P` is. For a "handful
   to a few dozen sectors" neighborhood (`R` on the order of a few
   `edge_pc`), that's a small, constant number of candidate shells
   regardless of galactic radius.

2. **Which slots within a candidate shell can possibly qualify**
   (`_slot_index_bounds_for_phi_range`). This is the part that keeps
   per-shell cost from scaling with `N_k`. `sector_position_pc`'s polar
   angle `phi_i = acos(1 - 2*(i+0.5)/N_k)` is *strictly monotonic* in the
   slot index `i` (as `i` runs `0..N_k-1`, the argument to `acos` runs from
   just under `+1` down to just over `-1`, and `acos` is strictly
   decreasing) — so a target range of `phi` values maps to exactly one
   contiguous range of slot indices, invertible in closed form:
   `i = N_k*(1 - cos(phi))/2 - 0.5`. The spherical law of cosines gives the
   maximum angular separation `alpha_max` (from `P`'s own direction, as
   seen from the origin) a point at radius `r_k` can have and still be
   within `R` of `P`:
   `cos(alpha_max) = (|P|^2 + r_k^2 - R^2) / (2*|P|*r_k)`. A short
   trigonometric argument (in the module docstring's
   `enumerate_sectors_within_radius`) shows `|phi_slot - phi_P| <= alpha_max`
   is a *necessary* condition for any slot regardless of its azimuthal
   angle `theta` — so the contiguous index range this maps to is a safe,
   exact upper bound on which slots could possibly qualify, even though it
   isn't tight in azimuth (two slots at the same `phi` but very different
   `theta` are treated alike by this prune; the exact per-slot distance
   check below is what actually discriminates on `theta`).

Every slot surviving both prunes gets its exact position
(`sector_position_pc`) and exact distance to `P` computed and compared to
`R` — so the final output is always exact, never an approximation; the
two prunes only decide *which* slots are worth that O(1) exact check, not
whether the check itself is skipped.

Two edge cases fall out of the same math rather than needing special
handling: `P` at the origin (no direction to prune slots by — but then
distance-to-`P` for any slot is just `r_k`, so a qualifying shell is
either wholly in or wholly out, handled before the per-slot loop even
starts), and `R` large enough that a whole candidate shell is guaranteed
inside it (`cos(alpha_max) <= -1`, i.e. `R >= |P| + r_k`, the shell's own
worst-case distance from `P`) — that shell's slots are all yielded
directly, skipping the angular prune (there's nothing left to prune).

### Big-O

Let `S` = number of candidate shells (`O(R/edge_pc)`, a small constant for
a local-neighborhood-sized `R`, independent of `|P|`), and let `B_k` = the
slot-index band width computed for shell `k` (bounded by, roughly, `N_k`
scaled by the fraction of that shell's surface subtended by the angular
cap of half-angle `alpha_max` — for `R << r_k`, `alpha_max ~ R/r_k`
radians, so `B_k` grows roughly with the actual *area* of the search
region on that shell, not with `N_k` itself). Total work is
`O(sum_k(B_k))`, which scales with the genuine output size (how many
sectors actually exist within `R` of `P` across the candidate shells) plus
a small constant-factor overhead from the band edges — **not** with
`sum_k(N_k)`, the total number of slots those shells hold. This is what
makes a local-neighborhood query near an outer shell (potentially
hundreds of millions of slots) exactly as cheap as one near the core: the
per-shell cost tracks `R`, not `N_k`.

The batch-mode case (`P` = origin, `R` = one shell's radius) is the
degenerate case of the same function: exactly one candidate shell
qualifies (or a small handful, if `R` is chosen slightly loose), the
origin special-case skips the angular prune entirely, and every one of
that shell's `N_k` slots is legitimately part of the answer — so batch
mode's cost is, correctly, `O(N_k)`: there is no way to return `N_k`
results in less than `O(N_k)` time, and no pruning is being asked to do
anything there.

### What was deliberately not built

No k-d tree, octree, or other general-purpose spatial index over the
address space, and no precomputed per-shell angular lookup structure
beyond the closed-form inversion above. Both were considered and rejected
as premature for a "handful to a few dozen sectors" neighborhood query,
per this track's brief: the two closed-form prunes above already remove
the only real scaling hazard (iterating a huge shell's full slot count),
and building a persistent index would add real complexity (index
construction, invalidation as new sectors are generated, storage) for a
problem this document's own §7 open question 8 already flags as
explicitly out of scope for now (persistent sector adjacency). If a future
need arises for large-`R` neighborhood queries (say, `R` spanning dozens
of shells rather than a handful of `edge_pc`), the same two prunes still
apply without modification — they degrade gracefully to "iterate more
shells, and wider slot bands within each" rather than breaking down, so
there is no cliff where this approach stops working; it simply does
progressively more of the genuinely-necessary work as the requested
neighborhood grows.

## 9. Sector prism vertices: exact per-shell Voronoi tessellation (addendum, revision 2)

**Status:** implemented (`stellarObjects/sectorGeometry.py`, the
`sector_vertices` table — schema v7; briefly a `sectors.vertices_pc` JSON
column in schema v6, reconsidered immediately in favor of a normalized
table, since this schema has no JSON-blob columns anywhere else). This
revises an earlier version of this same
addendum, which gave every sector a fixed 8-vertex cube and *nudged*
corners toward nearby neighbors' — that approach shipped, worked, and was
tested, but only ever reduced gaps (~35% aggregate improvement, measured
against real `galaxyGen.py` output), never eliminated them, because a
fixed 8-vertex/6-face shape cannot exactly reconcile a sector with more
real neighbors than it has faces — common on this Fibonacci-sphere
placement, which has no fixed "6 neighbors" the way a structured grid
would. Following a direct request for literal zero gaps, this revision
replaces corner-nudging with an **exact local spherical Voronoi
tessellation**, letting vertex/face count vary per sector instead of
staying fixed.

### The model

- **Lateral (same-shell) sharing is exact, not approximate.** For sector
  `P`, `stellarObjects.sectorGeometry.local_lateral_cell` finds `P`'s real
  same-shell geometric neighbors and computes the cell boundary as the
  exact 3D **circumcenter** of each pair of cyclically-adjacent neighbors
  together with `P` — the one point in 3D equidistant from all three. This
  is a plain geometric fact about three points, independent of which of
  them "does the computing," so two real neighbors, computing their own
  cells entirely independently, land on the identical floating-point value
  for their shared corner (verified directly: ~1e-16 agreement, i.e.
  floating-point noise, not an approximation residual). Vertex count
  varies per sector (typically 5-7, mean 6.0, matching standard Voronoi/
  Euler-formula theory for a near-uniform point set) — the direct,
  necessary consequence of insisting on exact gaps: a cube tiling of a
  sphere cannot be gap-free in general (the same reason a soccer ball
  needs pentagons mixed with hexagons; a plane tiles perfectly with
  squares, a sphere never does), so face count has to match each sector's
  own real neighbor count instead of staying fixed.
- **Radial (between-shell) coverage matches in area, not vertex-for-
  vertex.** Each lateral vertex is scaled along its own ray from the
  galactic origin to sit exactly on the shell's inner bound
  (`shell_index * edge_pc`) and outer bound (`(shell_index + 1) *
  edge_pc`) in turn (`prism_vertices`), giving every sector a radially-
  extruded prism. Shell `k`'s outer bound and shell `k+1`'s inner bound
  are the same sphere; each shell tiles that whole sphere completely and
  independently via its own sectors (Voronoi cells always partition their
  full surface), so there is no net gap in area between the two shells'
  sectors even though the two tessellations don't share edges with each
  other — a "non-conforming mesh interface," the same technique used
  where two independently-meshed regions meet in finite-element/CFD
  meshing. Shell 0 is a degenerate but correct special case: its inner
  bound is radius 0, so its sectors are wedges/cones from the galactic
  center rather than full prisms.
- No roll/orientation degree of freedom was added — §3's fixed convention
  (radial-outward local `+Z`, projected-galactic-north local `+X`) is
  still used as the tangent-plane basis the lateral cell's connectivity is
  worked out in, even though final vertex positions are exact 3D
  circumcenters rather than tangent-plane approximations.

### The performance problem this required solving, and how

A naive same-shell neighbor search — reusing
`galaxyGeometry.enumerate_sectors_within_radius` with a small physical
radius — is fast near a shell's poles but **catastrophically slow near its
equator for large outer shells**: that primitive prunes by polar angle
(`phi`) alone, and this placement's slot index is uniform in `cos(phi)`,
not `phi` itself, so a tiny physical radius maps to a huge slot-index range
right at the equator (measured: 100,000+ candidates to find ~6-8 true
neighbors, ~0.4-0.5 seconds per sector) — precisely where the disk/spiral
density model (this session's other major addendum) concentrates real
generation activity.

The fix exploits what this placement actually *is*: a Fibonacci sphere
built from the golden angle, which makes it a golden-ratio irrational
rotation in disguise. By the three-distance theorem, two slot indices land
at close azimuths essentially when their difference is a Fibonacci number
(golden-ratio continued-fraction convergents are exactly consecutive
Fibonacci numbers) — confirmed directly against real generated positions,
where true neighbor offsets were exactly `F_20` through `F_23`
(6765/10946/17711/28657) plus small integer combinations of adjacent pairs
(e.g. `76 = 2*F_10 - F_9`). Combined with polar angle changing
(approximately) linearly with index, the relevant offset scale works out
to `sin(phi) * sqrt(pi * N)` — shrinking away from the equator.
`sectorGeometry._same_shell_candidate_offsets` checks small integer
combinations of the few Fibonacci-number pairs nearest that scale, cutting
per-sector cost to ~0.3-0.4ms regardless of shell size (a ~1000x
improvement for the worst outer-equatorial case) — validated against the
guaranteed-correct brute-force search across 840 cases spanning the full
polar range and shell sizes from 3 to 211 million slots, with **zero
mismatches** once shells small enough that this asymptotic theory doesn't
apply cleanly (`shell_sector_count(k) <= 2000`, in practice only the
galactic-core-adjacent shells) fall back to plain brute force instead,
which is unconditionally correct and still cheap at that size.

### Two correctness bugs found by an end-to-end simulation, and fixed

Running a real (small, unfilled) galaxy through this module end-to-end —
sector layout, density, and vertices for every qualifying sector, with an
exhaustive check that every outer vertex is shared with another same-shell
sector's cell — surfaced two bugs neither unit test caught, both specific
to small/sparse shells (the galactic-core-adjacent shells that unit tests,
by convention, under-sampled in favor of large ones like shell 50):

1. **Naive tangent-plane projection understated distance.** The same-shell
   candidate projection originally used the raw chord vector's own
   components along the local axes (`dot(candidate - position, local_x)`
   etc.), which systematically *understates* how far away a candidate
   really is, worse the farther it is (a real 10.5 pc separation projected
   to as little as 1.5 pc on shell 1). Fixed with a proper **gnomonic**
   (central) projection (`_gnomonic_projection`): following the ray from
   the galactic origin through the candidate out to where it crosses the
   query sector's tangent plane, which grows monotonically with true
   angular separation and is the standard technique for reducing a
   spherical Voronoi problem to a planar one.
2. **The half-plane test after that projection used the wrong bisector.**
   Gnomonic projection maps the true spherical bisector between two
   co-radial points to a straight line in the projected `(u, v)` plane —
   but not to the *flat*-plane perpendicular bisector of `(0, 0)` and
   `(u, v)` (`u*x + v*y <= (u**2 + v**2) / 2`), which was the formula used.
   The two agree only in the small-angle limit (why large/dense shells
   were unaffected); on a small, sparse, fully-populated shell it was
   permissive enough to let two non-adjacent sectors' cells meet at a
   point a third, genuinely closer sector should have cut off first — a
   real, silent gap. The correct right-hand side, derived from equidistance
   in 3D between the query sector (radius `r_k`) and a co-radial candidate
   projected to `(u, v)`, is `r_k * (sqrt(r_k**2 + u**2 + v**2) - r_k)`.

Both fixes are covered by regression tests parametrized to include a small
shell (shell 1) alongside the large ones already under test, plus an
exhaustive "every vertex of every sector in a fully-populated small shell
is shared with another sector" check
(`test_local_lateral_cell_fully_tiles_a_small_shell_with_no_orphan_vertices`)
that reproduces the exact condition the simulation used to find these bugs.
Re-running the same simulation after both fixes: **0 unexplained gaps**
across all 31,255 outer vertices generated (previously 4,798, then 164 as
each fix landed).
