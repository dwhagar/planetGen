# Galaxy Disk/Spiral Density — Design Proposal

**Status:** proposal, not implemented. This is the design pass
`docs/design/galaxy-coordinate-system.md` section 7, question 2 ("Disk-density
envelope") calls for, and the remaining open item under `docs/TODO.md`'s
Phase 4 ("Galaxy thickness / shape ... still needs its own dedicated design
pass"). Track C (that document) defined *where sector addresses can exist*
(the shell/Fibonacci-sphere tiling) and deliberately left *which addresses
actually get content, and how richly* undecided — that's what this document
decides.

**Scope:** a density model over galaxy-frame position, and how `galaxyGen.py`
uses it to decide (a) whether a given `(shell_index, shell_slot_index)`
address gets a sector generated at all, and (b) how many systems that sector
gets. Explicitly **out of scope**: metallicity-vs-galactic-radius gradients
feeding star-type/composition selection (flagged in the coordinate doc's §7
Q2 as related but separate — that would change *what* gets generated inside
a system, not *whether/how much* gets generated at a galaxy position) and any
schema change (nothing here needs persisting — see §5).

## 0. Decisions carried over from the user

- **Gate + weight, not weight-only.** Low-density addresses (deep halo,
  between arms at large radius) are probabilistically skipped and never
  generated at all, not merely generated sparse. This is also what makes
  the coordinate doc's "~20x volume reduction from a disk envelope" argument
  real rather than theoretical (§4).
- **Default behavior**, not opt-in. `galaxyGen.py --shell`/`--center-sector`
  use this model by default; a `--uniform` flag reverts to today's flat,
  isotropic-sphere behavior (useful for tests and for comparing against the
  old behavior). Already-generated sectors are untouched either way — this
  only affects new generation going forward.
- **No hard galaxy edge.** No new "galaxy radius" constant gates generation
  outright. Density falls off exponentially with radius, so the occupancy
  gate (§3) already makes very distant shells generate next to nothing on
  their own, without an arbitrary cutoff radius/shell index to pick and
  maintain.

## 1. Model

Standard exponential-disk-plus-bulge-plus-spiral galaxy model, evaluated as
a pure function of galaxy-frame Cartesian position — no new stored state,
just arithmetic over the `(x, y, z)` a sector already has (or would have,
for a not-yet-generated candidate slot, via
`galaxyGeometry.sector_position_pc`).

### Cylindrical decomposition

Disk galaxies aren't naturally described in the coordinate doc's own
spherical `(r, theta, phi)` (defined from the galactic center in every
direction alike) — a disk needs a radius *within the plane* and a height
*off* the plane treated asymmetrically. New helper, local to this model
(kept out of `galaxyGeometry.py` on purpose — that module is the addressing
scheme, this is the population model layered on top, per the coordinate
doc's own "deliberately separates addressing from population" framing):

```
r_cyl = sqrt(x^2 + y^2)      # in-plane galactocentric radius
theta = atan2(y, x)          # in-plane azimuth (same convention as galaxyGeometry)
z     = z                    # height off the galactic plane
```

### Components

```
rho_bulge(r_3d)         = bulge_amplitude * exp(-r_3d / bulge_scale_radius_pc)
rho_disk_radial(r_cyl)  = exp(-r_cyl / disk_scale_length_pc)
f_z(z)                  = sech(z / disk_scale_height_pc)^2
arm_factor(r_cyl, theta) = 1 + arm_amplitude * cos(
                               arm_count * (theta - theta_arm(r_cyl))
                           )
theta_arm(r_cyl)        = spiral_reference_angle_rad
                           + ln(r_cyl / spiral_reference_radius_pc) / tan(pitch_angle_rad)

relative_density(x, y, z) =
    rho_bulge(r_3d)
    + rho_disk_radial(r_cyl) * f_z(z) * arm_factor(r_cyl, theta)
```

where `r_3d = galaxyGeometry.galactic_radius_pc((x, y, z))` (the existing
spherical-radius helper — the bulge is spherical, not flattened, so it uses
true 3D distance, not `r_cyl`).

- `sech(u) = 1 / cosh(u)`, so `f_z` is exactly 1 in-plane (`z = 0`) and
  decays smoothly and symmetrically above/below it — the standard vertical
  profile for an isothermal disk, giving a rounder falloff than a bare
  `exp(-|z|/h)`, which has an (unphysical) sharp point at `z = 0`.
- `arm_factor` is a logarithmic-spiral overdensity: at fixed `r_cyl`, it
  peaks (`1 + arm_amplitude`) exactly on the arm and troughs
  (`1 - arm_amplitude`) exactly between arms, `arm_amplitude` in `[0, 1)`
  keeping it always non-negative. It only modulates the disk term — real
  spiral structure is a disk-population phenomenon; it doesn't touch the
  bulge.
- **Normalization**: parameters are chosen so `relative_density` equals
  `1.0` at `(r_cyl, z, theta) = (R_sun_pc, 0, theta_arm(R_sun_pc) + pi/2)`
  — i.e. at Sol's own galactocentric radius (reusing
  `physical_constants.GALACTIC_CENTER_DISTANCE_LY`, converted via
  `utils.ly_to_pc`, rather than inventing a second "reference radius"
  constant), at the *inter-arm average* point (`cos(...) = 0`). This
  deliberately ties the new model to the exact same "1.0 = realistic local
  density" convention `sectorGen.py --density` already uses (§4) instead of
  introducing a second, unrelated normalization.

### Proposed defaults

Not rigorously derived — real-galaxy-scale numbers (Milky-Way-like, since
`GALACTIC_CENTER_DISTANCE_LY` already borrows Sol's real distance),
intended as a first pass to be tuned by eye once implemented (generate a
few example galaxies, look at the resulting sector distribution) rather
than fixed here:

| Constant | Proposed default | Basis |
|---|---|---|
| `disk_scale_length_pc` | ~2,800 pc | `R_sun_pc / disk_scale_length_pc ~= 2.8`, matching the real Milky Way's ratio |
| `disk_scale_height_pc` | ~350 pc | Between real thin-disk (~300 pc) and thick-disk (~900-1,000 pc) values the coordinate doc already cites; a single effective scale height, not a thin+thick blend (see §6 open questions) |
| `bulge_scale_radius_pc` | ~500 pc | Small next to `disk_scale_length_pc`, so the bulge only dominates very close to the core |
| `bulge_amplitude` | ~5.0 | Bulge center reads markedly denser than the disk's own peak (`rho_disk_radial(0) = 1`) |
| `arm_count` | 2 | Grand-design two-arm spiral — the most recognizably "spiral galaxy" shape |
| `pitch_angle_deg` | ~15 deg | Typical grand-design spiral pitch (real spirals run ~10-25 deg) |
| `arm_amplitude` | 0.4 | Arm/inter-arm contrast ratio `(1+0.4)/(1-0.4) = 2.33x` — visible structure, not a binary on/off |
| `spiral_reference_radius_pc` | = `disk_scale_length_pc` | Arbitrary but fixed anchor for `theta_arm`'s log-spiral formula |
| `spiral_reference_angle_rad` | 0.0 | Arbitrary fixed orientation — the galaxy's arms have to start pointing *somewhere* in the coordinate frame; no reason to prefer one angle over another |

All of the above become new constants in `physical_constants.py`, alongside
`GALACTIC_CENTER_DISTANCE_LY` (astrophysically-motivated numbers, not
generation-tuning knobs like `program_constants.DEFAULT_SECTOR_EDGE_LY`).

## 2. Where this lives in code

New module `src/stellarObjects/galaxyDensity.py`, mirroring
`galaxyGeometry.py`'s own shape: pure functions, no I/O, no randomness,
fully unit-testable against hand-computable values.

```
relative_stellar_density(position_pc) -> float   # >= 0, per §1
```

`position_pc` is `(x, y, z)` in parsecs — the same tuple shape
`galaxyGeometry.sector_position_pc`/`enumerate_sectors_within_radius`
already produce, so `galaxyGen.py` can pass either straight through with no
conversion.

## 3. Occupancy gate

For a candidate address whose position is `p`:

```
occupancy_probability = min(1.0, relative_stellar_density(p) ** density_gate_exponent)
```

`density_gate_exponent` (default `1.0`, a new `program_constants.py` knob
since it's a generation-shaping dial, not an astrophysical quantity) lets
the contrast between "generate" and "skip" be sharpened (`> 1`, starker
voids between arms) or softened (`< 1`) later without touching the density
model itself. Any `relative_stellar_density >= 1` (bulge, on-arm regions at
or inside the solar radius) always generates — the gate only ever thins out
sub-realistic-density regions, never regions denser than "typical."

**Deterministic, not RNG-driven.** The decision is a stable hash of the
address itself, not a draw from `galaxyGen.py`'s per-run
`random.seed(secrets.randbits(128))` state:

```
def _occupancy_roll(shell_index, shell_slot_index):
    digest = hashlib.sha256(f"{shell_index}:{shell_slot_index}".encode()).digest()
    return int.from_bytes(digest[:8], "big") / 2**64

occupied_by_gate = _occupancy_roll(shell_index, shell_slot_index) < occupancy_probability
```

This matters because a slot's fate must not depend on *when* or *how many
times* the command generating it happens to run — matching
`sector_position_pc`'s own "same address, same result, forever" guarantee
from the coordinate doc, and avoiding a surprise where re-running
`--shell K --limit 100` twice fills in a different-looking patch of the
shell each time.

### The efficiency point: pruning, not filtering

Evaluating the gate for every one of a shell's `N_k` slots is still `O(N_k)`
even though most distant/off-plane slots will fail it — for an outer shell
(hundreds of millions of slots, per the coordinate doc's own table), that's
still a lot of cheap-but-nonzero work, and it means the "disk cuts
addressable volume ~20x" argument from the coordinate doc's §3 never
actually saves anything at runtime.

Fix: reuse `galaxyGeometry._slot_index_bounds_for_phi_range` (already used
by `enumerate_sectors_within_radius` for the exact same kind of pruning) to
skip most of a shell's slot range outright. Given a shell's fixed radius
`r_k`, a chosen `z_cutoff` (a few `disk_scale_height_pc`, beyond which
`f_z(z)` is negligible enough to treat as exactly zero for generation
purposes — this is the one place the model is deliberately truncated rather
than left as an infinite tail), the polar-angle band that could possibly be
within `z_cutoff` of the plane at radius `r_k` is:

```
phi_disk_min = acos(min(1, z_cutoff / r_k))
phi_disk_max = pi - phi_disk_min
```

`run_shell_batch` computes `(i_min, i_max) =
_slot_index_bounds_for_phi_range(phi_disk_min, phi_disk_max, n_k)` once per
shell and only iterates that band (plus, separately, the handful of slots
near the poles that could still fall inside the *bulge*'s spherical
`r_3d`-based cutoff at small `r_k`, since the bulge isn't flattened) rather
than every slot `0..n_k-1`. For a shell with `r_k` many scale-heights out,
this band is a thin sliver of `n_k`, restoring the real Big-O win a disk
model is supposed to provide, exactly the way the coordinate doc's own
neighborhood-query pruning avoids ever touching a whole outer shell's slot
count.

## 4. Density weighting for generated sectors

`sectorGen.generate_sector` already isolates the exact hook point (line
374-376 today):

```python
if args.density is not None:
    args = copy.copy(args)
    args.num_systems = _sample_poisson_count(sector.expected_system_count() * args.density)
```

Add one new optional parameter, `density_multiplier=1.0`, used only in this
branch:

```python
def generate_sector(args, galactic_center_dist_ly=None, density_multiplier=1.0):
    ...
    if args.density is not None:
        args = copy.copy(args)
        args.num_systems = _sample_poisson_count(
            sector.expected_system_count() * args.density * density_multiplier
        )
```

`galaxyGen.py`'s `_generate_and_save_sector_at` passes
`density_multiplier=relative_stellar_density(position_pc)` (the *uncapped*
value — unlike the occupancy gate, richness isn't clamped at 1: a sector
right in the bulge or on a spiral arm crest should generate meaningfully
more systems than a "typical" one, not just be more likely to exist at
all).

This deliberately only touches the `--density` path, never
`--num-systems`: an explicit fixed system count is the caller overriding
realism outright, and galaxy position shouldn't second-guess that.
`sectorGen.py`'s own standalone CLI is entirely unaffected (default
`density_multiplier=1.0`, identical to today).

## 5. Schema / persistence

**No schema change.** `relative_stellar_density` is a pure function of a
position every generated sector already stores
(`center_x_pc`/`center_y_pc`/`center_z_pc`); anyone who wants "how dense was
it here" later can recompute it exactly, the same way `galactic_radius_pc`
being derivable didn't stop it from also being persisted for query
convenience (coordinate doc §1) — but there's no equivalent query need
here ("sectors above density X" isn't an anticipated query the way "sectors
within radius R" was), so nothing new is added to `sectors`.

## 6. Open questions carried forward (not decided here)

1. **Thin/thick disk split.** §1 uses one effective `disk_scale_height_pc`
   rather than a two-population blend (thin disk, most of the mass, small
   `h_z`; thick disk, a minority, large `h_z`). A blend would look more
   realistic at the vertical margins (a small population of stars well off
   the plane) at the cost of two more constants and a weighted-sum term.
   Deferred as a refinement, not required for a first working version.
2. **Metallicity/star-type gradients** (coordinate doc §7 Q2's other half)
   — whether star type/composition selection should itself vary with
   `r_cyl` (real disk galaxies are more metal-rich toward the core). Stays
   explicitly out of scope here; this document only decides *whether/how
   much* gets generated, not *what*.
3. **Tuning the proposed constants** (§1's table) is a "generate and look at
   it" exercise, not something to over-derive on paper — expect the first
   implementation pass to include a small visualization/inspection script
   (or reuse of the existing Sector Map work) to eyeball arm contrast and
   disk thickness before locking in defaults.

## 7. Testing plan

- `src/tests/test_galaxy_density.py` (new, mirroring
  `test_galaxy_geometry.py`'s style): `relative_stellar_density` is
  hand-checked at a handful of exact points (galactic center, the solar
  calibration point evaluating to `1.0`, on-arm vs. inter-arm at fixed
  `r_cyl`, far off-plane vs. in-plane at fixed `r_cyl`), plus monotonicity
  checks (density strictly decreases with `r_cyl` at fixed `theta`/`z` past
  the arm modulation's own period, and with `|z|` at fixed `r_cyl`).
- `_occupancy_roll` determinism: same `(shell_index, shell_slot_index)`
  always yields the same float across repeated calls/processes.
- `run_shell_batch` with `--uniform`: byte-identical sector count/positions
  to today's behavior (regression guard that the new default path didn't
  change the escape hatch).
- `run_shell_batch` without `--uniform`, on a shell straddling the disk
  plane: statistically, generated slots concentrate at low `|z|` and near
  arm phases; a distant/highly-inclined shell generates few-to-none.
