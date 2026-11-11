# The gas disk's outer edge (GEN.1520)

Status: built (Lane 14, 11-11). Decision: Boss, 3:04 PM EST 11-11. Evidence: `reports/gen1520-disk-size-by-mass-2026-11-11.md`.

## What it is

The cut-off radius of the gas disk a star's planets, belts and comets formed in. Before GEN.1520 it was 18 times the star's own snow line (49 au for the Sun); it now follows the star's mass at formation:

    R_edge = 100 au x (M / Msun)^0.5, times 10^(0.3 z), z a normal draw cut at 2 sigma

so the median is 45 au at 0.2 Msun, 71 au at 0.5 Msun, 100 au at 1 Msun and 118 au at 1.4 Msun, and any one star lies within a factor of 4 of that.

## Why

Observed dust disks have a median radius of 27 to 60 au (Hendler 2020), the CO gas is a median 2.5 times bigger (Sanchis 2021) and the planet-formation models start at 62 au for a Sun-like star (Burn 2021). The snow line scales with the star's light, which is the wrong variable: direct fits of disk size against star mass give M^0.2 to 0.6 and none against luminosity strong enough to matter. The constants are stand-ins DS1 to DS4 (`gap-filling-plan.md` section 6.274).

## Where it is read (one place, no stored column)

`physics/disk_size.edge_au(star)` draws the edge from a stream keyed by the star's mass, age and place, the way the snow line is drawn, so every caller gets the same value and no other draw moves. A white dwarf or other remnant uses its progenitor's mass (`mass_loss.member_formation_mass_kg`); a close pair uses the sum of its members.

- `generation/system.py`: where random planet slots end and the oligarch chain's outer limit (`_disk_outer_edge_au`, `_plan_architecture`).
- `generation/belt_birth.py`: the outer-disk belt and the scatter shares.
- `generation/comet.py`: the outer edge of the zone a comet's ice formed in.
- `physics/disk_drift.py` through `mixture_subject` and `mixture_giant`: where the pebbles end and which ices evaporate on the way in (GEN.1473). A larger edge is a colder one (the disk's temperature falls as the square root of the distance), so more of the volatile ices (CO2, CO, N2) have their lines inside the disk.

## What it changes

Planet, giant, belt and comet counts, the outer reaches of the planet lists, the comet reservoirs' ice and the C/O of giants. It is a lock item; the sweep before and after is in the report.

## Not done

A scatter that depends on the disk's mass, a disk age spread, and a different edge for the dust (the gas edge is what the model uses; the dust edge would be about 1/2.5 of it).
