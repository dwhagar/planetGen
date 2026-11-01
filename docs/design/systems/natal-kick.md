# Natal kick: how many neutron stars a galaxy ejects (GEN.554, slice 3)

Status: built (GEN.554 slice 3). The placement of the kicked stars (the thicker layer, the black hole core excess) is still `remnant_distribution.py`'s empirical factors; deriving it is a later slice. Code: `galaxy/natal_kick.py`, `galaxy/structure_draw.py` (`draw_neutron_star_unbound_share`), `galaxy/structure.py` (the reference structure's share), `galaxy/ejecta_field.py`.

## What changed

Before: each galaxy drew its unbound share of the neutron stars uniformly from 20% to 40% (`NEUTRON_STAR_UNBOUND_FRACTION_RANGE`, "[R1966]'s kick distribution against the escape speed", recalled), whatever its size.

After: the share is what the kick law sends past the galaxy's escape speed. A neutron star is born with the speed of its progenitor, the galaxy's circular speed `Vc` along the disc, plus a kick in a random direction. The kick is a Maxwellian of one-dimensional dispersion 265 km/s (`tuning.NEUTRON_STAR_KICK_SIGMA_KMS`, [R1966] as `multistar-and-compact-systems.md` records it), so the velocity is a three-dimensional Gaussian about `Vc` and the chance that it exceeds `v_esc` is the noncentral chi tail of three degrees of freedom, in closed form (`survival`). The massive stars that leave neutron stars are born in the thin disc, `R exp(-R / L)` per unit radius, so the galaxy's share is the mean over 32 equal shares of that law.

The escape speed is the galaxy's own: the halo's Navarro, Frenk and White potential from the mass chain (halo mass, virial radius, concentration) plus the stars' mass as a Plummer sphere of the disc's scale length.

```
v_esc^2 = 2 ( V200^2 c / m(c) * ln(1 + x) / x  +  G M_* / sqrt(R^2 + L^2) ),  x = R / r_s,  r_s = r200 / c
```

## Numbers

| Galaxy | Escape speed at 8.2 kpc | Neutron stars ejected |
|---|---|---|
| Milky Way masses (5e10 of stars, published halo 1.15e12) | 539 km/s (the measured value is about 550) | 27% |
| Median halo of 5e10 of stars (2e12, the reference structure) | 632 km/s | 16% |
| Over 200 drawn galaxies | | 0 to 74%, median about 50%: the small galaxies lose most |

The old 20% to 40% is a test target for the Milky Way's published masses (27%, inside it); the reference structure's 16% is 20% under the range's low end, and the gap is shown in `tests/test_natal_kick.py`. The nature of the draws changes: the share was one uniform draw per galaxy and is now a function of its stellar mass, halo, circular speed and disc length, so a seed still gives one value.

## What follows

- The neutron star count of every place (`remnant_density`) and the ejected count (`ejecta_field`) read the share from the structure, as in slice 1; bound plus ejected is still what the galaxy formed.
- The black hole share cap (GEN.563 card 5, "Keep the cap") is the Milky Way-like galaxy's, the same for every galaxy: with a derived share a galaxy that ejects more than 59% of its neutron stars would have a cap of 0 and no stellar black holes. The choice is on Boss's card (the cap at the Milky Way-like share is the recommended option, and the one built).

## The ejected stars' speed

A kick speed `s` above the escape speed `v` at the birth radius leaves the neutron star with `sqrt(s^2 - v^2)` far out. The speed
law is the density of `|V + k|` (the noncentral chi of three degrees of freedom, `_speed_density`), integrated above `v` over the
same birth radii, and the galaxy's structure holds its 16 equal-share quantiles (`Structure.neutron_star_ejecta_speeds_kms`,
`ejecta_speeds_kms`). The ejecta population's density takes their mean (`ejecta_field.mean_speed_kms`: 386 km/s for the
Milky Way-like structure, against the 400 of the old 300 to 500 range) and each object draws a speed between two quantiles.
A galaxy that ejects none has none. A black hole's speed keeps the range (`EJECTA_SPEED_RANGE_KMS`), a **[stand-in]** with the 2%.

## Stand-ins and open items

- **[recalled]** The dispersion 265 km/s [R1966] is recalled in the repository's notes and not re-read (Boss, 1:17 PM EDT 2026-10-21: use what is included); one Maxwellian stands for what is probably two components **[stand-in]**.
- The black holes' 2% (`BLACK_HOLE_UNBOUND_FRACTION`, one population synthesis) is not derived: a momentum-conserving kick of `sigma x 1.4 / M` (37 km/s for 10 solar masses) unbinds almost none. It stays a **[placeholder]**.
- **[stand-in]** Birth in the thin disc only: the bulge's and the halo's neutron stars are born at their own radii, not averaged here. Disc scale lengths a shape fixes by hand keep the seed's drawn length for this share. 
- The potential is spherical and static, and the progenitor's speed is the circular speed with no asymmetric drift.

## Tests

`tests/test_natal_kick.py`: the survival formula matches a drawn kick, falls with the threshold and rises with the progenitor's speed; the birth radii divide the disc into equal shares; the escape speed falls outward and rises with the halo, and is 490 to 600 km/s for the Milky Way's masses; the Milky Way's share is inside 20% to 40%; a faster kick ejects more and a deeper well fewer; a seed draws the share its structure carries; small galaxies lose more than large ones.
