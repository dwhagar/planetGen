# Neutron star density: what the black holes leave of the massive stars (GEN.554, slice 1)

Status: built (GEN.554 slice 1). Mass islands and the kick step are later slices (below). Code: `galaxy/remnant_density.py` (the black hole density of `black-hole-density.md`, now taking the remnant `kind`), `galaxy/structure.py` and `galaxy/structure_draw.py` (the galaxy's neutron star unbound share), `galaxy/ejecta_field.py`, `galaxy/profile.py`, and the scatter, reach and fill sites that read the kind factors.

## What changed

Before: every galaxy rolled neutron stars at a fixed 4.0e-3 per star (`tuning.PHENOMENON_RATE_PER_STAR`, from Boss's research density of 1e9 in 2.5e11, GEN.100), whatever its metals. A black hole needs a massive star that fails to explode (`black-hole-density.md`), and that star would otherwise have made a neutron star, so the two counts were not independent and the old fixed rate could not follow the black hole draw.

After: every star of 8 to 150 solar masses ends as a neutron star or a black hole (the IMF share, 0.637% of the stars born). A place's neutron stars per star are the black hole sum with `1 - F_p` in place of `F_p`:

```
r_NS(place) = IMF_SHARE x (1 - u) x sum over populations p [ w_p(place) x C_p x (1 - F_p(place)) ]
```

`F_p` is the failed-supernova share of `black-hole-density.md` (the same two drawn anchors per galaxy), `C_p` the same completeness share, and `u` the galaxy's neutron star unbound share (the 0.2 to 0.4 draw that `ejecta_field` already used, now a field of the galaxy's structure so the one value feeds both the bound count and the ejected count). The ejected count is `formed x u`, so bound plus ejected is what the galaxy formed. A galaxy that fails more massive stars (a higher black hole share) makes fewer neutron stars, and the black hole cap of the GEN.563 card ("Cap f so the sum fits") now holds by construction: neutron stars plus black holes formed never exceed the massive stars.

## Numbers

Middle draw, Milky Way-like shape: a standalone sector rolls 4.08e-3 neutron stars per star against Boss's research density of 4.0e-3 (2% apart); the galaxy mean bound is 3.95e-3. The old 4.0e-3 is a test target (within 25%), no longer an input. Black holes at the standalone rate are unchanged (5.8e-4 per star). Only the neutron star counts move, and only through the galaxy's drawn shares, so `generation-output` changes.

## Stand-ins and open items

- The shares are the black hole draws (`black-hole-density.md` stand-ins: uniform across the verified range, linear in Z). No new tunable.
- Magnetars still use `MAGNETAR_PER_NEUTRON_STAR` times the tuning rate (a stand-in); they follow the young stars, not the metals.
- Later slices: the neutron star mass function with its islands (GEN.554 mass part, card 7 of the black hole note) and the natal kick step (position, velocity from the progenitor's orbit). The cap stays under the recommended option until Boss's card answer says otherwise.

## Tests

`tests/test_neutron_star_density.py`: neutron stars and black holes add up to the formed remnants per population, bound plus ejected is what the galaxy formed, a higher black hole share leaves fewer neutron stars, a metal-rich place makes more, the solar-circle rate is within 25% of 4.0e-3, the profile's factor is the place's and within the bound, and only the stellar black holes and plain neutron stars follow the metals.
