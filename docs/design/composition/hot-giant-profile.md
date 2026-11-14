# The temperature-pressure profile of a hot giant (GEN.1577)

Status: built (new module, nothing calls it yet). Module `src/planetgen/physics/hot_giant_profile.py`, tests `src/tests/test_hot_giant_profile.py`. Built ahead by Lane 07 for Lane 14, who reviews it and wires it behind the cloud deck (`physics/hot_cloud_deck.py`, GEN.994).

## What it is for

The cloud deck of a hot giant (Cr, MnS, Na2S, ZnS, KCl) is where the planet's temperature profile meets each species' condensation curve (`hot_condensates.cloud_temperature_k`). The deck needs a profile `T(P)` over the pressures it searches, 1e-6 to 1e3 bar. The generator has no such profile above the 1 bar reference level of a giant (`surface_air` solves one temperature there), so this module gives one.

## The form

[R1715] equation 29, the analytic semi-grey irradiated atmosphere in radiative equilibrium, READ in the research (`atmosphere-composition-parts/03-infrared-greenhouse-thermal.md` section 5.1):

    T^4 = (3 T_int^4 / 4)(2/3 + tau) + (3 f T_irr^4 / 4)[2/3 + 1/(gamma sqrt 3) + (gamma/sqrt 3 - 1/(gamma sqrt 3)) exp(-gamma tau sqrt 3)]

with `tau = kappa_th P / g`. Inputs, all already built: the whole-globe equilibrium temperature (`climate_draws.equilibrium_temperature_k` with no internal flux, so `T_irr^4 = 4 T_eq^4` and `f = 1/4` returns it), the giant's own heat flux (`giant_draw.own_flux_w_m2`, GEN.816) as `T_int`, and the gravity. Two opacities: a constant thermal one and a visible one that is `gamma` times it.

Why this form: it is the research's own READ and coded-and-checked profile (the worked example below), it has the three features the deck needs (a cool skin, a deep isotherm, and an inversion when a strong visible absorber raises `gamma`), and it is a closed form with no solve, so it is cheap and draws nothing.

## Checked against

The research's worked example (HD 209458 b, kappa_th 1e-2 cm2/g, kappa_v 4e-3, T_int 100 K, T_eq 1469 K, g 1000 cm/s2, f 1/4): 1473 K at 66.7 mbar, 1647 K at 1 bar, 1648 K at 10 bar. The module reproduces 1472.6, 1647.4 and 1648.1 K (`test_the_worked_example_of_the_research`). The top equals the skin temperature; the bottom equals the deep isotherm; `gamma` large inverts the profile; the internal heat warms only the deep end.

## Stand-ins (registry 6.288)

| row | what | tier | why |
|---|---|---|---|
| HGP1 | `KAPPA_THERMAL_M2_KG` = 1e-3 m2/kg | T2, L0 | Guillot's own value for HD 209458 b; one value for every hot giant |
| HGP2 | `GAMMA_AT_2000_K` = 0.6, `gamma` as `0.6 sqrt(T_irr / 2000 K)` | T2, L0 | Guillot's `kappa_v = 6e-3 sqrt(T_irr / 2000 K)` over his `kappa_th` |
| HGP3 | constant opacity over 1e-6 to 1e3 bar | T6, L0 | above about 10 bar real collision-induced and electron opacity grow with pressure, so the profile underestimates the deep temperature; the lowest clouds (Cr, MnS) can sit there |

## What a later part adds

- A convective floor: the deep end follows the interior adiabat, not the isotherm (the radiative-convective joining of [R3783], section 5.3 of the research). The caller can floor the result with the adiabat of the giant's entropy until then.
- The two-band picket fence of [R3439] (section 5.2) for planets with a line-to-continuum split; its coefficients were not transcribed.
- Opacities from the atmosphere's own mixture (`line_opacity`, `greenhouse`) instead of HGP1 and HGP2.
- Wiring: Lane 14 calls `profile_function` for each hot giant and passes it to `hot_cloud_deck.deck`.
