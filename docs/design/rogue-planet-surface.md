# Rogue planet surface conditions

A rogue planet has no star. Its surface is set by its own heat alone,
leaking out of its interior and radiating into the 2.725 K cosmic
background. This note records the model `planetgen/physics/rogue_surface.py`
uses (schema v48), taken from Boss's research of 2026-10-01, and the
defaults picked where the research gives no number. The constants are
in `program_constants.py` under "Rogue planet surface conditions".

## The four steps

1. **Internal heat flux.**
   - A rocky rogue: radioactive decay, `F = (M_mantle / 4 pi R^2) sum
     C_i H_i e^(-lambda_i t)`, with Earth's mantle U, Th and K (Turcotte &
     Schubert) as they stood at the rogue's age, times an abundance drawn
     0.5-2x Earth's. Leftover formation heat is added in Earth's
     proportion (Urey ratio 0.67), so an Earth twin at 4.5 Gyr flows
     0.087 W/m^2, Boss's 0.08-0.1. A rogue that kept a moon adds tidal
     heating, 1e-4 to 1e-2 W/m^2.
   - A giant or brown dwarf: Kelvin-Helmholtz cooling, `L ~ t^-1.3`
     (Burrows et al.). In mass, a power law from Neptune (T_int 53 K) to
     Jupiter (5.4 W/m^2 at 4.5 Gyr), then from Jupiter to Burrows &
     Liebert's brown dwarf law at 13 Jupiter masses, then that law.
     T_eff is capped at 2,800 K: the law assumes an old, shrunken body
     and overshoots for very young, still-puffed-up brown dwarfs.
2. **Effective temperature.** `sigma T_eff^4 = F_int + sigma T_CMB^4`.
3. **Atmosphere.** A rocky rogue kept a primordial hydrogen envelope with
   chance 10% (terrestrial bin) or 50% (sub-Neptune bin), base pressure
   10-1,000 bar or 100-10,000 bar, and keeps it only if its Jeans escape
   parameter for H2 at T_eff is above 15 (at these temperatures it always
   is). Without one, any air freezes onto the surface; a rogue of 0.3
   Earth masses or more had air, and its nitrogen frost leaves a trace
   vapor pressure (Clausius-Clapeyron from N2's triple point: about 6 Pa
   at 40 K, nothing at 20 K).
4. **Surface.**
   - Envelope: down the dry adiabat `T = T_eff (P / 0.1 bar)^(2/7)` to its
     base (Stevenson: 30 K at 0.1 bar becomes ~417 K at 1,000 bar).
   - No envelope: the surface sits at T_eff.
   - Water: a rogue is water-rich with chance 35% (terrestrial) or 80%
     (sub-Neptune), water 1e-4 to 10% of its mass. Below 273 K at the
     top, an ice lid of thickness `D = (A / F) ln(273.15 K / T_top)`
     (k = A / T, A = 567 W/m) covers a liquid ocean, unless D is deeper
     than all its water: then it is frozen to the rock. Between melting
     and boiling (at the envelope's pressure, capped at the critical
     point) the ocean is liquid at the surface.
   - A giant has no surface: its temperature is given at 1 bar, on the
     gray profile `T^4 = (3/4) T_eff^4 (tau + 2/3)` above the photosphere
     and the adiabat below it, with one Rosseland opacity fixed so
     Jupiter's photosphere is at 0.3 bar (1 Jupiter mass at 4.5 Gyr:
     T_eff 99 K, ~140 K at 1 bar).

## Regimes

`surface_regime` is one of `bare-rock`, `frozen-atmosphere`,
`ice-shell-ocean`, `ice-world`, `hydrogen-envelope`, `gas-giant` and
`brown-dwarf`. `has_internal_heat` (shown as "Geologically Active") is no
longer a 40% roll: a rocky rogue is active when its heat flow reaches
0.04 W/m^2 (between Mars and Earth); a giant always is.

## Simplifications

- The water layer's depth uses liquid water's density and ignores
  high-pressure ice at the bottom of a deep ocean.
- The ice lid's base is at 273.15 K whatever the pressure.
- A rocky rogue's leftover heat stays in Earth's proportion to its
  radioactive heat at every age, rather than following its own cooling.
- A rogue's radius is not changed by its age.

## References

Abbot & Switzer 2011 (ApJL 735:L27); Burrows et al. 2001 (Rev. Mod.
Phys. 73:719); Burrows & Liebert 1993; Guillot 2005 (AREPS 33:493);
Lissauer & de Pater 2019; Pierrehumbert 2010; Stevenson 1999 (Nature
400:32); Turcotte & Schubert 2014 (Geodynamics, 3rd ed.).
