# Nightside collapse of a locked planet's air (GEN.1044, oddity W19)

Status: lock-free module, written 2026-10-30 (as the lane reported) by the Aux 3 lane. `src/planetgen/physics/nightside_collapse.py`; tests in `src/tests/test_nightside_collapse.py`; registry rows NC1 to NC6 (section 6.75 of the gap-filling plan). The wiring is GEN.1051, after the climate solver, with the bistable branch of W20.

## 1. What it answers

On a tidally locked rocky planet with thin air, the night side can get cold enough that the main gas freezes onto it, which lowers the pressure, which cools the night side more. The result is a world with its air gone to a nightside cap and a bare day side. The module takes the air's pressure and gas, the surface gravity, the starlight and the planet's mass and returns a `Collapse`: the share of the gas that stays, the pressure that stays, the pressure under which it collapses, whether it is collapsed and where the number came from.

## 2. Where the numbers come from

| Piece | Source | Registry |
|---|---|---|
| CO2 critical pressure, a cubic in the flux, for one and ten Earth masses | Wordsworth 2015 (arXiv 1412.5575), Eq. 43 and Table 2: 3D runs of a pure CO2 air, 0.2 to 3.0 times 1366 W/m2, 0.01 to 10 bar, albedo 0.2. 0.14 bar at Earth's flux. | NC1, NC2 |
| Night-side temperature under a thin air | The same paper, Eq. 10 and 25: `tau = kappa p / (g cos alpha)`, kappa 5e-5 m2/kg, the night side radiates `tau / (2 - tau)` of a quarter of the absorbed flux. 125 K at 0.1 bar. | NC3 |
| What stays | The vapour pressure of the gas (`physics/saturation.py`) at the night side's temperature | none (physics) |
| Calibration of the grey estimate to the fit | One factor on the night-side temperature (1.24 at Earth's flux) so the balance first loses gas at the fit's critical pressure; the retained share then rises continuously to 1 there | NC4 |
| Other gases | The same balance, uncalibrated; the paper leaves N2 and CO to later work | NC5 |

The 0.1 bar of the design note's table is the paper's own figure (collapse "for any stellar flux given p < 0.1 bar"); the 0.03 bar found by Joshi, Haberle and Reynolds 1997 was a grey-gas estimate and the paper explains the difference by windows in the CO2 spectrum and by a weaker coupling of the air to the night-side surface than grey gas gives. Leconte 2013 (listed in the work item) was not read: the paper found in the search is about land planets and moist bistability, not collapse.

## 3. What it does and does not claim

- Under the critical pressure the loss feeds on itself, so the air does not settle half way: what stays is at most the vapour pressure at the night side's temperature, a tiny share (0.1 bar of CO2 at Earth's flux keeps about a quarter, 0.05 bar keeps under one percent, Mars-like 600 Pa keeps none to the precision of the model). The retained share is the number to use; `collapsed` (under half) is a convenience.
- It is not asked for a planet that does not lock; the caller (GEN.1051) checks `regional_temperature.is_locked`.
- Air that is a mix of CO2 and a non-condensing gas is outside the fit. The paper says the pure-CO2 fit is conservative (a background gas lowers the collapse pressure). The module takes the main gas as the whole air.
- Outside 0.2 to 3.0 times Earth's flux the fit is held at its edge and the `Collapse` says `extrapolated`. At the paper's lowest flux the fit gives 2.9 bar while its text says collapse at every pressure it ran; the module follows the fit, and a Mars-like world is therefore marked extrapolated and collapsed.
- A night-side ocean, obliquity, eccentricity, clouds and the heat of condensation are not modelled.

## 4. Limits for the wiring step

The air a locked planet is given should be read as a pre-collapse air. Wiring should (a) call `collapse` with the main gas by mole fraction, (b) subtract the lost share from the air before the surface temperature is solved, and (c) store the cap as a night-side ice cover, so the day side's surface conditions come out bare and dry. Tests against a forced draw belong to the oddity harness (W-harness).
