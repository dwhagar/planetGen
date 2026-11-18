**Headline:** Carbonyl sulfide, carbon disulfide and sulfur trioxide now freeze out of a cold air like every other gas with a vapour-pressure curve, so a body at 8 K no longer keeps thousands of pascals of them.

### Fixed

- GEN.1794: `saturation` gains rows for OCS, CS2 and SO3: the NIST WebBook Antoine equation of the liquid (a new `antoine` liquid form, carried from the triple point to the critical temperature of the gas table) and a Clausius-Clapeyron solid below the row's own triple point (134.33 K, 161.11 K and 290.1 K; the gas table has no triple point for sulfur trioxide), with the sublimation enthalpy the vaporisation plus the fusion, since NIST lists none. Stand-ins SV2 to SV5. Three bodies of 0.003 to 0.0046 Earth masses at 8 K held 9,700 Pa of OCS because `freeze_out` leaves a gas with no curve alone; at 8 K the pressures are now under 1e-150 Pa.
- The gases the equilibrium solve can carry that still have no curve are named with a reason under `no_law` of `saturation_pressure.json` (`saturation.reason_without_law`; stand-in SV6): H, Cl, S, SO, PO, PS, PN and P2 are fragments of the hot equilibrium, and CH3Cl, HCOOH, N2O, NO2 and PH3 are traces with no curve filed yet. `tests/test_saturation.py` fails on a carried gas with neither a curve nor a reason.
- Carried over unchanged: the 10 kPa retention cap of the d6 row (GEN.1558 retires it separately).

generation-output: changed
generation-affects: fill
