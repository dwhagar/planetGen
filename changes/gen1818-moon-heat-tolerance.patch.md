**Headline:** The heat a moon's tide puts into a planet no longer holds the steps of a stiff moon history: its error tolerance is a share of the planet's spin energy, so the seed 9101 history of the aquifer tests finishes again on Train 52.

### Fixed

- The last state of the moon-orbit integration is the heat the moons' tides dissipate in the planet, an integral of a noisy power. Its absolute tolerance was a fixed 1e20 J, so at 1e27 J the error scale was `RTOL` of the heat (1e22 J), inside the noise of the rheology's power, and the steps of a stiff history were held there: on Train 52 the history of the 5e23 kg moon of seed 9101 ran every integrator out of evaluations (LSODA, three BDF runs and Radau) at a fifth to a third of its span with the orbits smooth, and the planet came out with a tide flux of 2.5e7 W/m2 again (the GEN.1815 symptom). The tolerance is now `HEAT_ATOL_SHARE` (1e-4) of the planet's rotational kinetic energy at the start of the stretch (stand-in tier T6, tweak level L1). The history finishes; the moon falls to its planet and the planet's tide flux is 8e-08 W/m2. Orbits, spins and fates are as before; the heat of a finished history moves by a few percent where it is a small share of the spin's energy (GEN.1818).

generation-output: changed
generation-affects: fill
