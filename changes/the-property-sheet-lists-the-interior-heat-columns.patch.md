**Headline:** A planet's or moon's page now shows the heat flowing out of its interior and the two draws behind it.

### Fixed
- **`test_property_sheet` failed on main (TEST.351).** GEN.699 added `internal_heat_flux_w_m2`, `mantle_eta0_pa_s` and `mantle_activation_kj_mol` to planets and moons without a row on the property sheet or an entry in `NOT_SHOWN`. They now have rows in the radiogenic group, after the Urey ratio.
TODO items resolved: TEST.351. Added: none.
