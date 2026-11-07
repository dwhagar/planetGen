**Headline:** The radiative-convective model can now be asked for the saturated clear-sky air the published runaway-greenhouse limits are worked for, and it runs away there at 300 to 315 W/m2 of outgoing flux (published 282 to 310): the missing plateau was the dry-aloft humidity profile of the Earth, not the water opacity.

### Added
- GEN.1194: `radiative_convective.solve(..., saturated=True, clear_sky=True)` holds the relative humidity at the surface's value up to the cold trap and leaves out the cloud's longwave effect (stand-in IR17). On 1 bar of nitrogen with 40 Pa of carbon dioxide the surface answers a watt with about 1 K up to 300 W/m2, with 1.3 K at 315, and with over 2 K past 320 (322 K at 300 W/m2, 343 K at 315, 375 K at 330). The Manabe and Wetherald profile that the generator keeps (IR7) gives 310 K for the same 330 W/m2, which is why the model looked as if it had no plateau.
- Two tests on the real tables hold the elbow and the finding that it belongs to the humidity, not to the opacity.

### Changed
- Nothing the generator makes: the greenhouse table is not refitted and the defaults of `solve` are as they were. The refit of the hot dense airs (GEN.1310) is where the option's findings would reach generated worlds, if the humidity of the hot ocean worlds is changed there.
