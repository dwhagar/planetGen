**Headline:** The methane and nitrogen clathrate curves are checked against published equilibrium pressures, not only against their own constants.

### Added
- A test holds the CH4 and N2 dissociation curves (`physics/clathrate_brine.py`) to the lower quadruple point of the hydrates, which PetroWiki gives as 2.56 MPa for methane and 14.3 MPa for nitrogen near 273 K (reference R5375, secondary, after Sloan and Koh). At the ice point the fits give 2.47 and 14.79 MPa, 3.6 percent under and 3.5 percent over; the test allows 5 percent, about 2 K of the page's "close to 273 K". The CO2 curve was already held to its measured quadruple point. No generator or stored-data change.

generation-output: unchanged
