# The atmosphere colour

Every planet or moon with a solid surface and air has an **atmosphere colour**: the colour of the light its air scatters
down to the ground, shown as a chip and a word ("pale blue", "cream", "orange") in the **Atmosphere Colour** row of its
description and as the glow ring on the System Map. Sky colour and atmosphere colour are the same thing here.

It is the *absolute* colour: it is not adapted to the star the way an eye would be, so the same air is pale blue under a
Sun-like star, cream under a red dwarf and orange under a very thick carbon dioxide atmosphere. The chip is a hue (its
brightest channel is always full); the brightness is kept separately, as the sky's luminance over the star's.

## What it is worked out from

| Input | Where it comes from |
|---|---|
| The star's spectrum | A blackbody of the star's effective temperature from 4500 K up; cooler stars are bent until their colour matches the BT-Settl atmosphere offset (`physics/photometry.py`, `physics/atmosphere_colour.py`) |
| The air | The stored mixture of the planet or moon (one row for each gas in `planet_gases` or `moon_gases`, the balance gas of its class included; see [atmosphere-data-model.md](atmosphere-data-model.md)) and its gravity (`physics/air_column.py`) |
| Rayleigh scattering | The mixture's molecular cross sections at each wavelength (`physics/rayleigh.py`); about 0.097 optical depth at 550 nm for 1 bar of Earth-like air |
| Ozone | The stored ozone column in Dobson units (`ozone_column_du`) and the Chappuis band cross section (`physics/data/ozone_chappuis.json`) |
| Haze, dust, cloud | Optical depth at 550 nm, single-scattering albedo, asymmetry and wavelength law of each layer (`atmosphere_colour.Aerosol`); stored with the body (`tau_haze`, `tau_dust`, `tau_cloud`, `cloud_cover`, `aerosol_draw`, `physics/aerosols.py`); the colour is worked out once at generation from them and stored as `atmosphere_colour` (`physics/body_colour.py`) |
| The ground | Its reflectance per wavelength; a stand-in albedo of 0.15 until the ground colour is stored |

## The method

1. The star's light goes through the air in layers: the gas (Rayleigh scattering and the gas table's absorption bands), the ozone
   layer at its Chapman-peak pressure, and the aerosol layers, as a plane-parallel atmosphere with the sun 30 degrees from the zenith.
2. Thin air (optical depth at 550 nm under 0.2) takes the single-scatter hue of the zenith sky; thick air (over 0.4) takes the
   delta-Eddington two-stream diffuse light at the ground; between them the two are blended so the colour never steps. The ozone and
   the ground's bounce enter as the ratio of two-stream runs with and without them.
3. The spectrum goes through the CIE 1931 colour matching functions (the [R4788] fit), the sRGB matrix, the
   brightest channel is scaled to 1, then the sRGB transfer curve gives the hex colour.
4. The word comes from the hue and saturation of the hex colour. When the aerosols alone move the colour by 0.03 or more in CIE xy
   the label says "(aerosol-dominated)".
5. **The thick-air limit.** Single scattering would saturate to red in very thick air; the two-stream light reddens smoothly instead,
   so 9 MPa of carbon dioxide gives an orange, never a clipped red. The direct beam dies as exp(-tau / mu0) and the sun disc is a faint
   deep red.
6. **Gas and ice giants** have no ground, so the colour is read at a reference level, the lower of 1 bar and the cloud top, with a
   cloud deck (or the dark deep gas) below (`physics/giant_colour.py`). Methane, ammonia and water use the gas table's coarse bins, so
   such a giant is marked partial quality until the band strengths are read. Descriptions show no row for a giant yet (it needs its
   stored cloud); the System Map glow uses the air above 1 bar.

## What it does not include

Sunsets, twilight and the day strip, the view of the sun's disc and the sky in other directions, and the colour of a planet seen from
space belong to the view from a planet and the planet's own appearance; they are not part of this colour. The glow ring is the limb
of the air, not the surface.

## Sources

Rayleigh depth: [R0575: Table 4], [R4138] (532 nm cross sections), the gas table's King factors. Ozone column:
[R2466]. Chappuis band: Brion, Daumont, Malicet and Jenouvrier data in the MPI-Mainz UV/VIS Spectral Atlas (CC BY 4.0, see
`DATA-LICENSE.md`). Two-stream: [R2199, R3013]. Colour: [R4788]; IEC
61966-2-1 (sRGB). Cool stars: BT-Settl chromaticity offsets. Dust: [R4743] (digitized, see the design note). Giants:
[R0024] for the cloud top; [R0622] as the comparison. Every number and the rule that picked it
is in [docs/design/composition/atmosphere-colour.md](../design/composition/atmosphere-colour.md).
