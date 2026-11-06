**Headline:** The ocean floors of water worlds and the ice lids above them now melt along the published melting curves of water ice, and a partly molten mantle loses viscosity by the laboratory's rate, so deep oceans, ice lids and the tides that heat them shift a little.

### Changed
- **GEN.1190, melting curves.** `water_layer.liquidus_pressure_pa` takes the melting pressure of ice III, V, VI and VII from the IAPWS release R14-08 equations in place of a table of nine points (a straight line between them, with ice VI's last slope standing in for ice VII above 2.2 GPa). The triple points move by up to 1.5 percent (the ice V and ice VI point from 620 to 632.4 MPa, the ice VI and ice VII point from 2,200 to 2,216 MPa), and above 355 K the curve is ice VII's own. The sheet's high-pressure ice rule reads 632 MPa to match. `hydrosphere.pressure_melted_shell_km` lowers the melting point of an ice lid's base along the ice Ih melting curve of the same release instead of 0.0074 K per bar, which fell short by almost a third at 200 MPa (252.3 K there, not 258 K).
- **GEN.1190, rheology.** The viscosity of a partly molten rocky mantle falls as exp(-30 melt) (the laboratory's diffusion-creep value) in place of exp(-26 melt); a partly molten ice takes exp(-22 melt), read from the material. Io's tide and the tide of every body with a mantle between its solidus and its critical melt fraction change.
- **GEN.1190, ice types.** Ammonia dihydrate's strength against water ice is the ratio of the 1993 creep law to the ice I law: 1.7 at 135 K, 0.91 at 143 K, 0.57 at 150 K, 0.18 at 170 K, 0.14 at 176 K, replacing three interpolated points. The phase of a water-ice body's centre follows the boundary lines of the water phase diagram (from a web page, tier Marginal Unverified) when a temperature is given, in place of recalled pressures. Nothing in generation reads the ice types yet.
- Stand-in rows TR1, TR4, TR5 and TR9 carry the read evidence (the Andrade ratio's tier falls to Marginal Unverified because its sources conflict for rock), IT6 and IT7 are replaced, and MC1 to MC5 are added in section 6.193.

generation-output: changed

TODO items advanced: GEN.1190.
