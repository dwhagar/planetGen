**Headline:** The dust that dims a galaxy's stars can now be worked out from the galaxy's own gas, so a small galaxy no longer inherits the Milky Way's dust.

### Changed

- `physics/dust_model.py` (VIEW.43) works out the dust of a galaxy from the gas of every ring of its disk (GEN.469), the metals of that gas (GEN.1269 and the dust-to-gas law of [R5102]) and the thickness a gas layer has in the disk's own gravity. It gives two tables, the extinction through the whole layer and the layer's scale height, against radius. The extinction per parsec is the column over twice the height times sech squared of the height above the plane. One constant comes from a measurement (2.21e21 hydrogen atoms per magnitude, [R5123]) and one is tuned (a half of the gas column dims the stars a catalogue can use).
- In a galaxy shaped like the Milky Way it gives 1.06 magnitudes a kpc in the plane at the Sun and 0.24 magnitudes through the layer, and about 39 magnitudes from the Sun to the Galactic Center; the fixed law gave 2.3. A galaxy of the test galaxy's size gives 0.14 magnitudes a kpc where the fixed law gave 9.6.
- `sky_photometry.dust_av_per_pc` takes the table form (`profile=`) beside the fixed law; `av_per_kpc_at_sun=0` still means no dust. Nothing reads the profile yet: the sky route, the glow band and the tier ledger take it in VIEW.44, so the sky is unchanged by this note.
- Five stand-in rows, DU1 to DU5 (section 6.188 of the registry), and five references ([R5123] to [R5127]).

generation-output: unchanged
