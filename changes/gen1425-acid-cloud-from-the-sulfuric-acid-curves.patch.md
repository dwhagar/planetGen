**Headline:** A world with sulfur dioxide, an oxidiser and water now has a sulfuric acid cloud whose base and strength come from the real acid and water curves of the air, not a fixed 367 K.

### Changed

- **The acid cloud base and strength come from the sulfuric acid curves (GEN.1425 part B, GEN.1424, GEN.1019).** The cloud sits at the warmest level of the air's dry adiabat where its own acid vapor and water condense (`sulfuric_vapor`), so a wetter air gets a weaker acid and a warmer base and a drier one a stronger acid and a cooler base. The fixed 367 K base (`sulfuric_cloud.json`) is gone. `aerosols.draw_aerosols` adds the acid cloud's depth (20 to 30) to the condensate cloud's, and its cover (the whole sky) is the larger of the two; it uses the stored `aerosol_draw`, so no column and no draw is added and airs without sulfur are unchanged to the digit. New epoch, `--affects fill`; no migration.
- **The acid medium of life Classes 5 and 9 reads the same strength (GEN.888, GEN.1424).** `acid_medium` takes the acid's strength from the Gibbs-energy model of `sulfuric_vapor` and no longer from the water-activity line fitted through the Venus cloud levels (three numbers gone from `acid_medium.json`). A very wet Venus-like air still has only a dilute acid; one with a tenth of that water keeps a cloud layer.
- **The cloud's ultraviolet absorber colours the sky (GEN.1425 V10).** A body with an acid cloud has the measured Gaussian absorber band (center 0.34 micron, width 0.14, moved inside its stated uncertainties by the stored `aerosol_draw`) in the cloud's single-scattering albedo, so the sky is a little yellower where the cloud is. Nothing is stored for it.
- **For the infrared model.** `sulfuric_cloud.infrared_cloud` gives `radiative_convective.solve(cloud=)` the depth, droplet radius, base temperature and strength of the air's own cloud; `cloud_infrared.layer_depth` takes the base temperature from the caller.

generation-output: changed
