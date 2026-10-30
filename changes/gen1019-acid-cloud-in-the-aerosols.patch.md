**Headline:** A world whose air holds sulfur dioxide, an oxidiser and water now has a sulfuric acid cloud, 20 to 30 deep, that covers its sky.

generation-output: changed

### Changed

- **Sulfuric acid cloud in the aerosols (GEN.1019).** `aerosols.draw_aerosols` adds `sulfuric_cloud.of_draw` to the condensate cloud: where the air has at least a part per million of sulfur dioxide, an oxidiser and a trace of water and cools to 367 K going up (the cold top found from the body's greenhouse factor, `aerosols.skin_of`), the stored cloud depth gains 20 to 30 and the cloud covers the whole sky; the droplet radius (0.3 to 3.85 micron) and the depth come from the stored `aerosol_draw`, so no column is added and no draw is added. Airs without sulfur are unchanged to the digit. New epoch, `--affects fill`; no migration. `sulfuric_cloud.of_body` gives the cloud and radius of a stored body for the infrared model (`radiative_convective.solve(..., cloud=)`).
