**Headline:** A cloud and a haze are painted in the colour of what they are made of: water white, sulfuric acid pale yellow, carbon dioxide ice grey-white, tholin orange, dust tan.

### Fixed
- MAP.418: every cloud was painted white and a haze only softened the picture. The look record now carries `cloudColor` (from the material that makes most of the cloud: the largest condensate by mass over particle density, or the sulfuric acid cloud where its depth is over the rest) and `hazeColor` (tholin and dust blended by their stored depths), worked from the stored air and aerosol values with no new column (`physics/aerosol_material.py`), and `static/solidpaint.js` paints the cloud in its colour and tints the haze deck toward its colour.
- MAP.416: the census of the data behind it (5,000 systems, `reports/map416/`) shows clouds are stored on 47 percent and haze on 54 percent of the bodies with air, so the thin look was the painting, not the data.

generation-output: unchanged
generation-affects: none
