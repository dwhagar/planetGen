**Headline:** The sulfuric acid cloud has an infrared opacity (Mie droplets from the HITRAN acid indices) and the minor gases have a far wing in carbon dioxide air, built and tested; no world changes yet.

### Changed

- **Infrared optics of the acid cloud (GEN.1020, lock-free step).** `cloud_infrared` turns the acid cloud's depth at 550 nm and droplet radius into an optical depth on each of the 130 infrared bands, from Mie spheres (`mie`) of 76 and 81 percent sulfuric acid (the HITRAN aerosol indices of Myhre and others 2003, tabled by `scripts/build_acid_cloud_optics.py`), spread between the cloud base and top in the layers of a column. `co2_far_wing` gains the far wing of water, sulfur dioxide and carbon monoxide lines in a carbon dioxide air. In the clear-sky Venus-like test air the sourced wing and the cloud reach 757 K (measured 735 K). Nothing calls them yet; `line_opacity` still has the flat IR3 wing until the wiring step. Registry rows CF1 to CF5, FW6 and FW7.
