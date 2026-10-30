**Headline:** A star's own spin sets how active it is, flares reach a world's dose by their size and the field's shield, and a world the flares scour loses its ozone, so UV, sky and life follow.

### Changed
- **A star's activity follows its own rotation (GEN.807).** A fast spinner keeps its saturated X-ray output longer and flares more, a slow one leaves saturation sooner (unsaturated activity falls as the Rossby number to the minus 2.7, Wright et al. 2011); the star's spin about the gyrochronology median is now drawn once and read by the activity, not only by the sheet.
- **The mixture draw reads the XUV exposure stored for the body's distance (GEN.806).** The exposure the atmosphere draw sees is the one `magnetism.update_exposure` stored, not a second computation.
- **The particle dose reads the flare slope and the magnetopause (GEN.805).** A shallow flare-energy slope puts more energy in the big flares and raises the dose; a field cuts it by its standoff.
- **A world the flares scour keeps a tenth of its ozone (GEN.808).** Where the flare irradiation index times the field's shield reaches 300, or the galaxy's lethal events pass 10 a gigayear, the ozone column is cut to 10 percent and the DNA-weighted UV at the ground follows (`ozone_loss_flag`, `uv_surface_index`).
- **The loss history follows the star's own activity (GEN.804).** The XUV decay the thermal escape integrates runs through the star's stored `log_lx_lbol` instead of the table's nominal decay.
- generation-output: changed

### Added
- **`physics/ion_escape.py` (GEN.804).** The stellar wind's ion loss through the open share of a dipole (`1 - sqrt(1 - 1/r_mp)`), tuned to the Venus Express and MAVEN rates, and `thermosphere.kept_shares(magnetopause_rp=...)` to follow it. Nothing passes the magnetopause yet: the field is drawn after the air until GEN.967 orders it first.
- **Stand-in rows WF1 to WF6** in the gap-filling plan (section 6.32) and the stand-in list.
