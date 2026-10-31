**Headline:** The physics of a white dwarf made by a close binary is in, ready to wire: the helium core a giant leaves when it fills its lobe.

### Changed
- **The helium-core white dwarf of a close binary (GEN.996, physics module only).** A star under 2.3 solar masses that fills its Roche lobe on the red giant branch leaves the helium core it then had, 0.18 to 0.45 solar masses, and not the 0.55 or more the single-star relation gives. `physics/helium_core_white_dwarf.py` holds Rappaport et al. 1995's core mass and radius relation (read second hand), the helium flash at a core of 0.45 solar masses and the lobe of Eggleton 1983. With the generator's own period law and mass-ratio law about 14 percent of a solar-type star's companions sit where the lobe is reached in that window, which makes about 6 percent of the white dwarfs of 1 to 1.5 solar mass progenitors helium-core, against 7.4 percent in the 20 pc sample; no rate is set. Nothing calls it yet: the stored tree records a white dwarf's mass from its progenitor, so the draw waits for the member record of GEN.1006 P3. Four stand-ins are registered (HC1 to HC4, section 6.160).

generation-output: unchanged

TODO items advanced: GEN.996 (module built, wiring waits for GEN.1006 P3).
