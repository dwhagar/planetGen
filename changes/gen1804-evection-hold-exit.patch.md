**Headline:** A moon the evection resonance holds no longer keeps an eccentricity its tide could not leave it, so moons with tides of a million watts a square metre are gone.

### Fixed

- **A caught moon was held at e 0.5 for the whole age, though its tide cut that in years.** The evection resonance (GEN.1750) held a moon at an eccentricity of up to 0.5 and kept it to the star's age; the exit from the resonance (the tide taking back more than the resonance gives) was documented as not modelled. A moon of an Earth-like planet at 3.8 planet radii in that state dissipated 8.7e7 W/m2 (Io is 2) at a closest approach inside its own fluid Roche limit. A caught moon is now held only while the planet's tide takes back less of the resonance's growth than it gives (`evection.damping_share`, `evection.sustained`) and the planet's spin can pay the tide the hold dissipates over the rest of the age (`evection.hold_power_w`, `evection.affordable`); a held moon keeps the ceiling less that share. The damping rate is the one that already wears a captured moon's orbit (`orbit_elements.circularisation_time_yr`).
- On 200 systems (224 moons): the strongest moon tide falls from 8.7e7 to 0.25 W/m2, moons over 1 W/m2 from 4 to 0, moons whose closest approach is inside their fluid Roche limit from 3 to 0; no moon is lost, 5 moons change their orbit (the released ones drift on out to 11 to 22 planet radii as the tide carries them), 10 stay held.

generation-output: changed
generation-affects: fill
