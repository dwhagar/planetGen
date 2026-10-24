**Headline:** Two small parts that let the generator call the atmosphere mixture draw are in: the description of a body it reads, and what an ocean takes out of the drawn air.

### Added
- **Body-to-mixture description (GEN.756).** `physics/mixture_subject.py` builds the quantities the mixture draw reads (mass, escape speed, age, insolation, equilibrium and layer temperatures, XUV exposure, snow-line and runaway-greenhouse ratios, zone and host type) from a generated planet or moon and its star; the climate loop passes in the surface temperature and pressure and whether an ocean can condense.
- **Ocean water cap and sulfur uptake (GEN.757).** `physics/mixture_sinks.py` holds the water vapour of a mixture over an ocean at its saturation pressure and dissolves sulfur dioxide and hydrogen sulfide by Henry's law. The generator does not call either yet (GEN.411 does), so no generated value changes.
