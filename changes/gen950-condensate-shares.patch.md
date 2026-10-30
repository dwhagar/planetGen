**Headline:** The water share of a disk's solids by temperature is now computed from Lodders 2003 (not yet used by the generator).

### Changed

- New library `physics/condensate_shares.py`: the condensation temperatures of water, ammonia and methane ice and the mass each condensate makes (Lodders 2003 Tables 10 and 11), the water share of the solids by disk temperature and by distance, and the ice lines of a star. Nothing in generation reads it yet; the wiring step (snow line temperature, ice boost, water range of an icy body) follows with the epoch train.
