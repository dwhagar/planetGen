**Headline:** How long a star stays X-ray saturated is now tabulated by mass and initial rotation from Johnstone et al. 2021 (not yet used by the generator).

### Changed

- New library `physics/saturation_time.py`: the age at which a star falls out of X-ray saturation by mass (0.1 to 1.2 solar masses) and initial rotation percentile, read off the published model figure. Cool red dwarfs stay saturated for gigayears, Sun-like stars for tens of millions of years. Nothing in generation reads it yet; the wiring step (the activity table and row p6) follows with the epoch train.
