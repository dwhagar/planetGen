**Headline:** The pulsar share no longer divides by zero at a galactic radius beyond about 580 kpc.

### Fixed
- **`pulsar_radial_factor` raised `ZeroDivisionError` beyond about 580 kpc (GEN.622).** The pulsar profile and the disk's own exponential are now divided in logarithms, so the factor reaches its cap of 2 (or its asymptote) instead of dividing by an underflowed zero. A test covers 579 kpc, 581 kpc and 10,000 kpc.
