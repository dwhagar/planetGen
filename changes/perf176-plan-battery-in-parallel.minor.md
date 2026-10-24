**Headline:** A galaxy plan takes about half as long, because the digest it stores with the galaxy is worked out alongside the plan.

### Changed
- **A new galaxy's plan no longer stops for the generator battery (PERF.176).** The plan records a digest of the generator with a new galaxy, and worked it out in the middle of storing the shape: 8 to 20 seconds of a 15 to 39 second plan. It is now worked out in a process of its own while the rest of the plan runs; the plan waits under a second for it. The digest stored is the same.
