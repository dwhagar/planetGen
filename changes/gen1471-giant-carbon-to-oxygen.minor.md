**Headline:** The carbon-to-oxygen ratio of gas giants now centers on the retrieved giants' 0.57 instead of 0.33, with the right share under 0.35 (GEN.1471).

### Changed

- A gas giant's C/O comes from the disk it formed in (GEN.1363), times a spread. Before, the spread was divided by 0.57, the center of the recipe's core, but the draw's own median is 0.51 (the core is cut at 0.8), so every giant came out about a sixth too low; and the recipe's spread was added in full on top of the disk's own scatter. It is now divided by the draw's median and taken to the half power.
- The envelope's C/O is set with a share of the planet's bulk heavy elements (0.06, stand-in DSC9), because the bulk metal is only an upper limit on the atmosphere's [R4350: Abstract] and gas giants show low water [R4663: Abstract]. Before, the whole bulk (18 times the star's for the median giant, 90 percent solids of C/O 0.44) went in, so the C/O was the solids' 0.44 times the star's whatever the amount.
- In plain numbers (G2V giants; the share was chosen on the first 347 giants of 1,000 systems and checked on 331 fresh ones of 1,000 more): median C/O 0.33 before, 0.57 on the first sample and 0.59 on the fresh one after, against 0.57 for the 69 retrieved giants; under 0.35: 53 percent before, 16 and 14 after, against 17; above 0.8: 7 percent before, 22 and 26 after, against 4 to 9. The upper tail is now too high: a giant that took mostly gas has the water-frozen gas C/O of 0.89 and the model has no drift and evaporation of water to spread it; that is noted as open and not tuned away.
- The metal inventory of the envelope is unchanged; only the carbon and oxygen split uses the share. Registry rows DSC7 (re-based), DSC8 and DSC9 in section 6.213.

generation-output: changed
