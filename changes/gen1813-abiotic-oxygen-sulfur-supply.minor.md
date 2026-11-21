**Headline:** An abiotic-oxygen world's sulfur and nitrogen now come from the share its surface holds, so the sulfate and sulfide sinks act on it (a 600 K ocean world no longer keeps H2S whatever the sinks do)

### Changed

- GEN.1813: row `j5` of `mixture_evolved` set its sulfur and nitrogen per carbon atom of the delivered inventory, which skipped `supply` (the share of each element the surface holds, with the sulfate and sulfide sinks of `sulfur_balance` and the nitrogen's melt share). They are now per carbon atom of the delivered inventory times the share the surface holds, as the other outgassing rows read it. The carbon is the ocean's, so its own share is not read. Measured on the two 600 K ocean worlds of Lane 07's trace (seeds 109 and 155 on the GEN.1769 slice 2 base): Pidoe II's H2S goes from 473 Pa to none and its N2 from 642 Pa to 22 Pa; Azehelha I holds no H2S.

generation-output: changed
generation-affects: fill
