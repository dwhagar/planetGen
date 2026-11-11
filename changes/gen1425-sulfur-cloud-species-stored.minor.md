**Headline:** Rocky planets and moons now carry the species of the Venus sulfur-cloud set their air, light and ground allow (acid droplets, sulfur dioxide, carbonyl sulfide and the other gases, sulfur monoxide, sulfur aerosol, iron chloride and ferric sulfates), and the property sheet names them.

### Changed

- The Venus chemistry set of GEN.1425 (V1 to V10) had its rule written (`sulfur_cloud_set.chemistry`) but only the acid cloud and the sky colour band reached a stored value. A rocky planet or moon with an air now stores `cloud_chemistry`, the ids the chain gives it, and the property sheet shows them as "Sulfur-cloud chemistry" with the name and the status of each (observed in a planet, predicted only).
- It reads what the body already has: its finished air, its absorbed light and the star's ultraviolet-C share, its ground (iron needs a rocky surface with no global ocean) and the acid cloud its aerosols were drawn with, so the stored cloud depth and the species agree. Species are stored as present or absent only; none has an amount, because none has a source for one. In 100 sampled systems (G2V, K5V, M2V, F5V; 811 rocky bodies, 470 with an air) 303 bodies store a species; V3 is on nearly all of them, V2 and V4 on 14, V5 on 7, and the full chain (V1, V6, V7, V10) on 2.
- V8 (iron-PAH from comet dust) and V9 (sulfur-chlorine products) are single-source hypotheses and stay off (`sulfur_cloud_set.INCLUDE_HYPOTHESES`, stand-in CG3); V8 also needs a dust supply per planet, which the generator does not track. A rogue planet has no star, so the gates that need light stay closed for it; the gas-only species can still apply.
- New column `cloud_chemistry` (VARCHAR(48), NULL) on planets, moons and rogue planets, migration 0244. No existing value changes: the draws and the air are as before, so this adds a column to the rows the fill and population stages make.

generation-output: changed
generation-affects: fill, population
