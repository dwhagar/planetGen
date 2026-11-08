**Headline:** A neutron star's page now shows the nuclear pasta layer at the bottom of its crust, the five shapes the nuclei are predicted to take, marked as predicted and never observed.

### Changed

- GEN.1427: `neutron_star_crust.PASTA` (data in `neutron_star_crust.json`, block `pasta`) holds the layer below the inner crust: from 0.2 of the normal nuclear density (0.165 nucleons per fm^3) to the bottom of the crust table, about 100 m thick, holding up to half the crust's mass, with the shapes in order of rising density (gnocchi, spaghetti, lasagna, anti-spaghetti, anti-gnocchi). The loader refuses a block with fewer than three shapes, a start above the crust-core boundary or a mass share above 1. The neutron star page's Crust table has a new row, "Nuclear pasta", that says it is predicted, not observed. Stand-ins PS1 to PS3 in section 6.232; sources [R5165], [R5166], [R5167].
- Nothing stored changes: the layer is a described fact read when the page is built, so the generated systems are the same as before.

generation-output: unchanged
