**Headline:** Every pair in a multi-star system now carries an interaction state worked out from Roche-lobe geometry.

### Changed
- New `planetgen/physics/roche.py` gives the Eggleton Roche-lobe radius, the filling factor of a star at its closest approach, and the separation at which a star just fills its lobe.
- Each pair the member-tree generator draws is classified as detached, semi-detached, contact, accreting, post-common-envelope or merger-product from the radii of its members at the system's age, and the mass gainer of a semi-detached pair is flagged. A pair never overfills its lobes by more than the model allows, and an eccentric orbit counts at its closest approach.
- Members carry the flags `mass_reversed`, `mass_gainer`, `merger_product`, `be`, `fragile` (the members of an ultra-wide system) and `theoretical`; a merger product gets a rejuvenated effective age while the real system age is kept.
- Three new member kinds exist in the model and the labels: stripped helium star, Thorne-Zytkow object and Wolf-Rayet star.
- The states are not stored in the database yet; that follows with the schema change (GEN.260). The exotic-star draws and the rate alignment stay with GEN.270.

TODO items resolved: none (GEN.269 stays open for storage on `scatter_nodes` and `orbit_nodes`) / Added: none
