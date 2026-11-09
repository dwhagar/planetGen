**Headline:** A new generator epoch no longer empties the map's stored nodes: they follow the stored rows, so only what was regenerated is fetched again.

### Changed

- OPS.133: the key the Galaxy Map's node cache follows (`tile_base` in `galaxy_content_state`) no longer holds the running generator epoch. It holds the epoch each scatter was drawn under (`scatter_epochs`) instead; the plan's own epoch was already in it through the stored shape. A node is read from stored rows, so an epoch bump with nothing regenerated changes no node. A regenerated sector bumps its own `modified_at` and so its own nodes only, and a scatter drawn again under a new epoch still empties every node.
- Nothing stored changes. The first update after this lands empties the stored nodes once, as every epoch did before.

generation-output: unchanged
