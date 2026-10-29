**Headline:** An admin can generate a whole block from the Galaxy Map: the block's info card has "Generate this block", which starts the same job as `planetgen sectors --block` on the Generate page.

### Changed
- **Generate a block from the map (UX.127).** The info card of a block (3, 9, 27, 81 or 243 sectors a side) shows "Generate this block" to a signed-in admin. It asks first when the block is large or its size is not known, then posts the Generate page's new Block mode (`block_key`, `m.shell.band.wedge`) and lands on the job's progress. While another job is running the button is disabled and the card says which job and links to it. The Generate page itself takes the same Block mode as a post.
