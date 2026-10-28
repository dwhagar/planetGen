**Headline:** The star glow now fades out much sooner, so blurred stars no longer leave huge faint halos over most of the screen.

### Fixed
- **Star glow halos grew huge when stars were blurred (MAP.330).** The glow of a block of unresolved stars began to fade only when the block was 100 pixels across and was drawn up to 220 pixels in radius, so faint halos covered most of the screen. It now begins to fade at a radius of 30 pixels, on a steeper curve, is nearly gone at four fifths of the way and is zero by 130 pixels, and a halo is never drawn larger than that.
