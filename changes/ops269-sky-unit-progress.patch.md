**Headline:** The bright-star backfill bar counts the sky tiers as they are drawn, and each sky share's counter is costed the way the stage was planned, so the ETA no longer reads in hundreds of days.

### Fixed
- OPS.269: a sky unit of the bright-star backfill is one piece, so its own count stood at 0 of 1 until the piece ended and the backfill bar credited nothing for hours of drawing (Boss's run: 1,088 of 7.00e8 after 4.5 hours, an ETA of 194 to 254 days). The sky draw now moves the unit's own count with the share of the tiers it has drawn, which the bar, the Queue page and the stage's ETA read.
- The counter of a sky share counted the visits of every ring's whole sphere, about three times the walk's visits that the stage's plan (and its bar total) costs; it is now costed by what the walk does, as the plan is.

generation-output: unchanged
