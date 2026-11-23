**Headline:** When the site is slow to answer after the update's reload, the update log now shows what the site's processes and Apache's error log say.

### Fixed
- OPS.299: after the reload the health check got no answer for about 60 s (Boss's update at 6:05 PM EST 11-23) and the log could not say why. About 25 s in, and when it gives up, the update now prints the site's processes (age, processor time, memory) and the last lines of Apache's error log, and an answer that came late says how late. The cause is not yet known; the next slow start names it.

generation-output: unchanged
