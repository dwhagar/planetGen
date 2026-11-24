**Headline:** The one-off system page and the system API take the star count and the cap and refuse what the command line refuses, in its words.

### Added
- ADM.128: `POST /api/systems` and a regeneration recipe take `stars` (exactly this many stars) and `max_stars` (the cap on the count drawn), as `planetgen system --stars --max-stars` does. A count under 1, a count over the cap or over 64, a cap outside 1 to 64, one star with a binary feature and two or more stars with `binary_system` off are a 400 with the command line's message. The one-off system page already had both fields; it now says the command line's words for a count under 1, over 64 and a cap outside 1 to 64 instead of its own range wording. Tests cover each refusal on all three surfaces and show that the three messages agree.

generation-output: unchanged
