**Headline:** The consolidation sweep's contradictions and prerequisite cycles are fixed in the plan.

### Changed

- MAP.189 and MAP.190 no longer need each other (MAP.190 first); GEN.271 no longer waits on GEN.253; GEN.327 waits for GEN.329, which owns `runaway_origin`; GEN.396 waits for GEN.407 and rides its epoch.
- GEN.400, GEN.404, GEN.405, GEN.370 and GEN.144 split the Sun-radius and central-mass work once: GEN.404 and GEN.405 no longer point at each other, GEN.406 (research) comes first and GEN.405 waits for it, GEN.370 follows GEN.400 and GEN.404.
- OPS.37 only shows the generator version and OPS.15 compares the battery digest; DB.58 accepts the moderator role; NAV.59, NAV.19 and TEST.180 lose text another item owns; DB.56, GEN.346 and GEN.143 use one name (`orbital_events`, written by GEN.197); MAP.121 steps as a camera move after MAP.227 and MAP.228.
- GEN.415, GEN.430 and GEN.446 wait for GEN.221 so the describer sentences are written once; GEN.440 lists the star kinds and flags it reads; MAP.236 follows MAP.240, UX.115 needs MAP.243, VIEW.14 needs GEN.417 and GEN.418, MAP.219 needs MAP.217, GEN.377 needs GEN.323, GEN.301 needs TEST.196.

TODO items resolved: none / Added: none
