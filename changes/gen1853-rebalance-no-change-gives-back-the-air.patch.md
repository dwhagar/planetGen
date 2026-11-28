**Headline:** A world whose carbon sink asks no change of its air now keeps exactly the air its doses were made on

### Fixed

- GEN.1853 (TEST.830): `biogenic_oxygen.rebalance_co2` stored the air again, solved again, when the sink asked no change, and the re-solved air ended a little apart (parts in a million: a few pascals of pressure, the trace gases) from the air it was given, while it answered "the air did not change". Since GEN.1821 (Train 58) the burial calls it for every world with life whose oxygen has not risen, so the doses, the ozone column and the ultraviolet index of those worlds (`life._rebuild` is only run when the air changed) were stale against their stored air; the radiation sweep named three of 858 bodies (classes E and G). A no-change answer now gives back every field the trials touched, so the air stands exactly as the doses were made on it. Stored output moves slightly for the worlds this reached.
