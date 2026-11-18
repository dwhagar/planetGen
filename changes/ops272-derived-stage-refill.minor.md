**Headline:** Bringing a galaxy to the running epoch now also refills the stored values worked out from other stored values (habitability scores, the effect ledger) when the code that wrote them has changed, instead of leaving them old or regenerating sectors.

### Added
- **Bring-to-epoch refills derived values.** A new derived stage, after the population stage, finds systems whose habitability score or effect ledger was written by older code (the existing stamps, GEN.1479) and refills them with the existing refill jobs. No sector is generated again and `sectors.modified_at` is not touched. The survey and the dry run report it.
