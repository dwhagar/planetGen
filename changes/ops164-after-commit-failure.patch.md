**Headline:** A failed after-commit step is now logged with its name and its traceback.

### Changed
- **After-commit failures name the step and carry the traceback (OPS.164).** These steps keep derived data such as caches and counts in step with the database, so a failure used to leave them drifting with only a one-line message. generation-output: unchanged
