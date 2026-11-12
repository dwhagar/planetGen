**Headline:** the message for a body that stores a quantity that is not a number no longer shows an item code.

### Changed

- `class_fit.classify` (GEN.1593) raised a `ClassError` whose text ended in the item code, which a reader sees; the code is gone from the message (it stays in the comment and the test), and `test_no_todo_codes_shown` passes again.
