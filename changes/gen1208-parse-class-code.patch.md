**Headline:** A written class code such as F2 now works in the class filters and the class page address: a letter alone still matches every digit, a letter and a digit match that one position.

### Changed

- GEN.1208: `class_code.parse_code` is read where a class code is typed or linked: the planet and moon class filters of the database queries and of the search page (`F` matches F, F1 to F4 as before; `F2` matches digit 2 of the family only, `M0` the ideal), and the class page address (`/classes/planet/F2` goes to the page of class F, the digit being a position of the family; an unknown code is still a 404). No stored value changes, and nothing that was matched before is matched differently.

generation-output: unchanged
