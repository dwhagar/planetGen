**Headline:** The plain-words check reads sixteen systems instead of ten so its corpus stays big enough, and two wording slips it found are fixed.

### Fixed
- The "Formerly" line on object pages no longer shows the epoch number; it says what changed, why and the date, in the story's words.
- A belt's sheet row reads "Largest body", not "Largest Body".

### Changed
- `test_plain_words` reads the Document of sixteen seeded systems. The ten it read had shrunk to 1,306 sentences and labels under the 1,500 floor, because the epochs of the last days redrew those systems (three now draw a lone star). The floor stays.
- Found with the extra systems and listed under its owner: the outer-disk belt origin sentence says "disc" and is stored, so fixing it rides a train.

generation-output: unchanged
TODO items resolved: TEST.576
