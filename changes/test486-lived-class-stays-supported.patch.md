**Headline:** A world's life class is always listed among the classes the world supports.

### Fixed
- **TEST.486: a cold world whose life used up its fuel kept its class but stored no supported classes.** When the air built by a class's own life (the hydrogen or donor it drew down) left the class scoring 0, the world kept the class and its stage while `life_class_scores` was empty and the page's "Supports" lines were missing. The class that lives on a world is now stored among the supported classes at the score it was chosen at, and counts for the best-supported class. Only those worlds change.

generation-output: changed
