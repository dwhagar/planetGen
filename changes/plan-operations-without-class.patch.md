**Headline:** the admin plan operations no longer call the class-editing functions the derived-class work removed, so moves, additions and regenerations work again

### Changed

A body's class has been derived from its conditions since epoch 121, but the admin plan operations still asked the stabiliser to pin a moved planet's class, set a class by hand, keep a class through a regeneration, and pick a class for a new moon. Every move, addition and regeneration of a planet failed with a TypeError or an AttributeError (14 plan tests). The operations now stabilise without a pin, a moved or added planet takes the class its new conditions give, a new moon is drawn and takes its own, and the class is no longer an aspect a plan can set or a regeneration can keep (a plan that names it is refused as unknown). No generated output changes.
