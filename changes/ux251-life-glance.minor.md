**Headline:** A planet or moon page now has a short Life paragraph in At a glance, and the Life section gives three words for each kind of life: how much, which way it is going and how secure it is.

### Added
- **Life paragraph in At a glance.** It gives the words for each kind of life that is there, for example "Microbes: abundant, stable, secure. Flora: occasional, declining, precarious.", then the life classes the conditions suit, best first. A world that suits a class but has no life yet says "Nothing lives here yet."; a world that could not have life shows none.
- **Table of the kinds of life in the Life section.** Microbes, flora, fauna and sentient life each get how much (with the score out of 9), the trend (growing, stable or declining as the star brightens), the security (secure, marginal, precarious or dying out) and, when it is known, how long the conditions last. Under it: the stage, animal grade, pace of evolution, age of the system and the age the stage was reached.
- **New API route.** `GET /api/planets/<id>/life-readings` and `GET /api/moons/<id>/life-readings` give the same words; they are worked out on request and stored nowhere.

### Changed
- **The summary no longer lists life scores.** The sentence "Its life scores are ..." is gone; the At a glance paragraph says it in words.
- **The lines that repeated the Life paragraph are off the page.** "Class 1a life lives on this planet, evolving at a standard pace" and "Its latest milestone is ..." are now in the table under the Life paragraph. The written description of a system (the text and Markdown forms) still has them, now as their own paragraphs.
- **Stage sentences name the whole class.** "Class 3a life has grown into colonies ..." instead of "Class 3 life ...".
