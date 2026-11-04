**Headline:** The map scripts no longer carry the light-theme code that black map backgrounds left unused. Nothing changes on screen.

### Changed
- **UX.255.** Every map is drawn on black, so the branches that chose colors for a pale background are gone. `mapcore.js` loses `isLightBackground`; `nebulalook.js` keeps one set of nebula looks (`NEBULA_LOOKS`); `superbubble.js` keeps one bubble look (`LOOK`); the galaxy glow layer, galaxy labels, sector point shader, region blobs, sector and stage views no longer take a light-background option, and the `.map-labels[data-light]` rules are removed from `style.css`. The JS tests for these now check the black-background look only. No generator output changes.

generation-output: unchanged

TODO items resolved: UX.255.
