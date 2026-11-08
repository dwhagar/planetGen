**Headline:** The System Map's 3D view keeps the info panel of the body you picked when the page is opened or reloaded on its address.

### Fixed
- UX.290: the diagram and the 3D view share one info panel (`#sysmap-info`), and the diagram sets itself up after the 3D view may already have picked the body named by `?object=`. When it finished last, it reset the panel to its "Click a star, planet..." hint, so a reload lost the panel (found as a failure of `test_system_page_3d_view_draws_switches_scale_and_keeps_the_diagram`, not as a fault of the time buttons: play, faster and now never touched the panel). The diagram now leaves the panel alone while it is hidden. A new browser test holds the diagram's last import back for 4 seconds so the order is the bad one every time.

generation-output: unchanged
