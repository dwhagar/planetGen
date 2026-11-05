**Headline:** The course planner can now ask which keep-out spheres a straight path touches.

### Added
- `db/obstacles.py` lists the keep-out spheres a straight segment of the galaxy touches, in order along it, each with its radius, basis, and whether an end of the segment lies inside it. Systems, black holes, neutron stars, rogue planets, quasars, and asteroid fields count; the search for each kind uses that kind's own widest radius. Nothing shows it yet; the path planner (NAV.26) is its caller.

generation-output: unchanged
