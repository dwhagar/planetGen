**Headline:** The TODO list records a failing galaxy-grid fuzz test as a bug to fix.

### Added
- **TEST.192, a float-cancellation failure in the galaxy-grid fuzz test (bug).** The Settings lane found that `test_skeleton_bound_is_an_upper_bound_at_every_angle` fails when the halo density is far above a tiny bound; the item carries the reproducer and goes to Core 1's bug phase, goal issue #1694.

TODO items resolved: none / Added: TEST.192
