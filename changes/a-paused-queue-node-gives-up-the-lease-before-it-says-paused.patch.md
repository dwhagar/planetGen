**Headline:** A paused work queue node now lets go of the lease before it says it is paused.

### Fixed
- **A node on the admin queue page could read "paused" while it still held the lease (TEST.288).** On a pause the work queue (`WorkQueue._standby`) first set its node to `paused` and only then gave up the lease, two separate writes; a reader that looked between them (the queue page, or a run waiting for the lease) saw a paused node that still blocked the other runs. Under load the gap is long enough that `test_a_paused_queue_node_stands_by_until_resumed[2]` failed on it, about once in two gates. The lease is now released first, and a new test with a slowed release pins the order.

TODO items resolved: TEST.288. Added: none.
