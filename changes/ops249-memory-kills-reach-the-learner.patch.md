**Headline:** The memory controller's learning now sees the job processes the kernel killed for memory (OPS.249).

### Fixed

- A job process that ended on SIGKILL, the kernel's kill for memory, is recorded as a `killed_memory` attempt (kind `killed`, as before). Since the old queue's unit path went (OPS.271), the launcher recorded it as `lost`, so `pressurelearn.read_record` found no kill and the controller's edges could never move down however many kills there were. A job that ends on any other signal stays `lost`. The retry on a grown claim, run alone, is as it was.
