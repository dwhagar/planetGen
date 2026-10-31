**Headline:** The galaxy's cluster field is now built ahead by the warm-up, in small saved parts, and the web only reads it.

### Changed
- The galaxy's cluster field (its globular clusters, satellites and young open clusters as one field the density model reads) took 20 to 60 seconds to build in whatever process first needed it, nearly all of it working out the mean system mass of each distinct age. It is now the first family of the warm-up's opening phase (`warm-cluster-field`, `galaxy/field_job.py`): one queued unit for every 1,000 ages works out its part and writes it as a file in the state folder, so every unit is a fraction of a second to a few seconds however busy the machine is; when the last part is in, the field is made from the rows and the parts and written as `galaxy-field/<signature>.field` by a rename. A warm-up that stopped carries on from the parts that are there. The signature names the stored rows and the source of the code that decides the numbers, so a changed row or formula is a new file, and the file is read through an unpickler that allows only the field's own classes.
- The web daemon only reads that file (`store._component_field`). A visitor's request for a galaxy of more than 1,500 components that has no file is refused at once with the 503 and `Retry-After` of the cache guard (cause `missing`) instead of building for a minute against the daemon's request limit; a worker, a command or the warm-up still builds it in its own process when no file is there, and a small galaxy is still built in the request (about a second).
- `ComponentField` takes the mean system masses by age as an optional argument, so the parts can hand them in; the numbers it gives are the same as before.

### Fixed
- A page that needed the cluster field after every generation run no longer waits a minute for it inside a web request.

generation-output: unchanged
TODO items resolved: PERF.225
