**Headline:** A new cloud session can now have MySQL 8, Apache, Redis and every test dependency ready from one pasted setup script.

### Added

- `examples/cloud-environment/setup-planetgen-env.sh`: a setup script for a Claude cloud environment that installs MySQL 8, Apache 2 with mod_wsgi, Redis, every Python dependency and the NLTK corpus, so a new session can run the test suite at once (`planetgen-env-up` starts the services). The README beside it says where to paste it; `docs/testing.md` links to both.

TODO items resolved: none / Added: none
