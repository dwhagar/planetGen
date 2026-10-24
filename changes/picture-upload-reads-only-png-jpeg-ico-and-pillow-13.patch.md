**Headline:** The picture upload now reads only PNG, JPEG and ICO files, and Pillow moves to 13.0.0, which fixes 37 reported flaws in its image decoders.

### Security
- **Pillow 12.3.0 to 13.0.0 in both lock files (SEC.34).** `pip-audit` listed 37 advisories for 12.3.0, all fixed in 13.0.0 (TIFF, GIF and PPM decoders reading out of bounds, a zero-width PPM, PDF appends). Only the pillow pin moved in `requirements.lock` and `requirements-server.lock`; the picture tests and the lock tests pass on 13.0.0.
- **The admin picture upload asks Pillow for PNG, JPEG and ICO only (SEC.35).** It used to let Pillow try every format it knows and reject all but three afterwards, so a TIFF, GIF or PPM header was parsed first. `Image.open` now gets `formats=` the three, so a file of any other kind is "That file is not a readable picture." and no other decoder runs. New cases in `test_web_seo.py` refuse a TIFF, a GIF and a PPM.
TODO items resolved: SEC.34, SEC.35. Added: none.
