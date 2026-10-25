**Headline:** An account can have a picture, and it may be a GIF, an animated GIF too, checked under frame and pixel caps and stored as a clean copy.

### Added
- **Your picture (SEC.36).** The account page has a "Your picture" section to upload or remove a PNG, JPEG, ICO or GIF, and the picture shows beside "Signed in as". `POST /api/auth/avatar` and `DELETE /api/auth/avatar` do the same for the API, and `GET /api/auth/me` returns the picture's address. Control schema v20 adds `admin_users.avatar` (the update adds the column).
- **The avatar rules.** At most 1 MiB, 1024 by 1024 pixels a frame, 100 frames, and 16 million pixels over all the frames. A GIF or PNG with several frames stays animated with its durations and loop count. The upload is refused in a plain sentence when it is over a cap, is a bomb header, or is another format.

### Security
- **A GIF is read only for an account picture, and never from the uploaded bytes.** Pillow is told which four formats to try, the size is read from the header before anything is decoded, a decompression-bomb warning is an error, `Image.MAX_IMAGE_PIXELS` is lowered to 25 million, and frames are unpacked one at a time under the caps. What is stored is Pillow's new encoding of the pixels, so comments, XMP, ICC profiles, EXIF and application extensions other than the loop block do not survive. It is served with its own content type and `nosniff`. The site logo, icons and preview picture keep the SEC.35 rule and still refuse a GIF.

TODO items resolved: SEC.36. Added: none.
