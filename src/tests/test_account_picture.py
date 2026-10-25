# tests/test_account_picture.py

"""
An account's picture (SEC.36): a GIF, animated or not, a PNG, a JPEG or an ICO that is checked under caps, then stored
as a new encoding of its pixels with nothing else the upload carried. The site's own pictures (logo, icons, preview)
keep the SEC.35 rule and still refuse a GIF.
"""

import base64
import io
import os
import re
import struct
import sys
import zlib

import pytest
from PIL import Image, PngImagePlugin

from planetgen.db import store
from planetgen.util import settings
from planetgen.web import assets, csrf, seo, seo_pages

from tests.test_web_seo import _configure
from tests.test_web_settings import config, site  # noqa: F401 -- fixtures


@pytest.fixture(autouse=True)
def _assets_folder(tmp_path, monkeypatch):
    monkeypatch.setenv("PLANETGEN_ASSETS_DIR", str(tmp_path / "assets"))
    settings.reset_cache()
    yield
    settings.reset_cache()
    seo.clear_cache()
    seo_pages.clear_caches()


def _gif(frames=1, size=(16, 16), comment=b"a secret comment", xmp=False, loop=3, duration=50):
    images = [Image.new("RGBA", size, (255 * (i % 2), 40 * i % 255, 90, 255)) for i in range(frames)]
    buffer = io.BytesIO()
    options = {"save_all": True, "append_images": images[1:], "duration": [duration * (i + 1) for i in range(frames)]}
    if frames > 1:
        options["loop"] = loop
    images[0].save(buffer, "GIF", comment=comment, **options)
    data = buffer.getvalue()
    if xmp:
        flags = data[10]
        offset = 13 + (3 * 2 ** ((flags & 7) + 1) if flags & 0x80 else 0)
        packet = b"<x:xmpmeta/>"
        block = b"\x21\xff\x0bXMP DataXMP" + bytes([len(packet)]) + packet + b"\x00"
        data = data[:offset] + block + data[offset:]
    return data


def _png(size=(16, 16), text=None):
    buffer = io.BytesIO()
    info = None
    if text:
        info = PngImagePlugin.PngInfo()
        info.add_text("Comment", text)
    Image.new("RGB", size, (10, 20, 30)).save(buffer, "PNG", pnginfo=info)
    return buffer.getvalue()


def _set(site, data):
    return site.post("/api/auth/avatar", json={"data": base64.b64encode(data).decode()})


def _stored(url):
    return assets.read(url.rsplit("/", 1)[1])


def _frames(data):
    with Image.open(io.BytesIO(data)) as image:
        found = []
        for index in range(getattr(image, "n_frames", 1)):
            image.seek(index)
            found.append(image.info.get("duration"))
        return found, image.info.get("loop")


def _accounts(site):
    conn = store.get_control_connection(site.config)
    try:
        return conn.execute("SELECT id, username, avatar FROM admin_users ORDER BY id").fetchall()
    finally:
        conn.close()


def test_a_plain_gif_is_accepted_and_served_as_a_gif(site, config):
    _configure(config)
    response = _set(site, _gif())
    assert response.status_code == 200, response.get_data(as_text=True)
    body = response.get_json()
    assert re.fullmatch(r"/site-assets/[0-9a-f]{16}\.gif", body["avatar"])
    assert (body["width"], body["height"], body["frames"]) == (16, 16, 1)
    served = site.get(body["avatar"])
    assert served.mimetype == "image/gif" and served.headers["X-Content-Type-Options"] == "nosniff"
    assert "immutable" in served.headers["Cache-Control"]
    assert served.get_data().startswith(b"GIF8")
    assert site.get("/api/auth/me").get_json()["avatar"] == body["avatar"]


def test_an_animated_gif_stays_animated_with_its_timing_and_loop(site, config):
    _configure(config)
    body = _set(site, _gif(frames=5, loop=4)).get_json()
    assert body["frames"] == 5
    durations, loop = _frames(site.get(body["avatar"]).get_data())
    assert durations == [50, 100, 150, 200, 250] and loop == 4


def test_what_is_stored_has_no_comment_and_no_xmp(site, config):
    _configure(config)
    upload = _gif(frames=2, xmp=True)
    assert b"a secret comment" in upload and b"XMP DataXMP" in upload
    stored = site.get(_set(site, upload).get_json()["avatar"]).get_data()
    assert b"a secret comment" not in stored and b"XMP Data" not in stored and stored != upload
    assert b"NETSCAPE2.0" in stored            # the loop block is the one extension kept


def test_a_png_text_chunk_a_jpeg_exif_tag_and_an_ico_are_cleaned_too(site, config):
    _configure(config)
    stored = site.get(_set(site, _png(text="hidden words")).get_json()["avatar"]).get_data()
    assert b"hidden words" not in stored and stored.startswith(b"\x89PNG")
    exif = Image.Exif()
    exif[0x010E] = "hidden description"
    jpeg = io.BytesIO()
    Image.new("RGB", (16, 16), (200, 0, 0)).save(jpeg, "JPEG", exif=exif)
    assert b"hidden description" in jpeg.getvalue()
    stored = site.get(_set(site, jpeg.getvalue()).get_json()["avatar"]).get_data()
    assert b"hidden description" not in stored and stored.startswith(b"\xff\xd8")
    ico = io.BytesIO()
    Image.new("RGBA", (32, 32), (0, 200, 0, 255)).save(ico, "ICO", sizes=[(16, 16), (32, 32)])
    body = _set(site, ico.getvalue()).get_json()
    assert body["avatar"].endswith(".png") and (body["width"], body["height"]) == (32, 32)


def _with_screen(data, width, height):
    """`data` (a GIF) with the size its header states changed, the pixels untouched."""
    return data[:6] + struct.pack("<HH", width, height) + data[10:]


def _png_claiming(width, height):
    def chunk(kind, body):
        return struct.pack(">I", len(body)) + kind + body + struct.pack(">I", zlib.crc32(kind + body) & 0xFFFFFFFF)
    return (b"\x89PNG\r\n\x1a\n" + chunk(b"IHDR", struct.pack(">IIBBBBB", width, height, 8, 0, 0, 0, 0))
            + chunk(b"IDAT", zlib.compress(b"\0" * 10)) + chunk(b"IEND", b""))


def _other(name):
    buffer = io.BytesIO()
    Image.new("RGB", (8, 8)).save(buffer, name)
    return buffer.getvalue()


@pytest.mark.parametrize("data, why", [
    pytest.param(b"", "empty", id="empty"),
    pytest.param(b"not a picture", "not a readable picture", id="text"),
    pytest.param(_other("TIFF"), "not a readable picture", id="tiff"),
    pytest.param(_other("PPM"), "not a readable picture", id="ppm"),
    pytest.param(_other("BMP"), "not a readable picture", id="bmp"),
    pytest.param(_gif(frames=101, size=(8, 8)), "more than 100 frames", id="too-many-frames"),
    pytest.param(_gif(size=(1025, 1025)), "at most 1,048,576 pixels", id="too-many-pixels"),
    pytest.param(_with_screen(_gif(), 5000, 5000), "at most 1,048,576 pixels", id="header-over-the-pixel-cap"),
    pytest.param(_with_screen(_gif(), 6000, 6000), "decompression bomb", id="header-in-the-warning-range"),
    pytest.param(_with_screen(_gif(), 60000, 60000), "decompression bomb", id="gif-bomb-header"),
    pytest.param(_png_claiming(60000, 60000), "decompression bomb", id="png-bomb-header"),
    pytest.param(_gif(frames=70, size=(512, 512)), "too large to unpack", id="too-many-pixels-over-all-frames"),
    pytest.param(_gif(frames=2, size=(8, 8)) + b"\0" * 1_100_000, "most allowed", id="too-many-bytes"),
])
def test_a_picture_over_a_cap_or_of_another_format_is_refused_in_a_sentence(site, config, data, why):
    _configure(config)
    response = _set(site, data)
    assert response.status_code == 422, response.get_data(as_text=True)
    assert why in response.get_json()["error"]
    assert not os.path.exists(assets.assets_dir()) or not os.listdir(assets.assets_dir())
    assert site.get("/api/auth/me").get_json()["avatar"] is None


def test_the_caps_are_tested_before_the_frames_are_decoded(monkeypatch):
    data = _gif(frames=101, size=(8, 8))
    unpacked = []
    real = Image.Image.convert

    def counting(self, *args, **kwargs):
        if sys._getframe(1).f_code.co_filename == assets.__file__:     # the checker's own unpacking, not Pillow's
            unpacked.append(1)
        return real(self, *args, **kwargs)

    monkeypatch.setattr(Image.Image, "convert", counting)
    with pytest.raises(assets.AssetError):
        assets.prepare_avatar(data)
    assert len(unpacked) == assets.AVATAR_MAX_FRAMES            # the 101st frame is never unpacked


@pytest.mark.parametrize("kind", ["favicon_png", "apple_touch_icon", "default_image"])
def test_a_gif_is_still_refused_for_the_site_pictures(site, config, kind):
    _configure(config)
    response = site.post(f"/api/admin/images/{kind}", json={"data": base64.b64encode(_gif(size=(1200, 1200))).decode()})
    assert response.status_code == 422
    assert "not a readable picture" in response.get_json()["error"]


def test_a_new_picture_replaces_the_old_file_and_removing_deletes_it(site, config):
    _configure(config)
    first = _set(site, _png((16, 16))).get_json()["avatar"]
    second = _set(site, _png((17, 17))).get_json()["avatar"]
    assert first != second and site.get(first).status_code == 404 and site.get(second).status_code == 200
    assert _set(site, _png((17, 17))).get_json()["avatar"] == second          # the same picture, the same file
    assert site.delete("/api/auth/avatar").get_json() == {"avatar": None}
    assert site.get(second).status_code == 404 and site.get("/api/auth/me").get_json()["avatar"] is None


def test_a_file_two_accounts_share_stays_until_the_last_one_drops_it(site, config):
    _configure(config)
    url = _set(site, _png()).get_json()["avatar"]
    name = url.rsplit("/", 1)[1]
    conn = store.get_control_connection(site.config)
    try:
        conn.execute("INSERT INTO admin_users (username, password_hash, must_change_credentials, avatar) "
                     "VALUES ('second', 'x', 0, ?)", (name,))
        conn.commit()
    finally:
        conn.close()
    assert site.delete("/api/auth/avatar").status_code == 200
    assert site.get(url).status_code == 200                                    # 'second' still has it
    conn = store.get_control_connection(site.config)
    try:
        from planetgen.admin import auth as adminAuth
        adminAuth.set_avatar(conn, [a["id"] for a in _accounts(site) if a["username"] == "second"][0], None)
        adminAuth.forget_avatar_file(conn, name)
    finally:
        conn.close()
    assert site.get(url).status_code == 404


def test_a_picture_needs_a_signed_in_browser_session(site, config):
    _configure(config)
    anonymous = site.application.test_client()
    assert anonymous.post("/api/auth/avatar", json={"data": base64.b64encode(_png()).decode()}).status_code == 401
    assert anonymous.delete("/api/auth/avatar").status_code == 401
    assert site.post("/api/auth/avatar", json={"data": "***"}).status_code == 422
    assert site.post("/api/auth/avatar", json={"data": "AAAA", "extra": 1}).status_code == 400


def test_the_account_page_uploads_shows_and_removes_the_picture(site, config):
    _configure(config)
    page = site.get("/account").get_data(as_text=True)
    assert "Your picture" in page and 'class="avatar"' not in page
    response = site.post("/account/avatar", data={
        csrf.FIELD_NAME: site.csrf, "file": (io.BytesIO(_gif(frames=3)), "me.gif")}, content_type="multipart/form-data")
    assert response.status_code == 200
    page = response.get_data(as_text=True)
    assert "Your picture is saved." in page and re.search(r'<img class="avatar" src="/site-assets/[0-9a-f]{16}\.gif"', page)
    assert "Remove my picture" in page
    page = site.post("/account/avatar", data={csrf.FIELD_NAME: site.csrf, "remove": "1"}).get_data(as_text=True)
    assert "Your picture is removed." in page and 'class="avatar"' not in page


def test_the_account_page_says_why_a_picture_was_refused(site, config):
    _configure(config)
    response = site.post("/account/avatar", data={
        csrf.FIELD_NAME: site.csrf, "file": (io.BytesIO(_gif(frames=101, size=(8, 8))), "big.gif")},
        content_type="multipart/form-data")
    assert response.status_code == 422 and "more than 100 frames" in response.get_data(as_text=True)
    response = site.post("/account/avatar", data={csrf.FIELD_NAME: site.csrf}, content_type="multipart/form-data")
    assert response.status_code == 400 and "Choose a file to upload." in response.get_data(as_text=True)


def test_the_picture_is_not_offered_before_the_first_password_is_changed(site, config):
    _configure(config)
    conn = store.get_control_connection(site.config)
    try:
        conn.execute("UPDATE admin_users SET must_change_credentials = 1")
        conn.commit()
    finally:
        conn.close()
    assert "Your picture" not in site.get("/account").get_data(as_text=True)
    response = site.post("/account/avatar", data={csrf.FIELD_NAME: site.csrf, "remove": "1"})
    assert response.status_code in (302, 303)
