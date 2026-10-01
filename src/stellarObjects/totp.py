# stellarObjects/totp.py

"""
Time-based one-time codes (RFC 6238, the six-digit codes of any
authenticator app) for admin two-factor sign-in (SEC.26), in the standard
library only: HMAC-SHA1, 30-second steps, six digits -- the settings every
authenticator app supports.

    secret = new_secret()                    # base32, stored per admin
    uri = provisioning_uri(secret, "boss")   # for the QR code
    step = verify(secret, "123456", last_step=0)
    # an int (store it as the new last_step) or None

A code is accepted for the current step and one step either side (clock
drift), and never twice: `verify` refuses any step at or before
`last_step`, so a code seen over someone's shoulder can't be replayed.
"""

import base64
import hashlib
import hmac
import secrets
import struct
import time
from urllib.parse import quote

STEP_SECONDS = 30
DIGITS = 6
DRIFT_STEPS = 1
ISSUER = "planetGen"


def new_secret():
    """A random 160-bit secret, base32 without padding (32 characters)."""
    return base64.b32encode(secrets.token_bytes(20)).decode("ascii").rstrip("=")


def _key(secret):
    secret = secret.strip().replace(" ", "").upper()
    return base64.b32decode(secret + "=" * (-len(secret) % 8))


def code_at(secret, step):
    """The code for time step `step` (RFC 4226 HOTP with counter `step`)."""
    digest = hmac.new(_key(secret), struct.pack(">Q", step), hashlib.sha1).digest()
    offset = digest[-1] & 0x0F
    value = struct.unpack(">I", digest[offset:offset + 4])[0] & 0x7FFFFFFF
    return str(value % 10 ** DIGITS).zfill(DIGITS)


def current_step(now=None):
    return int((time.time() if now is None else now) // STEP_SECONDS)


def normalize_code(code):
    """`code` without spaces or dashes, or `None` if that isn't six digits."""
    if not isinstance(code, str):
        return None
    code = code.replace(" ", "").replace("-", "")
    return code if len(code) == DIGITS and code.isdigit() else None


def verify(secret, code, last_step=0, now=None):
    """
    The time step `code` matches (within `DRIFT_STEPS` of now and after
    `last_step`), or `None`.
    """
    code = normalize_code(code)
    if code is None or not secret:
        return None
    now_step = current_step(now)
    for step in range(now_step - DRIFT_STEPS, now_step + DRIFT_STEPS + 1):
        if step > (last_step or 0) and hmac.compare_digest(code_at(secret, step), code):
            return step
    return None


def provisioning_uri(secret, username, issuer=ISSUER):
    """The `otpauth://` URI an authenticator app reads from the QR code."""
    label = quote(f"{issuer}:{username}", safe="")
    return (f"otpauth://totp/{label}?secret={secret}&issuer={quote(issuer, safe='')}"
            f"&algorithm=SHA1&digits={DIGITS}&period={STEP_SECONDS}")


def qr_svg(text):
    """
    `text` as a QR code in SVG (dark modules on a white, quiet-zoned
    square), drawn with the vendored `qrcodegen` (Project Nayuki, MIT).
    """
    from . import qrcodegen
    qr = qrcodegen.QrCode.encode_text(text, qrcodegen.QrCode.Ecc.MEDIUM)
    border = 4
    size = qr.get_size() + border * 2
    parts = []
    for y in range(qr.get_size()):
        for x in range(qr.get_size()):
            if qr.get_module(x, y):
                parts.append(f"M{x + border},{y + border}h1v1h-1z")
    return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {size} {size}" width="{size * 5}" '
            f'height="{size * 5}" role="img" aria-label="QR code for an authenticator app" '
            f'shape-rendering="crispEdges"><rect width="100%" height="100%" fill="#fff"/>'
            f'<path d="{"".join(parts)}" fill="#000"/></svg>')
