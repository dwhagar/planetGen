# tests/test_db_login_server.py

"""
OPS.49: `db.login.connect` against a small fake MySQL 8 server that has no TLS and wants caching_sha2_password's full
authentication (the first login after a server start): pymysql asks for the RSA public key, sends the encrypted
password, and the server answers OK. pymysql 1.2.1 to 1.2.3 drop that OK (their `caching_sha2_password_auth` lacks a
`return`) and crash; PyMySQL PR 1275 adds the `return`. The login must work with either, so both versions of the step
are built from the installed source and the login runs with each, without pinning a pymysql.
"""

import inspect
import re
import socket
import struct
import threading

import pymysql
import pytest
from cryptography.hazmat.primitives import serialization
from cryptography.hazmat.primitives.asymmetric import rsa
from pymysql import _auth

from planetgen.db import login

OK = b"\x00\x00\x00\x02\x00\x00\x00"


def _step_source(returns):
    """pymysql's `caching_sha2_password_auth` as installed, rewritten to end with (or without) `return pkt`."""
    source = inspect.getsource(_auth.caching_sha2_password_auth).rstrip()
    source = re.sub(r"\n\s+return pkt$", "", source)
    return source + ("\n    return pkt\n" if returns else "\n")


def _step(returns):
    namespace = {}
    exec(_step_source(returns), vars(_auth), namespace)
    return namespace["caching_sha2_password_auth"]


def _read(sock):
    header = b""
    while len(header) < 4:
        chunk = sock.recv(4 - len(header))
        if not chunk:
            raise EOFError
        header += chunk
    length, sequence = struct.unpack("<I", header[:3] + b"\0")[0], header[3]
    payload = b""
    while len(payload) < length:
        payload += sock.recv(length - len(payload))
    return sequence, payload


def _write(sock, sequence, payload):
    sock.sendall(struct.pack("<I", len(payload))[:3] + bytes([sequence]) + payload)


def _serve(listener, public_pem):
    connection, _ = listener.accept()
    with connection:
        salt = b"abcdefgh" + b"ijklmnopqrst"
        _write(connection, 0, b"\x0a" + b"8.0.46\0" + struct.pack("<I", 7) + salt[:8] + b"\0"
               + struct.pack("<HBHHB", 0xA201, 45, 2, 0x38, 21) + b"\0" * 10 + salt[8:] + b"\0"
               + b"caching_sha2_password\0")
        _read(connection)                               # the client's handshake response
        _write(connection, 2, b"\x01\x04")              # the password isn't cached: full authentication
        assert _read(connection)[1] == b"\x02"          # the client asks for the RSA public key
        _write(connection, 4, b"\x01" + public_pem)
        assert len(_read(connection)[1]) == 256         # the password, RSA encrypted
        _write(connection, 6, OK)
        while True:
            try:
                sequence, payload = _read(connection)
            except EOFError:
                return
            if payload[:1] == b"\x01":                  # COM_QUIT
                return
            _write(connection, sequence + 1, OK)


@pytest.fixture
def server():
    key = rsa.generate_private_key(public_exponent=65537, key_size=2048)
    pem = key.public_key().public_bytes(serialization.Encoding.PEM, serialization.PublicFormat.SubjectPublicKeyInfo)
    listener = socket.socket()
    listener.bind(("127.0.0.1", 0))
    listener.listen(1)
    thread = threading.Thread(target=_serve, args=(listener, pem), daemon=True)
    thread.start()
    yield listener.getsockname()[1]
    listener.close()
    thread.join(5)


def _login(port):
    return login.connect(host="127.0.0.1", port=port, user="planetgen", password="secret", ssl_disabled=True,
                         connect_timeout=5)


def test_the_unfixed_step_alone_cannot_finish_this_login(server, monkeypatch):
    """Guards the test: without the hook the step that returns nothing ends in the NoneType error."""
    monkeypatch.setattr(_auth, "caching_sha2_password_auth", _step(returns=False))
    with pytest.raises(AttributeError, match="is_auth_switch_request"):
        pymysql.connect(host="127.0.0.1", port=server, user="planetgen", password="secret", ssl_disabled=True,
                        connect_timeout=5)


@pytest.mark.parametrize("returns", [False, True], ids=["pymysql without the return", "pymysql with the return"])
def test_the_login_works_whether_or_not_pymysql_returns_the_servers_reply(server, monkeypatch, returns):
    monkeypatch.setattr(_auth, "caching_sha2_password_auth", _step(returns))
    connection = _login(server)
    assert connection.open
    connection.close()


def test_the_login_works_with_the_installed_pymysql(server):
    connection = _login(server)
    assert connection.open
    connection.close()
