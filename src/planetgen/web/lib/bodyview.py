# planetgen/web/lib/bodyview.py

"""
The rendered example on an object's page (UX.171): the look record (`web/maps/bodylook.py`, MAP.271) of a star,
planet or moon, handed to `static/bodyview.js` as `data-look`, which paints the same texture the system map paints
and shades it as a lit sphere on a canvas.

A belt or a comet has no look record yet, so its page has no example. The caption says how the picture was worked
out (`look_note`, a stand-in until the colour account of UX.158 / MAP.271 replaces it).
"""

from planetgen.web.maps.bodylook import look_json, look_note, planet_look, star_look


def view(kind, row, stars):
    """`{"look", "note"}` for a star, planet or moon row of a system's detail (`stars` are the system's star rows, for the
    temperature of the star a body's air is seen under), or `None` for a kind with no picture."""
    if kind == "star":
        look = star_look(row)
    elif kind in ("planet", "moon"):
        planet = row["planet"] if kind == "moon" else row
        host = next((s for s in stars if s["id"] == planet.get("star_id")), stars[0] if stars else None)
        look = planet_look(dict(row, _star_temperature_k=host.get("temperature_k") if host else None), kind)
    else:
        return None
    return {"look": look_json(look), "note": look_note(look)}
