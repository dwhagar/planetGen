# planetgen/api/schemas.py

"""
Request bodies as Pydantic models (ADM.21).

Each write route validates its JSON body by handing it to `parse_body`
with the model for that route. The models are strict (a string is not a
number, a bool is not an int), refuse fields they don't name, and report
every problem at once: a failed body is a 400 whose `error` joins them
and whose `errors` lists `{"field", "message"}` for each, so a form can
mark every field in one round trip. The limits are the ones the routes
always had (`routes.MAX_NAME_LENGTH` and the rest).
"""

import math
from typing import Annotated, List, Literal, Optional, Union

from pydantic import (
    BaseModel, ConfigDict, Field, StrictBool, StrictFloat, StrictInt, StrictStr, ValidationError, field_validator,
    model_validator,
)

from planetgen import tuning
from planetgen.admin import throttle
from planetgen.api.common import ApiError, is_http_url
from planetgen.db import store
from planetgen.generation import limits as generationLimits
from planetgen.population import facilities as facility_rules

MAX_NAME_LENGTH = 255
"""int: `sectors.name`/`star_systems.name` are VARCHAR(255)."""

MAX_WIKI_URL_LENGTH = 2048
"""int: `sectors.wiki_url` is VARCHAR(2048)."""

MAX_SECTOR_EDGE_LY = 1e9
"""float: An upper bound on a sector's `edge_ly` -- far beyond any real
sector, well short of overflowing the unit conversion."""

SYSTEM_NAME_MAX_LENGTH = store.SYSTEM_NAME_MAX_LENGTH
"""int: A system's or star's name leaves room for its planets' numerals."""

Number = Annotated[float, Field(strict=True, allow_inf_nan=False)]
"""A JSON number: an int or a float, not a bool, not NaN or infinity."""


class Body(BaseModel):
    """Base of every request model: strict types, unknown fields refused."""

    model_config = ConfigDict(extra="forbid", strict=True)


def _field_name(location):
    return ".".join(str(part) for part in location)


def errors_of(exc):
    """`[{"field", "message"}]` for a Pydantic `ValidationError`, with the
    unknown and missing fields each folded into one line."""
    unknown, missing, found = [], [], []
    for error in exc.errors():
        field = _field_name(error["loc"])
        if error["type"] == "extra_forbidden":
            unknown.append(field)
        elif error["type"] == "missing":
            missing.append(field)
        elif error["type"] == "value_error":
            found.append({"field": field, "message": str(error["ctx"]["error"])})
        else:
            found.append({"field": field, "message": f"'{field}' is invalid: {error['msg']}"
                          if field else error["msg"]})
    out = []
    if unknown:
        out.append({"field": ", ".join(sorted(unknown)), "message": f"unrecognized field(s): {', '.join(sorted(unknown))}"})
    if missing:
        out.append({"field": ", ".join(sorted(missing)), "message": f"missing required field(s): {', '.join(sorted(missing))}"})
    return out + found


def parse_body(model, body, message=None):
    """
    `body` checked against `model`. Returns the model instance. `message`
    replaces the joined text of the problems (the sign-in routes answer
    every bad body the same way); `errors` still lists them.

    Raises:
        ApiError: 400, whose message joins every problem and whose
            `errors` lists each one.
    """
    try:
        return model.model_validate(body)
    except ValidationError as exc:
        problems = errors_of(exc)
        raise ApiError(message or "; ".join(problem["message"] for problem in problems), errors=problems) from None


def given(instance):
    """The fields the caller sent, as a plain dict."""
    return instance.model_dump(exclude_unset=True)


def _name(value):
    if not value.strip():
        raise ValueError("'name' must be a non-empty string")
    if len(value) > MAX_NAME_LENGTH:
        raise ValueError(f"'name' must be at most {MAX_NAME_LENGTH} characters")
    return value


def _wiki_url(value):
    if value is not None and not (len(value) <= MAX_WIKI_URL_LENGTH and is_http_url(value)):
        raise ValueError("'wiki_url' must be an http or https URL or null")
    return value


class SectorCreate(Body):
    """`POST /api/sectors`."""

    name: StrictStr
    edge_ly: Number = Field(gt=0, le=MAX_SECTOR_EDGE_LY)

    _check_name = field_validator("name")(_name)


class SectorUpdate(Body):
    """`PATCH /api/sectors/<id>`: any non-empty subset."""

    name: StrictStr = None
    edge_ly: Number = Field(default=None, gt=0, le=MAX_SECTOR_EDGE_LY)
    wiki_url: Optional[StrictStr] = None

    _check_name = field_validator("name")(_name)
    _check_wiki = field_validator("wiki_url")(_wiki_url)

    @model_validator(mode="after")
    def _something(self):
        if not self.model_fields_set:
            raise ValueError("body must include at least one field to update")
        return self


class Rename(Body):
    """`PATCH /api/{stars,planets,moons}/<id>`: the new name, trimmed."""

    name: StrictStr

    @field_validator("name")
    @classmethod
    def _trim(cls, value):
        if not value.strip():
            raise ValueError("'name' must be a non-empty string")
        return " ".join(value.split())


class SystemRename(Rename):
    """A star's or system's new name: shorter, so its planets' names fit."""

    @field_validator("name")
    @classmethod
    def _short(cls, value):
        if len(value) > SYSTEM_NAME_MAX_LENGTH:
            raise ValueError(f"'name' must be at most {SYSTEM_NAME_MAX_LENGTH} characters")
        return value


class BodyRename(Rename):
    """A planet's or moon's new name."""

    @field_validator("name")
    @classmethod
    def _short(cls, value):
        if len(value) > MAX_NAME_LENGTH:
            raise ValueError(f"'name' must be at most {MAX_NAME_LENGTH} characters")
        return value


class SystemRecipe(Body):
    """The recipe fields `POST /api/systems` and a regeneration take
    (`SystemConfig` minus its per-orbit `slots`); `None` lets the
    generator decide."""

    markdown: StrictBool = None
    habitable_world: Optional[StrictBool] = None
    asteroid_belt: Optional[StrictBool] = None
    large_star: Optional[StrictBool] = None
    moons: Optional[StrictBool] = None
    max_planets: Optional[StrictBool] = None
    planets: Optional[StrictBool] = None
    intelligent_life: Optional[StrictBool] = None
    binary_system: Optional[StrictBool] = None
    wide_binary: Optional[StrictBool] = None
    star_type: Optional[StrictStr] = None
    name: Optional[StrictStr] = None
    age: Literal["young", "old", None] = None
    num_orbits: Optional[StrictInt] = None

    @field_validator("name")
    @classmethod
    def _recipe_name(cls, value):
        if value is not None and len(value) > SYSTEM_NAME_MAX_LENGTH:
            raise ValueError(f"'name' must be at most {SYSTEM_NAME_MAX_LENGTH} characters")
        return value

    @field_validator("num_orbits")
    @classmethod
    def _orbits(cls, value):
        if value is not None and not 0 < value <= generationLimits.MAX_NUM_ORBITS:
            raise ValueError(f"'num_orbits' must be a positive integer up to {generationLimits.MAX_NUM_ORBITS}, or null")
        return value


class SystemCreate(SystemRecipe):
    """`POST /api/systems`: a recipe, optionally placed in a sector."""

    sector_id: Optional[StrictStr] = Field(default=None, max_length=40)
    position: Optional[list[Number]] = None

    @field_validator("position")
    @classmethod
    def _three(cls, value):
        if value is not None and len(value) != 3:
            raise ValueError("'position' must be three finite numbers (light-years from the sector's center)")
        return value

    @model_validator(mode="after")
    def _needs_sector(self):
        if self.position is not None and self.sector_id is None:
            raise ValueError("'position' needs 'sector_id'")
        return self


class RegenerateRecipe(SystemRecipe):
    """A regeneration's recipe: it keeps the system's name."""

    @model_validator(mode="after")
    def _keeps_name(self):
        if "name" in self.model_fields_set:
            raise ValueError("'regenerate' keeps the system's name; rename with 'name' instead")
        return self


class SystemPatch(Body):
    """`PATCH /api/systems/<id>`: a rename, a regeneration, or both."""

    name: StrictStr = None
    regenerate: RegenerateRecipe = None
    drop_facilities: StrictBool = False

    @field_validator("name")
    @classmethod
    def _system_name(cls, value):
        return SystemRename.model_validate({"name": value}).name

    @model_validator(mode="after")
    def _something(self):
        if "name" not in self.model_fields_set and "regenerate" not in self.model_fields_set:
            raise ValueError("send 'name', 'regenerate', or both")
        return self


class NeighborhoodRequest(Body):
    """`POST /api/sectors/<id>/generate-neighborhood` (every field optional)."""

    radius_ly: Optional[Number] = Field(default=None, gt=0, le=generationLimits.MAX_GENERATE_RADIUS_LY)
    estimate_only: StrictBool = False


class FacilityCreate(Body):
    """`POST /api/facilities`."""

    name: StrictStr
    kind: StrictStr
    placement: StrictStr
    host_type: StrictStr
    host_id: StrictStr
    distance_km: Optional[Number] = Field(default=None, gt=0)
    phase_deg: Optional[Number] = None
    offset_ly: Optional[list[Number]] = None
    description: Optional[StrictStr] = Field(default=None, max_length=4000)

    _check_name = field_validator("name")(_name)

    @field_validator("kind")
    @classmethod
    def _kind(cls, value):
        if value not in tuning.FACILITY_KINDS:
            raise ValueError(f"'kind' must be one of: {', '.join(sorted(tuning.FACILITY_KINDS))}")
        return value

    @field_validator("placement")
    @classmethod
    def _placement(cls, value):
        if value not in facility_rules.PLACEMENTS:
            raise ValueError(f"'placement' must be one of: {', '.join(facility_rules.PLACEMENTS)}")
        return value

    @field_validator("host_type")
    @classmethod
    def _host_type(cls, value):
        if value not in facility_rules.HOST_TYPES:
            raise ValueError(f"'host_type' must be one of: {', '.join(facility_rules.HOST_TYPES)}")
        return value

    @field_validator("offset_ly")
    @classmethod
    def _offset(cls, value):
        if value is not None and len(value) != 3:
            raise ValueError("'offset_ly' must be three finite numbers [x, y, z] in light-years")
        return value


class WikiUpload(Body):
    """`POST /api/{systems,sectors}/<id>/wiki`."""

    backend: Literal["wikijs", "mediawiki"]
    path: Optional[StrictStr] = None

    @field_validator("path")
    @classmethod
    def _strip(cls, value):
        return value.strip() if value else None

    @model_validator(mode="after")
    def _wikijs_needs_a_path(self):
        if self.backend == "wikijs" and not self.path:
            raise ValueError("'path' is required for the 'wikijs' backend")
        return self


class DropFacilities(Body):
    """The optional body of a body delete or regeneration."""

    drop_facilities: StrictBool = False


class ChangeClass(Body):
    """`POST /api/{planets,moons,belts}/<id>/class`."""

    planet_class: StrictStr = Field(alias="class")
    force: StrictBool = False

    @field_validator("planet_class")
    @classmethod
    def _known(cls, value):
        if value not in tuning.PLANET_CLASSES:
            raise ValueError(f"'class' must be one of {', '.join(sorted(tuning.PLANET_CLASSES))}")
        return value


class ChangeStar(Body):
    """`POST /api/systems/<id>/star`."""

    star_type: StrictStr
    drop_facilities: StrictBool = False


class Lenient(BaseModel):
    """Base of the sign-in and admin-account models: strict types, but
    unknown fields are ignored as these routes always did."""

    model_config = ConfigDict(extra="ignore", strict=True)


MAX_API_KEY_LIFETIME_DAYS = 3650
"""int: The longest an API key may be made to last (API.9)."""

MAX_API_KEY_LABEL_LENGTH = 128
"""int: The longest label an API key takes."""


class Login(Lenient):
    """`POST /api/auth/login`."""

    username: StrictStr
    password: StrictStr

    @field_validator("username")
    @classmethod
    def _username(cls, value):
        value = value.strip()
        if not value:
            raise ValueError("username and password are required")
        return value

    @field_validator("password")
    @classmethod
    def _password(cls, value):
        if not value:
            raise ValueError("username and password are required")
        return value


class LoginCode(Lenient):
    """`POST /api/auth/login/code`: the pending token and a code."""

    pending: StrictStr
    code: StrictStr

    @field_validator("code")
    @classmethod
    def _code(cls, value):
        if not value.strip():
            raise ValueError("pending and code are required")
        return value


class ChangeCredentials(Lenient):
    """`POST /api/auth/change-credentials` (blank means missing)."""

    current_password: Optional[StrictStr] = None
    new_username: Optional[StrictStr] = None
    new_password: Optional[StrictStr] = None


class PasswordCheck(Lenient):
    """A two-factor change's re-check of the password."""

    current_password: Optional[StrictStr] = None


class TotpConfirm(Lenient):
    """`POST /api/auth/totp/confirm`."""

    code: Optional[StrictStr] = None


class TotpDisable(Lenient):
    """`POST /api/auth/totp/disable`."""

    current_password: Optional[StrictStr] = None
    code: Optional[StrictStr] = None


class ApiKeyCreate(Lenient):
    """`POST /api/auth/api-keys`."""

    label: StrictStr
    scopes: Optional[List[StrictStr]] = None
    expires_days: Optional[Union[StrictInt, StrictFloat]] = None

    @field_validator("label")
    @classmethod
    def _label(cls, value):
        value = value.strip()
        if not value:
            raise ValueError("'label' is required")
        if len(value) > MAX_API_KEY_LABEL_LENGTH:
            raise ValueError(f"'label' must be at most {MAX_API_KEY_LABEL_LENGTH} characters")
        return value

    @field_validator("scopes")
    @classmethod
    def _scopes(cls, value):
        if value is None:
            return None
        from planetgen.admin import auth as adminAuth
        try:
            return list(adminAuth.check_scopes(value))
        except ValueError as exc:
            raise ValueError(f"'scopes': {exc}")

    @field_validator("expires_days")
    @classmethod
    def _expires(cls, value):
        if value is not None and not (0 < value <= MAX_API_KEY_LIFETIME_DAYS):
            raise ValueError(f"'expires_days' must be more than 0 and at most {MAX_API_KEY_LIFETIME_DAYS}")
        return value


class NamingKeyChange(Body):
    """`POST /api/admin/naming-key`: a key, or `draw`."""

    key: Optional[StrictStr] = None
    draw: StrictBool = False

    @model_validator(mode="after")
    def _one(self):
        from planetgen.names import naming_key

        if self.draw:
            if self.key is not None:
                raise ValueError("send 'key' or 'draw', not both")
            return self
        if self.key is None:
            raise ValueError("send 'key' (8 hex digits) or 'draw': true")
        self.key = naming_key.parse_key(self.key)
        return self


class LockoutLift(Lenient):
    """`POST /api/admin/lockouts/lift`: one lockout, or `all`."""

    scope: Optional[StrictStr] = None
    subject: Optional[StrictStr] = None
    all: Optional[StrictBool] = None

    @model_validator(mode="after")
    def _one(self):
        if self.all is True:
            return self
        if self.scope not in throttle.SCOPES:
            raise ValueError("'scope' must be 'ip' or 'user'")
        if not self.subject or not self.subject.strip() or len(self.subject) > throttle.MAX_SUBJECT_LENGTH:
            raise ValueError("'subject' is required")
        return self


class WorkAction(Lenient):
    """`POST /api/admin/work/<id>/control`."""

    action: Optional[StrictStr] = None

    @model_validator(mode="after")
    def _known(self):
        from planetgen.queue import work

        if self.action not in work.ACTIONS:
            raise ValueError("'action' must be pause, resume or cancel")
        return self


class QueueAction(Lenient):
    """`POST /api/admin/work/queue`."""

    action: Optional[StrictStr] = None

    @model_validator(mode="after")
    def _known(self):
        if self.action not in ("pause", "resume"):
            raise ValueError("'action' must be pause or resume")
        return self
