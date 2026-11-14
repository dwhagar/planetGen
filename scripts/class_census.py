#!/usr/bin/env python3
"""
The class census as a repeatable script (TEST.719): draw a seeded sample of systems, store what the planet classifier saw in each body, and print the tables of
reports/class-census-2026-11-11.md in one markdown format so runs of different trains compare line by line.

    python scripts/class_census.py collect OUT.jsonl SET [--workers 4] [--forced-young] [--caps-off all|d34] [--class-m TIER] [--count N]
    python scripts/class_census.py report DIR [--caps-off-dir DIR]

`collect` writes one JSON line per system (the stars, the belts, one record per planet and moon with its stored values, the classifier's quantities and the
conditions of the habitable-family classes M, R, P, O and H). `SET` is one of `SETS`: `natural`, `hab`, `stars`, `young`, `youngforced`, `capsoff`, `sunlike`. A run
that stopped continues where it ended (a system already in OUT.jsonl is skipped). It takes the host's memory and load gate (PERF.254, as `scripts/test.py`
does) for its workers before it draws.

`sunlike` is the Sun-like sample of the Class M meter (`scripts/class_m_meter.py`): `--count` systems (default 400) of the star types F7V to K2V in rotation,
seeds from 30000. `--class-m TIER` makes every system under that Class M abundance tier (`SystemConfig.ABUNDANCE_TIER`); the records carry it as `class_m`.
While it draws, `collect` prints the running counts (systems done of total, candidates, Class M, Blue).

`--forced-young` draws every random primary star with an age from 0.5 to 50 Myr, log-uniform, so the disk-era stars can be counted (the GEN.1553 disk-by-age row);
it changes the draw of this process only. `--caps-off` switches the pressure caps of the atmosphere recipe row d6 off in this process (`all` drops every cap,
`d34` only the tenuous 1 kPa and moderate 10 kPa caps): a measurement of what the caps hold, never a setting.

`report` reads every `*.jsonl` and `*.jsonl.gz` in DIR and prints the tables A to J (K with `--caps-off-dir`, the same seeds drawn with the caps off). It needs no library beyond
numpy. The same files give the same tables.
"""

import argparse
import bisect
import collections
import contextlib
import glob
import gzip
import json
import logging
import math
import multiprocessing
import os
import pathlib
import re
import sys
import time
import traceback

import numpy

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

SETS = ("natural", "hab", "stars", "young", "youngforced", "capsoff", "sunlike")
"""tuple: The sample names `collect` takes (`jobs_for`)."""

SOLAR_LUMINOSITY_W = 3.828e26
EARTH_MASS_KG = 5.972e24
AGE_CAP_GY = 13.79
"""float: The age of the galaxy a star's age is held to in the stars table."""

FAMILY = ("M", "R", "P", "O", "H")
"""tuple: The classes of the habitable family whose conditions every body record carries."""

BODY_ATTRS = """name body_type planet_type zone distance eccentricity mass radius gravity density period rotation_period_hours axial_tilt_deg
surface_temperature atmospheric_pressure hydrosphere ocean_fraction land_fraction water_mass_fraction ice_mass_fraction bond_albedo greenhouse_factor flux_variation
energy_flux_w_m2 internal_heat_flux_w_m2 tidal_flux_w_m2 envelope_mass_fraction core_mass_earth metal_mass_fraction composition_origin moon_origin parent_kind
planet_class class_match runner_up_class runner_up_match class_variant life_class life_stage life_microbial life_flora life_fauna life_sentient life_complexity
temp_hottest_k temp_coldest_k temp_livable_share surface_dose_msv_yr born_beyond_snow_line formation_distance_au hill_radius ocean_depth_km ocean_ph
relative_humidity cloud_cover tau_cloud tau_haze tau_dust groundwater aquifer_depth_km ground_type lid_regime life_system_age_gyr tier_air tier_pressure tier_temperature tier_water tier_radiation tier
xuv_flux_earth uv_surface_index sky_suns days_per_year solar_day_hours water_activity porosity crust_thickness_km core_h_fraction
phi4 mantle_redox ring_kind ring_mass_kg magnetic_moment_a_m2 dipole_class""".split()
SUN_TYPES = ("F7V", "F9V", "G0V", "G2V", "G4V", "G6V", "G8V", "K0V", "K1V", "K2V")
"""tuple: The star types of the `sunlike` sample, in rotation (the set of reports/class-m-blue-levers-2026-11-14.md)."""

SUNLIKE_SEED = 30000
SUNLIKE_COUNT = 400

STAR_ATTRS = """type yerkes_class mass luminosity radius temperature age lifespan metallicity_fe_h snow_line_au disk_kind disk_inner_au disk_outer_au disk_mass_solar a_crit_au
age_ceiling_gy initial_mass_sol""".split()


# ---------------------------------------------------------------------------------------------------------------- collect

def apply_caps_off(mode):
    """Switches the pressure caps of the atmosphere recipe off in this process (a measurement; nothing in the repository changes)."""
    from planetgen.physics import surface_air
    original = surface_air._retention_cap_pa
    if mode == "all":
        surface_air._retention_cap_pa = lambda mixture, subject: None
        return

    def without_tenuous_and_moderate(mixture, subject):
        cap = original(mixture, subject)
        return None if cap is not None and (abs(cap - 1.0e3) < 1.0 or abs(cap - 1.0e4) < 1.0) else cap
    surface_air._retention_cap_pa = without_tenuous_and_moderate


class ForcedYoung:
    """Every random primary star gets an age of 0.5 to 50 Myr, log-uniform, drawn from a generator the system's seed resets (GEN.1546 removed the age raise for planets, so nothing sends a star with temperate bodies back to billions of years); this process only."""

    def __init__(self):
        import random
        self.random = random.Random(0)

    def install(self):
        from planetgen.generation import star
        from planetgen.physics import stellar_evolution

        def young(*args, **kwargs):
            return 10 ** self.random.uniform(math.log10(5e-4), math.log10(0.05))
        star.sample_star_age_gy = young
        stellar_evolution.sample_star_age_gy = young

    def seed(self, seed):
        self.random.seed(seed)


CLASS_M = None
"""str: The Class M abundance tier every system is drawn under (`--class-m`), or None; set before the pool forks."""

FORCED = None
"""ForcedYoung: Set in the main process before the pool forks, so every worker holds the patched functions."""


def jsonable(value):
    if isinstance(value, bool) or value is None or isinstance(value, (str, int)):
        return value
    if isinstance(value, float):
        return value if math.isfinite(value) else repr(value)
    if isinstance(value, (list, tuple)):
        return [jsonable(item) for item in value]
    if isinstance(value, dict):
        return {str(key): jsonable(item) for key, item in value.items()}
    return repr(value)[:80]


class Capture(logging.Handler):
    def __init__(self):
        super().__init__(logging.WARNING)
        self.records = []

    def emit(self, record):
        self.records.append(f"{record.levelname}:{record.name}:{record.getMessage()[:200]}")


def body_record(body, is_moon, host_index, parent_name):
    from planetgen.physics import atmosphere, class_fit, habitability
    record = {"is_moon": is_moon, "host": host_index, "parent": parent_name}
    try:
        star = getattr(body, "star", None)
        record["host_type_str"] = getattr(star, "type", None)
        record["host_proxy"] = type(star).__name__
        record["host_hz"] = jsonable(getattr(body, "habitable_zone", None))
    except Exception:  # noqa: BLE001 -- a record without its host is still a record
        pass
    for attribute in BODY_ATTRS:
        try:
            record[attribute] = jsonable(getattr(body, attribute, None))
        except Exception as exc:  # noqa: BLE001
            record[attribute] = f"ERR {type(exc).__name__}"
    try:
        record["gas_pa"] = jsonable(atmosphere.gas_partials_pa(body))
    except Exception as exc:  # noqa: BLE001
        record["gas_pa"] = f"ERR {type(exc).__name__}"
    try:
        phi4 = record.get("phi4")
        record["overall_tier"] = habitability.tier(phi4) if isinstance(phi4, (int, float)) and not isinstance(phi4, bool) else None
    except Exception as exc:  # noqa: BLE001
        record["overall_tier"] = f"ERR {type(exc).__name__}"
    try:
        found = class_fit._found(body)
        record["found"] = jsonable(dict(found))
        verdict = class_fit.classify(body)
        record["fit"] = [verdict.letter, verdict.runner_up, verdict.match, verdict.runner_up_match, verdict.digit]
        record["scores"] = jsonable({key: round(value, 4) for key, value in class_fit.scores(body).items()})
        conditions = {}
        for row in class_fit.CLASSES:
            if row.letter in FAMILY and class_fit.passes_gates(row, found):
                conditions[row.letter] = [[c["quantity"], "min" if "min" in c else ("lt" if "lt" in c else "eq"), c.get("min", c.get("lt", c.get("equals"))),
                                           round(class_fit.membership(c, found), 3)] for c in row.when]
        record["family_conditions"] = conditions
    except Exception as exc:  # noqa: BLE001
        record["fit_error"] = f"{type(exc).__name__}: {exc}"[:200]
    return record


def draw_one(job):
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem
    from planetgen.util import draw
    tag, seed, kwargs = job
    if FORCED is not None:
        FORCED.seed(seed)
    capture = Capture()
    root = logging.getLogger()
    root.addHandler(capture)
    started = time.time()
    out = {"tag": tag, "seed": seed, "kwargs": kwargs}
    try:
        config = SystemConfig()
        config.PLANETS = True
        for key, value in kwargs.items():
            setattr(config, key, value)
        if CLASS_M is not None:
            config.ABUNDANCE_TIER = CLASS_M
            out["class_m"] = CLASS_M
        with draw.bound(seed):
            system = StarSystem(system_config=config)
        stars = []
        for star in system.stars:
            row = {}
            for attribute in STAR_ATTRS:
                try:
                    row[attribute] = jsonable(getattr(star, attribute, None))
                except Exception as exc:  # noqa: BLE001
                    row[attribute] = f"ERR {type(exc).__name__}"
            try:
                row["hz"] = jsonable(star.habitable_zone)
            except Exception:  # noqa: BLE001
                row["hz"] = None
            stars.append(row)
        out["stars"] = stars
        out["binary_type"] = jsonable(system.binary_type)
        out["wide_binary"] = jsonable(system.wide_binary)
        out["habitable_request"] = jsonable(getattr(system, "habitable_request", None))
        try:
            status = system.habitable_request_status
            out["habitable_request_status"] = jsonable(status() if callable(status) else status)
        except Exception:  # noqa: BLE001
            pass
        try:
            out["unmet"] = jsonable(system.unmet_requirements)
        except Exception:  # noqa: BLE001
            pass
        index = {id(star): number for number, star in enumerate(system.stars)}
        belts, bodies = [], []
        for planet in system.all_planets:
            if planet.body_type == "a":
                belts.append({key: jsonable(value) for key, value in vars(planet).items()
                              if isinstance(value, (int, float, str, bool)) or value is None})
                continue
            host = index.get(id(getattr(planet, "star", None)), -1)
            bodies.append(body_record(planet, False, host, None))
            for moon in planet.moons:
                bodies.append(body_record(moon, True, host, getattr(planet, "name", None)))
        out["belts"] = belts
        out["bodies"] = bodies
        out["comet_count"] = len(getattr(system, "comets", []) or [])
    except Exception:  # noqa: BLE001 -- the census counts a failed system, it does not stop on one
        out["error"] = traceback.format_exc()[-1500:]
    out["warnings"] = capture.records[:20]
    out["seconds"] = round(time.time() - started, 1)
    root.removeHandler(capture)
    return out


def jobs_for(name, count=None):
    """The `(tag, seed, config overrides)` of a sample; the seeds are fixed so two builds draw the same systems."""
    jobs = []
    if name == "natural":
        for base, count in ((1000, 300), (2000, 150), (3000, 150)):
            jobs += [("natural", base + i, {}) for i in range(count)]
    elif name == "hab":
        for star, base in (("G2V", 4000), ("K3V", 5000), ("F5V", 6000), ("M2V", 7000)):
            count = {"G2V": 60, "K3V": 50}.get(star, 25)
            jobs += [(f"hab-{star}", base + i, {"STAR_TYPE": star, "HABITABLE_WORLD": True}) for i in range(count)]
    elif name == "stars":
        for star, base in (("F5V", 8000), ("G2V", 8100), ("K5V", 8200), ("M2V", 8300)):
            jobs += [(f"star-{star}", base + i, {"STAR_TYPE": star}) for i in range(30)]
    elif name == "youngforced":
        jobs += [("young-forced", 20000 + i, {}) for i in range(100)]
    elif name == "capsoff":
        jobs += [("off-natural", 1000 + i, {}) for i in range(150)]
        jobs += [("off-hab-G2V", 4000 + i, {"STAR_TYPE": "G2V", "HABITABLE_WORLD": True}) for i in range(40)]
    elif name == "sunlike":
        jobs += [("sunlike", SUNLIKE_SEED + i, {"STAR_TYPE": SUN_TYPES[i % len(SUN_TYPES)]}) for i in range(count or SUNLIKE_COUNT)]
    elif name == "young":
        jobs += [("young", 9000 + i, {"AGE": "young"}) for i in range(40)]
    else:
        raise SystemExit(f"class_census: no sample called {name!r}; one of {', '.join(SETS)}")
    return jobs


@contextlib.contextmanager
def host_gate(workers, per_worker_mb=1500):
    """The memory and load gate (PERF.254) for `workers` draw workers, as `scripts/test.py` takes it; without Redis the run goes on its worker count alone."""
    try:
        from planetgen.queue import memorygate, redisqueue
        conn = redisqueue.connect()
        conn.ping()
    except Exception as exc:  # noqa: BLE001
        print(f"class_census: no memory gate ({type(exc).__name__}); the draw runs on its worker count alone", flush=True)
        yield
        return
    claim = memorygate.Claim("class-census", float(workers * per_worker_mb), float(workers), memorygate.LANE_OTHER)

    def waiting(decision, waited):
        if int(waited) % 60 < memorygate.POLL_SECONDS:
            print(f"class_census: waits at the memory gate: {decision.reason}", flush=True)

    with memorygate.hold(conn, claim, on_wait=waiting):
        yield


def done_before(path):
    done = set()
    if os.path.exists(path):
        with open(path) as handle:
            for line in handle:
                try:
                    record = json.loads(line)
                    done.add((record["tag"], record["seed"]))
                except (ValueError, KeyError):
                    pass
    return done


class Progress:
    """The running counts of a draw: systems done of total, rocky candidates (`abundance_tier.candidate`'s test on a stored record), Class M and Blue worlds."""

    def __init__(self, total):
        self.total, self.done, self.errors, self.candidates, self.class_m, self.blue = total, 0, 0, 0, 0, 0
        self.started = time.time()

    def add(self, record):
        self.done += 1
        self.errors += bool(record.get("error"))
        for body in record.get("bodies", []):
            self.candidates += body.get("body_type") == "t" and not body.get("is_moon") and body.get("zone") == "e"
            self.class_m += body.get("planet_class") == "M"
            self.blue += body.get("overall_tier") == "Blue"

    def line(self):
        elapsed = time.time() - self.started
        eta = f", about {elapsed / self.done * (self.total - self.done) / 60:.0f} min left" if self.done else ""
        return (f"{self.done} of {self.total} systems ({100.0 * self.done / max(1, self.total):.0f}%), {self.candidates} candidates, "
                f"{self.class_m} Class M, {self.blue} Blue, {self.errors} errors, {elapsed / 60:.1f} min{eta}")


def draw_sample(out, set_name, workers, count=None, every=25):
    """Draws the systems of `set_name` not yet in `out` and appends one line each, printing the `Progress` line every `every` systems (and at the end)."""
    jobs = jobs_for(set_name, count)
    done = done_before(out)
    jobs = [job for job in jobs if (job[0], job[1]) not in done]
    print(f"{set_name}: {len(jobs)} systems to do ({len(done)} already)", flush=True)
    progress = Progress(len(jobs))
    with host_gate(workers):
        with multiprocessing.get_context("fork").Pool(workers, maxtasksperchild=20) as pool, open(out, "a") as handle:
            for record in pool.imap_unordered(draw_one, jobs):
                handle.write(json.dumps(record) + "\n")
                handle.flush()
                progress.add(record)
                if progress.done % every == 0:
                    print(progress.line(), flush=True)
    print("finished:", progress.line(), flush=True)


def collect(args):
    global FORCED, CLASS_M
    if args.caps_off:
        apply_caps_off(args.caps_off)
    if args.forced_young:
        FORCED = ForcedYoung()
        FORCED.install()
    CLASS_M = args.class_m
    draw_sample(args.out, args.set, args.workers, args.count)


# ---------------------------------------------------------------------------------------------------------------- report

def number(value):
    """A float, or None for a missing or text value."""
    if isinstance(value, bool):
        return float(value)
    if isinstance(value, (int, float)) and math.isfinite(value):
        return float(value)
    return None


def spectral(text):
    """The spectral letter of a star's type text; `WD` for a white dwarf."""
    if not isinstance(text, str):
        return "?"
    found = re.match(r"([OBAFGKMLTY])\d?", text)
    if text.startswith("D") or "white dwarf" in text.lower():
        return "WD"
    return found.group(1) if found else text[:3]


class Census:
    """The records of a directory of `collect` files: `systems`, `bodies` (with the host's star values joined and the classifier's `found` quantities as `q_...`),
    `belts` and `stars`."""

    def __init__(self, directory):
        self.systems, self.bodies, self.belts, self.stars = [], [], [], []
        for path in sorted(glob.glob(os.path.join(directory, "*.jsonl")) + glob.glob(os.path.join(directory, "*.jsonl.gz"))):
            with (gzip.open(path, "rt") if path.endswith(".gz") else open(path)) as handle:
                for line in handle:
                    try:
                        self._add(json.loads(line))
                    except ValueError:
                        continue

    def _add(self, record):
        tag, seed = record["tag"], record["seed"]
        stars = record.get("stars", [])
        self.systems.append({"tag": tag, "seed": seed, "error": record.get("error"), "n_belts": len(record.get("belts", []))})
        hosts = {}
        for index, star in enumerate(stars):
            row = dict(star, tag=tag, seed=seed, idx=index, sp=spectral(star.get("type")))
            self.stars.append(row)
            hosts[index] = row
        for belt in record.get("belts", []):
            self.belts.append(dict(belt, tag=tag, seed=seed))
        for body in record.get("bodies", []):
            row = {"tag": tag, "seed": seed}
            for key, value in body.items():
                if key not in ("found", "family_conditions", "scores", "fit", "gas_pa", "life_living"):
                    row[key] = value
            for key, value in (body.get("found") or {}).items():
                row["q_" + key] = value
            fit = body.get("fit")
            if fit:
                row["fit_letter"], row["fit_runner"], row["fit_match"], row["fit_runner_match"], row["fit_digit"] = fit
            row["fam"] = body.get("family_conditions")
            host = hosts.get(body.get("host"))
            row["circumbinary"] = body.get("host") == -1
            row["host_type"] = host.get("type") if host else None
            row["host_lum"] = number(host.get("luminosity")) if host else None
            row["host_age"] = number(host.get("age")) if host else None
            row["sp"] = host["sp"] if host else self._proxy_class(body.get("host_type_str"))
            row["mass_e"] = number(row.get("q_mass_earth"))
            self.bodies.append(row)

    @staticmethod
    def _proxy_class(text):
        """A circumbinary body takes the first letter of its proxy's name (`Binary (K4V/M1V)`)."""
        found = re.search(r"\(([OBAFGKM])", text) if isinstance(text, str) else None
        return found.group(1) if found else "?"


def md(header, rows, floatfmt="{:.3g}"):
    """A markdown table: `header` the column names (the first is the row label), `rows` lists of cells."""
    if not rows:
        return "(none)"
    out = ["| " + " | ".join(str(h) for h in header) + " |", "|" + "|".join(["---"] * len(header)) + "|"]
    for row in rows:
        out.append("| " + " | ".join(floatfmt.format(v) if isinstance(v, float) else str(v) for v in row) + " |")
    return "\n".join(out)


def crosstab(items, row_key, col_key, percent=False):
    """Counts of `items` by a row value and a column value, as a markdown table (rows and columns sorted); `percent` gives each row's share."""
    counts = collections.defaultdict(collections.Counter)
    columns = set()
    for item in items:
        row, column = row_key(item), col_key(item)
        if column is None:
            continue                                  # a body with no class (a belt member) is in no column
        counts[row][column] += 1
        columns.add(column)
    columns = sorted(columns, key=str)
    rows = []
    for row in sorted(counts, key=str):
        total = sum(counts[row].values())
        rows.append([row] + [round(100.0 * counts[row][c] / total, 1) if percent else counts[row][c] for c in columns])
    return md([""] + [str(c) for c in columns], rows, "{:.1f}")


def heading(text):
    print("\n### " + text + "\n")


def value_counts(items, key):
    return collections.Counter(str(key(item)) for item in items)


def class_of(body):
    return body.get("planet_class") or "(asteroid)"


def sets_of(census):
    bodies = census.bodies
    return {"natural": [b for b in bodies if b["tag"] == "natural"],
            "habitable request": [b for b in bodies if b["tag"].startswith("hab")],
            "by star type": [b for b in bodies if b["tag"].startswith("star-")],
            "all": bodies}


FUNNEL_STAGES = (
    ("all planets and moons", lambda q: True),
    ("rocky (terrestrial by the class rule, not asteroid)", lambda q: q("body_type") == "t" and q("planet_class") is not None),
    ("air at 10 Pa or more", lambda q: q("pressure_pa") >= 10),
    ("gravity 0.75 to 2 g", lambda q: 0.75 <= q("gravity_g") < 2),
    ("mean temperature 258 to 343 K (liquid-water span)", lambda q: 258 <= q("temperature_k") < 343),
    ("air 50 to 250 kPa", lambda q: 5e4 <= q("pressure_pa") < 2.5e5),
    ("mean temperature 275 to 305 K (M window)", lambda q: 275 <= q("temperature_k") < 305),
    ("liquid water at the surface (phase_water = 1)", lambda q: q("phase_water") == 1),
    ("ocean cover 20 to 80 percent", lambda q: 0.2 <= q("ocean_fraction") < 0.8),
    ("oxygen 16 to 35.5 kPa", lambda q: 16000 <= q("o2_partial_pa") < 35500),
    ("carbon dioxide under 1 kPa", lambda q: q("co2_partial_pa") < 1000),
    ("reducing gases under 5 percent, H2 under 30, SO2 under 3", lambda q: q("reducing_share") < .05 and q("share_h2") < .3 and q("share_so2") < .03),
    ("bands: hottest < 310, coldest >= 250, night >= 245, wet bulb < 28",
     lambda q: q("hottest_region_k") < 310 and q("coldest_region_k") >= 250 and q("night_temperature_k") >= 245 and q("wet_bulb_c") < 28),
)
"""tuple: The funnel's stages, each a test of one body's classifier quantities (a missing quantity fails it) that every stage before it passed."""


def funnel(bodies):
    """How many bodies pass each stage of `FUNNEL_STAGES` and all before it."""
    counts = [0] * len(FUNNEL_STAGES)
    for body in bodies:
        def q(name, body=body):
            if name == "planet_class":
                return body.get("planet_class")
            if name == "body_type":
                return body.get("q_body_type")
            value = number(body.get("q_" + name))
            return math.nan if value is None else value

        for number_, (_, test) in enumerate(FUNNEL_STAGES):
            try:
                passed = bool(test(q))
            except TypeError:
                passed = False
            if not passed:
                break
            counts[number_] += 1
    return counts


def hard_failures(body, letter="M"):
    """The conditions of class `letter` the body fails outright (membership 0), as text."""
    fam = body.get("fam")
    if not isinstance(fam, dict) or letter not in fam:
        return None
    return [c[0] + {"min": ">=", "lt": "<", "eq": "="}[c[1]] + f"{c[2]:g}" for c in fam[letter] if c[3] <= 0.0]


def cut(values, edges):
    """The index of the half-open interval (edges[i], edges[i+1]] holding `values`, or None."""
    if values is None or not math.isfinite(values) or values <= edges[0] or values > edges[-1]:
        return None
    return bisect.bisect_left(edges, values) - 1


def percentile(values, share):
    return float(numpy.percentile(values, share))


def label(edges, index):
    return f"({edges[index]:g}, {edges[index + 1]:g}]"


def table_a(census):
    heading("A. Samples")
    rows = []
    for tag in sorted({s["tag"] for s in census.systems}):
        systems = [s for s in census.systems if s["tag"] == tag]
        bodies = [b for b in census.bodies if b["tag"] == tag]
        rows.append([tag, len(systems), sum(not b.get("is_moon") for b in bodies), sum(bool(b.get("is_moon")) for b in bodies),
                     sum(s["n_belts"] for s in systems), sum(bool(s["error"]) for s in systems)])
    print(md(["set", "systems", "planets", "moons", "belts", "errors"], rows))
    print("total systems", len(census.systems), "bodies", len(census.bodies))


def table_b(census):
    heading("B. Class counts (all bodies)")
    groups = sets_of(census)
    counters = {name: value_counts(bodies, class_of) for name, bodies in groups.items() if name != "all"}
    classes = sorted(set().union(*counters.values()), key=lambda c: (-counters["natural"][c], c))
    total = sum(counters["natural"].values()) or 1
    rows = [[c] + [counters[name][c] for name in counters] + [round(100.0 * counters["natural"][c] / total, 2)] for c in classes]
    print(md(["class"] + list(counters) + ["natural %"], rows))


def table_c(census):
    groups = sets_of(census)
    natural, stars = groups["natural"], groups["by star type"]
    heading("C. Natural sample: class by host star class (rows: star class; planets and moons)")
    print(crosstab(natural, lambda b: b["sp"], lambda b: b.get("planet_class")))
    heading("C2. Class by zone, natural")
    print(crosstab(natural, lambda b: b.get("zone"), lambda b: b.get("planet_class")))
    heading("C3. By-star-type sets: class share (%) by star type")
    print(crosstab(stars, lambda b: b["tag"], lambda b: b.get("planet_class"), percent=True))
    heading("C4. By-star-type sets: class share (%) by zone")
    print(crosstab(stars, lambda b: f"{b['tag']} / {b.get('zone')}", lambda b: b.get("planet_class")))


def table_d(census):
    heading("D. M family")
    for letter in "MOHPRL":
        print(letter, sum(b.get("planet_class") == letter for b in census.bodies), "of", len(census.bodies))
    rows = [[b["tag"], b["seed"], b.get("name"), b.get("host_type_str"), b.get("zone"), b.get("planet_class"), b.get("class_match"), b.get("q_gravity_g"),
             b.get("q_temperature_k"), b.get("q_pressure_pa"), b.get("q_ocean_fraction"), b.get("q_o2_partial_pa"), b.get("q_co2_partial_pa")]
            for b in census.bodies if b.get("planet_class") in tuple("MOHPRL")]
    print(md(["tag", "seed", "name", "host", "zone", "class", "match", "g", "K", "Pa", "ocean", "O2 Pa", "CO2 Pa"], rows))


def table_e(census):
    groups = sets_of(census)
    heading("E. Funnel")
    columns = {name: funnel(bodies) for name, bodies in groups.items()}
    print(md(["stage"] + list(columns), [[stage] + [columns[name][i] for name in columns] for i, (stage, _) in enumerate(FUNNEL_STAGES)]))
    heading("E2. Funnel by star class (all sets, all bodies)")
    by_class = collections.defaultdict(list)
    for body in census.bodies:
        by_class[body["sp"]].append(body)
    columns = {k: funnel(v) for k, v in sorted(by_class.items()) if len(v) > 40}
    print(md(["stage"] + list(columns), [[stage] + [columns[name][i] for name in columns] for i, (stage, _) in enumerate(FUNNEL_STAGES)]))
    heading("E3. Conditions of M failed, bodies with a condition list (all sets)")
    every, only, count = collections.Counter(), collections.Counter(), 0
    for body in census.bodies:
        hard = hard_failures(body)
        if hard is None:
            continue
        count += 1
        every.update(hard)
        if len(hard) == 1:
            only[hard[0]] += 1
    print("bodies", count)
    print(md(["condition", "bodies failing", "the only failure"], [[c, n, only[c]] for c, n in every.most_common()]))
    heading("E4. Nearest misses")
    ranked = []
    for body in census.bodies:
        fam = body.get("fam")
        if isinstance(fam, dict) and "M" in fam:
            hard = [c for c in fam["M"] if c[3] <= 0.0]
            ranked.append((len(hard), -sum(c[3] for c in fam["M"]), body, hard))
    ranked.sort(key=lambda item: (item[0], item[1]))
    rows = []
    for count_, _, body, hard in ranked[:15]:
        fails = "; ".join(f"{c[0]} {'>=' if c[1] == 'min' else '<' if c[1] == 'lt' else '='} {c[2]:g} (has {body.get('q_' + c[0], '?')})" for c in hard)
        rows.append([body["tag"], body["seed"], body.get("name"), body.get("host_type_str") or body.get("host_type"), body.get("zone"),
                     round(number(body.get("distance")) or 0.0, 3), body.get("planet_class"), count_, fails])
    print(md(["tag", "seed", "name", "star", "zone", "au", "class", "hard", "fails"], rows))


def table_f(census):
    def mass_by_flux(hosts, edges, title):
        heading(title)
        groups = collections.defaultdict(list)
        for body in census.bodies:
            luminosity = body.get("host_lum")
            if (body.get("is_moon") or body.get("body_type") == "a" or body["sp"] not in hosts or not luminosity or body.get("planet_type") != "terrestrial"):
                continue
            distance, mass = number(body.get("distance")), number(body.get("q_mass_earth"))
            if distance is None or mass is None:
                continue
            index = cut(distance / math.sqrt(luminosity / SOLAR_LUMINOSITY_W), edges)
            if index is not None:
                groups[index].append(mass)
        rows = [[label(edges, i), len(v), round(percentile(v, 25), 2), round(percentile(v, 50), 2), round(percentile(v, 75), 2), round(max(v), 2)]
                for i, v in sorted(groups.items())]
        print(md(["Sun-equivalent au", "count", "25%", "50%", "75%", "max"], rows))

    mass_by_flux("FGK", [0, 0.1, 0.2, 0.4, 0.7, 1.0, 1.5, 2.5, 5, 10, 100], "F. Rocky planet mass by Sun-equivalent flux distance, FGK hosts")
    mass_by_flux("M", [0, 0.4, 0.7, 1.0, 1.5, 2.5, 5, 100], "F2. Same for M hosts (distance scaled the same way)")


def air_bodies(census):
    return [b for b in census.bodies if b.get("q_body_type") == "t" and (number(b.get("q_pressure_pa")) or 0.0) >= 10]


def table_g(census):
    heading("G. Pressure pile-ups")
    air = air_bodies(census)

    def pinned(body):
        log = math.log10(body["q_pressure_pa"])
        return abs(log - round(log)) < 5e-4, round(log)

    at = collections.Counter(round_ for is_pinned, round_ in map(pinned, air) if is_pinned)
    print("air bodies", len(air), " at 1 kPa:", at[3], " at 10 kPa:", at[4], " at 100 kPa:", at[5])
    print(dict(sorted(value_counts([b for b in air if pinned(b)[0]], class_of).items(), key=lambda item: (-item[1], item[0]))))
    temperate = [b for b in air if 258 <= (number(b.get("q_temperature_k")) or -1) < 343 and 0.75 <= (number(b.get("q_gravity_g")) or -1) < 2]
    print("temperate (258-343 K) rocky air bodies with 0.75-2 g:", len(temperate), "pinned at 1 or 10 kPa:",
          sum(1 for b in temperate if pinned(b)[0] and pinned(b)[1] in (3, 4)), "reach 50 kPa:", sum(b["q_pressure_pa"] >= 5e4 for b in temperate))
    edges = [0, 100, 1e3 * 1.001, 1e4 * 1.001, 5e4, 2.5e5, 1e7, 1e13]
    counts = collections.Counter(cut(b["q_pressure_pa"], edges) for b in air)
    print(md(["pressure (Pa)", "bodies"], [[label(edges, i), counts[i]] for i in range(len(edges) - 1)], "{:g}"))


def table_h(census):
    heading("H. Oxygen chain on temperate rocky air bodies (258 to 343 K, air 10 Pa or more)")
    warm = [b for b in air_bodies(census) if 258 <= (number(b.get("q_temperature_k")) or -1) < 343]
    photosynthetic = ("photosynthesis", "complex_cells", "multicellularity")
    sea = [b for b in warm if b.get("q_phase_water") == 1 and b.get("life_stage") in photosynthetic]
    print(len(warm), "bodies;", "with liquid sea", sum(b.get("q_phase_water") == 1 for b in warm), "; life stage:", dict(value_counts(warm, lambda b: b.get("life_stage") or "none")))
    print("o2 >= 16 kPa:", sum((number(b.get("q_o2_partial_pa")) or 0) >= 16000 for b in warm), "; o2 > 0:", sum((number(b.get("q_o2_partial_pa")) or 0) > 0 for b in warm),
          "; photosynthesis or later and sea:", len(sea))
    print(md(["tag", "seed", "name", "host", "class", "g", "K", "Pa", "ocean", "O2 Pa", "CO2 Pa", "life stage", "host age Gy"],
             [[b["tag"], b["seed"], b.get("name"), b.get("host_type_str"), b.get("planet_class"), b.get("q_gravity_g"), b.get("q_temperature_k"), b.get("q_pressure_pa"),
               b.get("q_ocean_fraction"), b.get("q_o2_partial_pa"), b.get("q_co2_partial_pa"), b.get("life_stage"), b.get("host_age")] for b in sea]))


def table_i(census):
    heading("I. Stars: age cap, disks, young")
    stars = census.stars
    print("stars", len(stars), f"age >= {AGE_CAP_GY} Gy:", sum((number(s.get("age")) or 0) >= AGE_CAP_GY for s in stars))
    rows = []
    for letter in sorted({s["sp"] for s in stars}):
        ages = [number(s.get("age")) for s in stars if s["sp"] == letter and number(s.get("age")) is not None]
        rows.append([letter, len(ages), round(percentile(ages, 50), 2) if ages else "", sum(a >= AGE_CAP_GY for a in ages)])
    print(md(["class", "n", "median age Gy", "at the cap"], rows))
    bodies = collections.Counter((b["tag"], b["seed"]) for b in census.bodies)
    young = [s for s in stars if (number(s.get("age")) or 1.0) < 0.1]
    edges = [0, 3, 10, 30, 100]
    rows = []
    for index in range(len(edges) - 1):
        found = [s for s in young if cut(s["age"] * 1000, edges) == index]
        rows.append([label(edges, index) + " Myr", len(found), sum(s.get("disk_kind") is not None for s in found), sum(bodies[(s["tag"], s["seed"])] > 0 for s in found)])
    print(md(["age", "n", "with disk", "with bodies"], rows))
    print(md(["tag", "seed", "type", "age Gy", "disk", "disk Msun", "disk outer au", "bodies"],
             [[s["tag"], s["seed"], s.get("type"), s.get("age"), s.get("disk_kind"), s.get("disk_mass_solar"), s.get("disk_outer_au"), bodies[(s["tag"], s["seed"])]]
              for s in young[:40]]))
    print("disk kinds:", dict(value_counts(stars, lambda s: s.get("disk_kind"))))


def table_j(census):
    heading("J. Other value checks")
    bodies = census.bodies
    temperature = lambda b: number(b.get("q_temperature_k"))
    print("surface T pinned 2000 K:", sum(temperature(b) == 2000 for b in bodies), "; 375 K:", sum(temperature(b) is not None and round(temperature(b)) == 375 for b in bodies),
          "; T<10 K:", sum(temperature(b) is not None and temperature(b) < 10 for b in bodies))
    print("bond albedo 0.4898:", sum(number(b.get("bond_albedo")) is not None and round(number(b["bond_albedo"]), 4) == 0.4898 for b in bodies), "of", len(bodies))
    print("ocean_fraction == 1 on ice hydrosphere:", sum(b.get("q_ocean_fraction") == 1 and b.get("q_hydrosphere") in ("ice", "ice-covered ocean") for b in bodies))
    belts = census.belts
    earth = [(number(b.get("mass_kg")) or 0.0) / EARTH_MASS_KG for b in belts]
    print("belts", len(belts), "median Me", round(percentile(earth, 50), 3) if earth else "", ">1 Me:", sum(e > 1 for e in earth), ">10 Me:", sum(e > 10 for e in earth),
          "largest body > Earth radius:", sum((number(b.get("largest_body_km")) or 0.0) > 6371 for b in belts))
    with_planets = {(b["tag"], b["seed"]) for b in bodies if not b.get("is_moon")}
    print("systems with no planets (excluding belts):", len(census.systems) - len(with_planets))
    print("circumbinary bodies:", sum(b["circumbinary"] for b in bodies))


def table_k(census, off):
    heading("K. Caps off, identical seeds")
    base = [b for b in census.bodies if (b["tag"] == "natural" and 1000 <= b["seed"] <= 1149) or (b["tag"] == "hab-G2V" and 4000 <= b["seed"] <= 4039)]
    rows = []
    for name, group in (("caps on", base), ("caps off", off.bodies)):
        air = [b for b in group if b.get("q_body_type") == "t" and (number(b.get("q_pressure_pa")) or 0) >= 10]
        temperate = [b for b in air if 258 <= (number(b.get("q_temperature_k")) or -1) < 343 and 0.75 <= (number(b.get("q_gravity_g")) or -1) < 2]
        rows.append([name, len(group), len(air), sum(b["q_pressure_pa"] >= 5e4 for b in air), len(temperate), sum(b["q_pressure_pa"] >= 5e4 for b in temperate),
                     " ".join(f"{letter}:{sum(b.get('planet_class') == letter for b in group)}" for letter in "MOHPRL")])
    print(md(["", "bodies", "air", "air >= 50 kPa", "temperate 0.75-2 g", "of those >= 50 kPa", "classes"], rows))
    print("funnel caps on:", dict(zip((s for s, _ in FUNNEL_STAGES), funnel(base))))
    print("funnel caps off:", dict(zip((s for s, _ in FUNNEL_STAGES), funnel(off.bodies))))


def report(args):
    census = Census(args.directory)
    for table in (table_a, table_b, table_c, table_d, table_e, table_f, table_g, table_h, table_i, table_j):
        table(census)
    if args.caps_off_dir:
        table_k(census, Census(args.caps_off_dir))


def _tiers():
    from planetgen.generation import abundance_tier
    return abundance_tier.TIERS


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)
    drawer = sub.add_parser("collect", help="draw a sample and store one JSON line per system")
    drawer.add_argument("out")
    drawer.add_argument("set", choices=SETS)
    drawer.add_argument("--workers", type=int, default=max(1, min(4, (os.cpu_count() or 1))))
    drawer.add_argument("--forced-young", action="store_true", help="every random primary star gets an age of 0.5 to 50 Myr")
    drawer.add_argument("--caps-off", choices=("all", "d34"), help="switch the pressure caps off in this process")
    drawer.add_argument("--class-m", choices=_tiers(), help="draw every system under this Class M abundance tier")
    drawer.add_argument("--count", type=int, help=f"systems of the sunlike sample (default {SUNLIKE_COUNT})")
    drawer.set_defaults(run=collect)
    reporter = sub.add_parser("report", help="print the tables of a directory of collected files")
    reporter.add_argument("directory")
    reporter.add_argument("--caps-off-dir", help="a directory drawn with --caps-off on the same seeds (table K)")
    reporter.set_defaults(run=report)
    args = parser.parse_args(argv)
    args.run(args)


if __name__ == "__main__":
    main()
