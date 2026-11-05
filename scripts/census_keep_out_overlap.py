#!/usr/bin/env python3
"""
Counts how many keep-out spheres overlap in a generated galaxy (NAV.56). The answer decides whether the path planner's cluster
fallback (NAV.26, docs/design/orbits/course-avoidance.md section 5) ever fires: if no two spheres overlap, a course never has to go
around two at once and the fallback is skipped; if some do, the fallback is built.

It reads the galaxy's stored positions and takes each object's keep-out radius the way `db/query.keep_out_radius` does (NAV.54):

* a star system: its galactic Hill radius, the larger of the stars' and the binary's `system_perimeter_km`;
* a black hole, a neutron star or a rogue planet that was placed in the galaxy: its galactic Hill radius against the galaxy's mass
  inside its distance from the center (read from a table of the enclosed mass, 64 radii, log-spaced), never smaller than the object;
* a supermassive black hole and a quasar: the sphere of influence (and a quasar's light); an asteroid field: its stored radius.

Planets, moons, belts and the radiation radii inside a system are left out: they sit inside their star's sphere, and a course that
leaves or enters a system is exempt from that star's sphere by containment.

A pair overlaps when the distance of the centers is under the sum of the radii; it is nested when one center lies inside the other
sphere. A cluster is a set of spheres joined by overlaps. The report counts the spheres, the pairs (by the kinds of the two), the
spheres in at least one overlap, the nested pairs and the clusters by size.

    python scripts/census_keep_out_overlap.py --database planetgen [--json out.json]

It only reads the database.
"""

import argparse
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

KM_PER_PC = 3.085677581491367e13
ENCLOSED_GRID_POINTS = 64
"""int: Radii at which the galaxy's enclosed mass is read; the rest is interpolated in log-log."""


def _parse(argv):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--database", required=True, help="the galaxy's database")
    parser.add_argument("--json", default=None, help="also write the numbers to this file")
    return parser.parse_args(argv)


def census(centers, radii, kinds):
    """
    The overlap census of spheres.

    Args:
        centers (array): `(n, 3)` centers, parsecs.
        radii (array): `(n,)` radii, parsecs.
        kinds (list of str): The kind of each sphere.

    Returns:
        dict: `spheres`; `by_kind` (spheres of each kind); `pairs` (overlapping pairs); `pairs_by_kinds` (`"a+b"` sorted); `nested`
            (pairs where one center lies inside the other sphere); `in_overlap` (spheres in at least one pair); `clusters` (a
            `{size: count}` for the sets of two or more spheres that overlap one another, directly or through others);
            `largest_cluster`.
    """
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    from scipy.spatial import cKDTree

    centers = np.asarray(centers, dtype=float).reshape(-1, 3)
    radii = np.asarray(radii, dtype=float)
    n = len(radii)
    by_kind = {}
    for kind in kinds:
        by_kind[kind] = by_kind.get(kind, 0) + 1
    report = {"spheres": n, "by_kind": dict(sorted(by_kind.items())), "pairs": 0, "pairs_by_kinds": {}, "nested": 0,
              "in_overlap": 0, "clusters": {}, "largest_cluster": 0}
    if n < 2:
        return report
    tree = cKDTree(centers)
    found = np.array(sorted(tree.query_pairs(2.0 * float(radii.max()), output_type="ndarray").tolist()), dtype=int).reshape(-1, 2)
    if not len(found):
        return report
    first, second = found[:, 0], found[:, 1]
    distance = np.linalg.norm(centers[first] - centers[second], axis=1)
    overlap = distance < radii[first] + radii[second]
    first, second, distance = first[overlap], second[overlap], distance[overlap]
    report["pairs"] = int(len(first))
    for a, b in zip(first.tolist(), second.tolist()):
        key = "+".join(sorted((kinds[a], kinds[b])))
        report["pairs_by_kinds"][key] = report["pairs_by_kinds"].get(key, 0) + 1
    report["pairs_by_kinds"] = dict(sorted(report["pairs_by_kinds"].items()))
    report["nested"] = int(np.count_nonzero((distance < radii[first]) | (distance < radii[second])))
    report["in_overlap"] = int(len(np.unique(np.concatenate([first, second]))))
    if len(first):
        graph = coo_matrix((np.ones(len(first)), (first, second)), shape=(n, n))
        _count, labels = connected_components(graph, directed=False)
        sizes = np.bincount(labels)
        for size in sizes[sizes > 1].tolist():
            report["clusters"][size] = report["clusters"].get(size, 0) + 1
        report["clusters"] = dict(sorted(report["clusters"].items()))
        report["largest_cluster"] = int(sizes.max())
    return report


def _enclosed_table(conn, radii_pc):
    """`(grid_pc, mass_kg)`: the galaxy's mass inside each of `radii_pc` (a log-spaced grid over their span)."""
    from planetgen.db import query
    low, high = max(float(np.min(radii_pc)), 1e-3), max(float(np.max(radii_pc)), 2e-3)
    grid = np.geomspace(low, high * 1.0001, ENCLOSED_GRID_POINTS)
    return grid, np.array([query._enclosed_mass_kg(conn, float(r)) for r in grid])


def _hill_pc(galactic_radius_pc, mass_kg, grid, enclosed):
    """The galactic Hill radius, parsecs: `a (m / 3 M)^(1/3)` with M the galaxy's mass inside `a` (log-log interpolation)."""
    central = np.exp(np.interp(np.log(np.maximum(galactic_radius_pc, 1e-3)), np.log(grid), np.log(enclosed)))
    return galactic_radius_pc * np.cbrt(mass_kg / (3.0 * central))


def load_spheres(conn):
    """`(centers, radii_pc, kinds)` of every placed object that has a keep-out sphere, as the module docstring says."""
    from planetgen.galaxy import keepout
    from planetgen.physics import constants

    centers, radii, kinds = [], [], []

    def add(kind, xyz, radius_pc):
        for point, radius in zip(xyz, radius_pc):
            if radius and radius > 0.0 and all(v is not None for v in point):
                centers.append(point)
                radii.append(float(radius))
                kinds.append(kind)

    rows = conn.execute(
        "SELECT s.galaxy_x_pc AS x, s.galaxy_y_pc AS y, s.galaxy_z_pc AS z, "
        "GREATEST(COALESCE(s.binary_system_perimeter_km, 0), COALESCE(MAX(st.system_perimeter_km), 0)) AS perimeter_km "
        "FROM star_systems s LEFT JOIN stars st ON st.star_system_id = s.id WHERE s.galaxy_x_pc IS NOT NULL "
        "GROUP BY s.id, s.galaxy_x_pc, s.galaxy_y_pc, s.galaxy_z_pc").fetchall()
    add("system", [(r["x"], r["y"], r["z"]) for r in rows], [(r["perimeter_km"] or 0.0) / KM_PER_PC for r in rows])

    sun = repr(constants.SOLAR_MASS_TO_KG)
    compact = {
        "black_hole": conn.execute(
            f"SELECT center_x_pc AS x, center_y_pc AS y, center_z_pc AS z, galactic_radius_pc AS r, mass_solar * {sun} AS mass_kg, "
            "event_horizon_radius_km AS size_km, mass_solar, mass_class FROM black_holes "
            "WHERE star_id IS NULL AND center_x_pc IS NOT NULL").fetchall(),
        "neutron_star": conn.execute(
            f"SELECT center_x_pc AS x, center_y_pc AS y, center_z_pc AS z, galactic_radius_pc AS r, mass_solar * {sun} AS mass_kg, "
            "radius_km AS size_km FROM neutron_stars WHERE star_id IS NULL AND center_x_pc IS NOT NULL").fetchall(),
        "rogue_planet": conn.execute(
            "SELECT center_x_pc AS x, center_y_pc AS y, center_z_pc AS z, galactic_radius_pc AS r, mass_kg, radius_km AS size_km "
            "FROM rogue_planets WHERE center_x_pc IS NOT NULL").fetchall(),
    }
    every_radius = [r["r"] for rows_ in compact.values() for r in rows_ if r["r"]]
    grid, enclosed = _enclosed_table(conn, every_radius or [1.0, 2.0])
    for kind, rows_ in compact.items():
        hole = [r for r in rows_ if kind == "black_hole" and r["mass_class"] == "supermassive"]
        rest = [r for r in rows_ if not any(r is h for h in hole)]
        add("supermassive_black_hole", [(r["x"], r["y"], r["z"]) for r in hole],
            [keepout.sphere_of_influence_km(r["mass_solar"]) / KM_PER_PC for r in hole])
        if rest:
            radius = np.array([r["r"] or 0.0 for r in rest], dtype=float)
            hill = _hill_pc(radius, np.array([r["mass_kg"] for r in rest], dtype=float), grid, enclosed)
            floor = np.array([r["size_km"] for r in rest], dtype=float) / KM_PER_PC
            add(kind, [(r["x"], r["y"], r["z"]) for r in rest], np.where(radius > 0.0, np.maximum(hill, floor), floor))
    rows = conn.execute("SELECT center_x_pc AS x, center_y_pc AS y, center_z_pc AS z, black_hole_mass_solar AS m, luminosity_w AS l "
                        "FROM quasars WHERE center_x_pc IS NOT NULL").fetchall()
    add("quasar", [(r["x"], r["y"], r["z"]) for r in rows],
        [keepout.quasar_keep_out(r["m"], r["l"]).radius_km / KM_PER_PC for r in rows])
    rows = conn.execute("SELECT center_x_pc AS x, center_y_pc AS y, center_z_pc AS z, radius_ly FROM asteroid_fields "
                        "WHERE center_x_pc IS NOT NULL").fetchall()
    add("asteroid_field", [(r["x"], r["y"], r["z"]) for r in rows], [keepout.asteroid_field_keep_out(r["radius_ly"]).radius_km / KM_PER_PC
                                                                    for r in rows])
    return np.array(centers, dtype=float).reshape(-1, 3), np.array(radii, dtype=float), kinds


def render(report):
    """The report as text."""
    lines = [f"{report['spheres']:,} keep-out spheres: " + ", ".join(f"{n:,} {k}" for k, n in report["by_kind"].items())]
    lines.append(f"{report['pairs']:,} overlapping pairs, {report['nested']:,} of them nested (one center inside the other sphere); "
                 f"{report['in_overlap']:,} spheres ({100.0 * report['in_overlap'] / max(report['spheres'], 1):.2f}%) are in at least one overlap")
    for key, count in report["pairs_by_kinds"].items():
        lines.append(f"  {key}: {count:,}")
    if report["clusters"]:
        lines.append("clusters of overlapping spheres, by size: " + ", ".join(f"{size}: {count:,}" for size, count in report["clusters"].items())
                     + f"; the largest holds {report['largest_cluster']:,}")
    else:
        lines.append("no two spheres overlap: the planner's cluster fallback never fires")
    return "\n".join(lines)


def main(argv=None):
    args = _parse(argv)
    from planetgen.db import store
    conn = store.get_connection(store.MySQLConfig(database=args.database), ensure_schema=False)
    try:
        centers, radii, kinds = load_spheres(conn)
    finally:
        conn.close()
    report = census(centers, radii, kinds)
    print(render(report))
    if args.json:
        with open(args.json, "w", encoding="utf-8") as handle:
            json.dump(report, handle, indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
