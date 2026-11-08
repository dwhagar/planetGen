import json, sys
tree = sys.argv[1]
sys.path.insert(0, tree + "/src"); sys.path.insert(0, tree + "/src/tests")
from multiprocessing import Pool
import class_sweep_cli
from planetgen.physics import atmosphere
def run(seed):
    rows = []
    for body, is_moon in class_sweep_cli.sample_bodies(100, seed, 1.0):
        p = atmosphere.gas_partials_pa(body)
        if not p: continue
        rows.append({"seed": seed, "moon": is_moon, "m": body.mass, "r": body.radius, "T": body.surface_temperature, "P": body.atmospheric_pressure,
                     "o2": p.get("O2", 0) / 1e3, "co2": p.get("CO2", 0) / 1e3, "ch4": p.get("CH4", 0) / 1e3, "cls": body.planet_class,
                     "hyd": getattr(body, "hydrosphere", None), "sink": getattr(body, "_co2_sink", None) is not None})
    return rows
if __name__ == "__main__":
    with Pool(4) as pool:
        out = [r for rs in pool.map(run, range(1, 9)) for r in rs]
    json.dump(out, open(sys.argv[2], "w")); print(len(out))
