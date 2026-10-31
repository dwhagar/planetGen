# planetgen/physics/life_air.py

"""
The gases a life class puts into the air besides oxygen (GEN.891, GEN.857, GEN.858, GEN.859)
================================================================================================

`physics/biogenic_oxygen.py` gives the air the oxygen of a class that splits water. This module gives it the rest of the class's ledger
(`physics/data/life_ledger.json`, docs/design/habitability/life-classes.md section 4): the methane that anoxygenic, hydrogen-eating and
(before the oxygen rise) oxygenic life makes, and the hydrogen the hydrogen-eaters draw down.

**Methane.** A class's flux is its fuel supply (a share of Earth's, `physics/life_supply.py`; for the oxygenic class the organic carbon its
photosynthesis buries, `biogenic_oxygen.productivity`) times the area over Earth's times the flux Earth's supply gives
(`earth_donor_methane_tmol_yr`: the food web of the anoxygenic photoferrotrophs, 10 Tmol a year, inside the 9 to 67 Tmol of the published
Archean model; `earth_h2_tmol_yr` over `h2_per_methane` for the hydrogen-eaters). The air holds the flux times the life of methane in an air with
no free oxygen (`methane_lifetime_yr`), at most `methane_max_share` of its molecules, and the carbon of that methane comes out of the air's
own carbon dioxide, of which (with the methane) `carbon_left_share` stays carbon dioxide at least: methane cannot be made from carbon the world does not have.

**Hydrogen and hydrogen sulfide.** An air where hydrogen-eating life lives holds no more than `hydrogen_max_share` of hydrogen, and one where
sulfur eaters or acid eaters live no more than `h2s_max_share` of hydrogen sulfide (GEN.894). What the life makes of the gas it takes is the
row's `products`: the acid eaters of Classes 5 and 9 turn the hydrogen sulfide into the same moles of sulfur dioxide, the hydrogenators of Class 7
turn two thirds of a mole of methane out of each mole of hydrogen (the acetylene they also take is not an air gas).

**Stage.** A class's gases begin when its producers do (the second step of the schedule) and a class whose methane ends at the oxygen rise
(`until: "oxygen"`) gives none to a world whose air has risen: the oxygen has burned it (`biogenic_oxygen.burn`). Before the producers, and
for a world that has no life, the air is the abiotic one.

What this does not do: the haze's cooling of the ground (the haze itself, its ultraviolet shield and its colour are redrawn from the new air
by `generation.life`), the sulfate, sulfur and silica the sulfur and acid classes leave in the water and rock (not air gases), and the hydrosphere: a world that the new air
moves across the freezing point keeps its stored state of water, as it does after the oxygen. Every number is a stand-in with its tier and
tweak level in the ledger file.
"""

import json
import os
from dataclasses import dataclass

from planetgen.physics import biogenic_oxygen, constants, elements, habitability_world, life_supply

DATA_FILE = os.path.join(elements.DATA_DIRECTORY, "life_ledger.json")
TMOL = 1.0e12
FUELS = ("burial", "donor_supply", "h2_supply")
"""tuple: What a methane row's fuel may be: the carbon buried (`burial`) or a supply of `physics/life_supply.py`."""
STATUSES = ("built", "pending")
DRAW_DOWN_CAPS = {"H2": "hydrogen_max_share", "H2S": "h2s_max_share"}
"""dict: Each gas a class may draw down and the ledger number that caps its share of the air's molecules where that life lives."""
HISTORY_MINIMUM_KPA = 1.0e-4
"""float: A gas the Life and air chart draws is above this partial pressure, kPa, at some age (0.1 Pa)."""


class LifeLedgerError(ValueError):
    """A ledger file that is refused."""


def _check_entry(owner, entry):
    if entry.get("tier") not in tuple(f"T{n}" for n in range(1, 10)) or entry.get("tweak") not in tuple(f"L{n}" for n in range(5)) or not entry.get("basis"):
        raise LifeLedgerError(f"{owner}: every number and class needs a tier T1..T9, a tweak L0..L4 and a basis, got {entry!r}")


def _load():
    with open(DATA_FILE, encoding="utf-8") as handle:
        document = json.load(handle)
    for name, entry in document["numbers"].items():
        _check_entry(name, entry)
        if "value" not in entry:
            raise LifeLedgerError(f"{name}: a number needs a value")
    for number, row in document["classes"].items():
        _check_entry(f"class {number}", row)
        if row.get("status") not in STATUSES:
            raise LifeLedgerError(f"class {number}: the status is one of {STATUSES}")
        methane = row.get("methane")
        if row["status"] == "built" and methane and methane.get("fuel") not in FUELS:
            raise LifeLedgerError(f"class {number}: the methane fuel is one of {FUELS}")
        for gas in row.get("draws_down", ()):
            if gas not in DRAW_DOWN_CAPS:
                raise LifeLedgerError(f"class {number}: a gas drawn down is one of {tuple(DRAW_DOWN_CAPS)}, got {gas!r}")
        for gas, product in (row.get("products") or {}).items():
            if gas not in row.get("draws_down", ()) or product.get("per") not in document["numbers"] or not product.get("gas"):
                raise LifeLedgerError(f"class {number}: a product is of a gas the class draws down, a gas and a ledger number per mole, got {gas!r}: {product!r}")
    return document


DOCUMENT = _load()
NUMBERS = DOCUMENT["numbers"]
CLASSES = DOCUMENT["classes"]
VERSION = DOCUMENT["version"]


@dataclass(frozen=True)
class Living:
    """
    A life class that lives on a body (GEN.1213): the class (`life_classes.LifeClass`), its `share` of the body's food production (0 to 1; the shares of the
    classes on one body add up to 1, and a body with one class gives it 1), and the age, Gyr, its producers began (`producers_gyr`; `None` reads the body's own
    `life_age_photosynthesis_gyr`, the headline class's second step).
    """

    life_class: object
    share: float = 1.0
    producers_gyr: float | None = None

    def started(self, body):
        """The age, Gyr, this class's producers began on `body`, or `None`."""
        return self.producers_gyr if self.producers_gyr is not None else getattr(body, "life_age_photosynthesis_gyr", None)


def number(name):
    """The value of the ledger's number `name`."""
    return NUMBERS[name]["value"]


def row_of(life_class):
    """The ledger row of a `life_classes.LifeClass` (by its class number)."""
    return CLASSES[str(life_class.number)]


def methane_flux_mol_yr(body, life_class, supplies=None):
    """
    The methane, mol a year, the class makes on `body` (0 for a class with no methane row, a pending row, a body without radius or a fuel it
    does not have). `supplies` is `life_supply.supplies` of the body, worked out if not given.
    """
    row = row_of(life_class)
    methane = row.get("methane")
    if row["status"] != "built" or not methane or not getattr(body, "radius", None) or not getattr(body, "gravity", None):
        return 0.0
    area_ratio = biogenic_oxygen.area_m2(body) / biogenic_oxygen.EARTH_AREA_M2
    fuel = methane["fuel"]
    if fuel == "burial":
        return number("methane_per_buried_carbon") * biogenic_oxygen.productivity(body, biogenic_oxygen.area_m2(body))
    if supplies is None:
        supplies = life_supply.supplies(habitability_world.values_from_body(body))
    share = supplies.get(fuel, 0.0)
    if fuel == "donor_supply":
        return share * area_ratio * number("earth_donor_methane_tmol_yr") * TMOL
    return share * area_ratio * number("earth_h2_tmol_yr") * TMOL / number("h2_per_methane")


def methane_target_mol(flux_mol_yr, moles):
    """The methane, mol, an air of `moles` (formula to mol) holds at that flux: the flux over its life, within `methane_max_share` of the
    air's molecules (making methane out of carbon dioxide keeps their number) and within the carbon the air has (`carbon_left_share` of its
    carbon dioxide and methane stays carbon dioxide)."""
    steady = flux_mol_yr * number("methane_lifetime_yr")
    share = number("methane_max_share")
    carbon = moles.get("CO2", 0.0) + moles.get("CH4", 0.0)
    return min(steady, share * sum(moles.values()), (1.0 - number("carbon_left_share")) * carbon)


def _oxygenic(body, living):
    """The `Living` of the oxygenic classes whose producers have begun by the star's age, as `(age_gyr, share)` pairs."""
    star_age = getattr(getattr(body, "star", None), "age", None)
    pairs = []
    for one in living:
        started = one.started(body)
        if one.life_class.releases_oxygen and started is not None and star_age is not None and started <= star_age:
            pairs.append((started, one.share))
    return pairs


def standing_methane_mol(body, living, age_gyr):
    """The methane, mol, the classes of `living` that have begun by `age_gyr` and whose methane the oxygen ends hold in the oxygen-free air of `body` (the
    steady amount of their summed flux, `methane_target_mol`); the oxygen of the rise has to burn it first (GEN.1213)."""
    moles = biogenic_oxygen.air_moles(body)
    if not moles:
        return 0.0
    supplies = None
    flux = 0.0
    for one in living:
        row = row_of(one.life_class)
        started = one.started(body)
        if row["status"] != "built" or not row.get("methane") or row["methane"].get("until") != "oxygen" or started is None or started > age_gyr:
            continue
        if row["methane"]["fuel"] != "burial" and supplies is None:
            supplies = life_supply.supplies(habitability_world.values_from_body(body))
        flux += one.share * methane_flux_mol_yr(body, one.life_class, supplies)
    return methane_target_mol(flux, moles) if flux > 0.0 else 0.0


def oxygen_state(body, living):
    """
    The oxygen the classes of `living` have made on `body` by the star's age (GEN.1213): `(rise_gyr, target_mol)`, `(None, 0.0)` when none has. The oxygenic
    classes' producers carry their summed share of the food production, so the source is that share of the body's (`biogenic_oxygen.productivity`); the
    store the first oxygen fills holds, besides the ocean's iron and the crust's reductants, the oxygen the methane of the classes before the rise takes
    (`standing_methane_mol`, two moles of oxygen to a mole of methane). The oxygen rises at the first producers' age plus the time the surplus takes to fill
    that store (`biogenic_oxygen.rise_age_gyr`), and then the air holds the balance of the summed source and the sinks.
    """
    pairs = _oxygenic(body, living)
    star_age = getattr(getattr(body, "star", None), "age", None)
    if not pairs or star_age is None:
        return None, 0.0
    rise = biogenic_oxygen.rise_age_gyr(body, pairs)
    if rise is None or rise > star_age:
        return None, 0.0
    extra = 2.0 * standing_methane_mol(body, living, rise)
    if extra > 0.0:
        rise = biogenic_oxygen.rise_age_gyr(body, pairs, extra)
        if rise is None or rise > star_age:
            return None, 0.0
    profile = biogenic_oxygen.profile_of(body, sum(share for _, share in pairs), extra)
    if profile.delay_gyr is None:
        return None, 0.0
    return rise, profile.steady_mol


def _producing(body, living, rose):
    """Whether the class has producers by the star's age and its gases have not ended (`rose`: the oxygen has risen on this world)."""
    methane = row_of(living.life_class).get("methane") or {}
    started = living.started(body)
    star_age = getattr(getattr(body, "star", None), "age", None)
    if started is None or star_age is None or star_age < started:
        return False
    return not (rose and methane.get("until") == "oxygen")


def apply_gases(body, living, rose):
    """
    Gives the air of `body` the gases of the ledger rows of the classes in `living` (a list of `Living`) at the star's age (GEN.1213): the methane their
    fluxes hold together (each class's flux times its share of the food production, so the fuel they share is counted once; the carbon dioxide it is made
    from taken out) and the hydrogen, hydrogen sulfide and the like each draws down, class by class. `rose` says the oxygen of a class that splits water has
    risen on this world (`biogenic_oxygen.apply` returned true): a row that ends there gives no methane. Applying it twice changes nothing more. Returns
    whether the air changed; the doses, aerosols, colour and score that read the air are the caller's to make again (`generation.life`).
    """
    if getattr(body, "body_type", "t") == "g":
        return False
    moles = biogenic_oxygen.air_moles(body)
    if not moles:
        return False
    changed = False
    new = dict(moles)
    oxygen_free = moles.get("O2", 0.0) <= number("oxygen_max_share") * sum(moles.values())
    flux = 0.0
    supplies = None
    for one in living:
        row = row_of(one.life_class)
        if row["status"] != "built" or not row.get("methane") or not _producing(body, one, rose):
            continue
        if row["methane"]["fuel"] != "burial" and supplies is None:
            supplies = life_supply.supplies(habitability_world.values_from_body(body))
        flux += one.share * methane_flux_mol_yr(body, one.life_class, supplies)
    if oxygen_free and flux > 0.0:
        have = moles.get("CH4", 0.0)
        wanted = methane_target_mol(flux, moles)
        pascal_per_mol = body.atmospheric_pressure / sum(moles.values())
        if wanted > have * number("apply_rise") and (wanted - have) * pascal_per_mol >= biogenic_oxygen.MINIMUM_PA:
            new["CH4"] = wanted
            new["CO2"] = moles["CO2"] - (wanted - have)
            changed = True
    for one in living:
        row = row_of(one.life_class)
        if row["status"] != "built" or not row["draws_down"] or not _producing(body, one, False):
            continue
        for gas in row["draws_down"]:
            share = number(DRAW_DOWN_CAPS[gas])
            limit = share / (1.0 - share) * (sum(new.values()) - new.get(gas, 0.0))
            have = new.get(gas, 0.0)
            if have <= limit * (1.0 + 1e-9):
                continue
            new[gas] = limit
            changed = True
            product = (row.get("products") or {}).get(gas)
            if product:      # what the life makes of the gas it takes: sulfur dioxide from hydrogen sulfide, methane from hydrogen
                new[product["gas"]] = new.get(product["gas"], 0.0) + (have - limit) * number(product["per"])
    if not changed:
        return False
    new = {formula: amount for formula, amount in new.items() if amount > 0.0}
    air = biogenic_oxygen.partials_pa(body, new)
    if not air:
        return False
    biogenic_oxygen.store_air(body, air)
    return True


HISTORY_GASES = (("O2", "Oxygen"), ("CO2", "Carbon dioxide"), ("CH4", "Methane"))
"""tuple: The gases the Life and air chart draws (`air_history`): the ones the oxygen and the methane rows of the ledger move."""


def _kpa(body, moles):
    """The partial pressures, kPa, of `moles` on `body` (formula to kPa)."""
    return {formula: pascal / 1000.0 for formula, pascal in biogenic_oxygen.partials_pa(body, moles).items()}


def air_history(body, living, stage_ages):
    """
    The air of `body` through the age of its star, as far as the life on it changed it (UX.237, GEN.1213, docs/design/habitability/life-classes.md section 6.2):
    worked backwards from the air it holds now, on demand, with nothing stored, as the sum of the classes in `living` (a list of `Living`, each with its
    share of the food production and the age its producers began). One step in time for each age something changed it: before the first producers the
    air as the life found it (the air now less what the ledger put in: the methane taken out, its carbon back in the carbon dioxide, and the oxygen taken out);
    from each class's producers on, the air with the methane the classes begun hold together (the summed flux, `methane_flux_mol_yr`, `methane_target_mol`)
    and no oxygen; and from the oxygen rise on (`oxygen_state`: the rise of the oxygenic classes' summed source, only if the star is old enough), the air
    now. The sulfur gases, the hydrogen the hydrogen-eaters draw down and the reduced gases the rise burned are not in it (the reductants of an abiotic air
    are not stored).

    Args:
        body: A rocky planet or moon with the life of `living` and the stored air.
        living (list): The `Living` classes of the body.
        stage_ages (dict): The age, Gyr, of each milestone by key (`evolution.MILESTONE_KEYS`), `None` for one never reached.

    Returns:
        dict or None: `age_gyr` (the star's age), `producers_gyr` (the first producers' age), `rose_gyr` (the age of the oxygen rise or `None`), `classes` (a
            list of `{"life_class", "variant", "share", "producers_gyr"}`, one for each class), `gases` (a list of `{"formula", "name", "steps": [[age_gyr,
            kPa], ...]}`, each value holding until the next age, only the gases that are ever above `HISTORY_MINIMUM_KPA`), `stages` (the `stage_ages` the
            star is old enough for, sorted by age) and `speculative` (the ledger's tier and tweak level for a class with `tweak` L3 or L4, else `None`).
            `None` when the body is a giant, has no air, no producers, or the air never differs from one age to another (nothing to draw).
    """
    star_age = getattr(getattr(body, "star", None), "age", None)
    begun = sorted((one.started(body), index) for index, one in enumerate(living) if one.started(body) is not None)
    if getattr(body, "body_type", "t") == "g" or not begun or not star_age or begun[0][0] >= star_age:
        return None
    now = biogenic_oxygen.air_moles(body)
    if not now:
        return None
    rise, _ = oxygen_state(body, living)
    if rise is not None and now.get("O2", 0.0) <= 0.0:
        rise = None
    supplies = None
    flux_of = []
    for one in living:
        row = row_of(one.life_class)
        fuel = (row.get("methane") or {}).get("fuel")
        if fuel and fuel != "burial" and supplies is None:
            supplies = life_supply.supplies(habitability_world.values_from_body(body))
        flux_of.append(one.share * methane_flux_mol_yr(body, one.life_class, supplies) if row["status"] == "built" and row.get("methane") else 0.0)
    if rise is not None:
        base = dict(now)
        base.pop("O2", None)
    else:
        base = dict(now)
    found = dict(base)                          # the air as the life found it: the methane's carbon back in the carbon dioxide
    methane_here = found.pop("CH4", 0.0)
    found["CO2"] = found.get("CO2", 0.0) + methane_here
    airs = [(0.0, found)]
    ages = sorted({age for age, _ in begun if age < star_age and (rise is None or age < rise)})
    for position, age in enumerate(ages):
        if rise is None and position == len(ages) - 1:
            airs.append((age, now))             # everything has begun and the oxygen has not risen: the air now
            continue
        flux = sum(flux_of[index] for started, index in begun if started <= age)
        air = dict(found)
        wanted = methane_target_mol(flux, air) if flux > 0.0 else 0.0
        if wanted > 0.0 and wanted * body.atmospheric_pressure / sum(air.values()) >= biogenic_oxygen.MINIMUM_PA:
            air["CH4"] = wanted
            air["CO2"] = max(air.get("CO2", 0.0) - wanted, 0.0)
        airs.append((age, air))
    if rise is not None:
        airs.append((rise, now))
    gases = []
    for formula, name in HISTORY_GASES:
        steps = [[age, round(_kpa(body, moles).get(formula, 0.0), 6)] for age, moles in airs]
        if max(value for _, value in steps) >= HISTORY_MINIMUM_KPA:
            gases.append({"formula": formula, "name": name, "steps": steps})
    if not gases or all(len({value for _, value in gas["steps"]}) == 1 for gas in gases):
        return None
    tuned = [row_of(one.life_class) for one in living if row_of(one.life_class)["tweak"] in ("L3", "L4")]
    return {"age_gyr": star_age, "producers_gyr": begun[0][0], "rose_gyr": rise,
            "classes": [{"life_class": one.life_class.number, "variant": one.life_class.variant or None, "share": round(one.share, 6), "producers_gyr": one.started(body)}
                        for one in living],
            "gases": gases,
            "stages": sorted(({"key": key, "age_gyr": age} for key, age in stage_ages.items() if age is not None and age <= star_age), key=lambda stage: stage["age_gyr"]),
            "speculative": {"tier": tuned[0]["tier"], "tweak": tuned[0]["tweak"]} if tuned else None}
