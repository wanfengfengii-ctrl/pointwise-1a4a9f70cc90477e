"""Read sealed meter evidence; this module does not infer missing readings."""
import csv
import json
from fractions import Fraction
from pathlib import Path


def rational(value):
    if not isinstance(value, str):
        raise ValueError("rational values must be strings")
    return Fraction(value)


def load(directory):
    root = Path(directory)
    meta = json.loads((root / "meta.json").read_text(encoding="utf-8"))
    days = [rational(x) for x in meta["days"]]
    if len(days) != 3 or not days[0] < days[1] < days[2]:
        raise ValueError("two increasing day windows required")

    def table(name):
        with (root / name).open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle))
        if not rows or len({r["record"] for r in rows}) != len(rows):
            raise ValueError("empty table or duplicate record: " + name)
        return rows

    meters = table("meters.csv")
    ids = {r["meter"] for r in meters}
    if len(ids) != len(meters) or ids != {meta["trunk"], *meta["branches"]}:
        raise ValueError("meter directory mismatch")
    for row in meters:
        row["cap"] = rational(row["max_m3ph"])
        if row["cap"] < 0:
            raise ValueError("negative cap")
    readings = table("readings.csv")
    seen = set()
    for row in readings:
        row["t"] = rational(row["hour"])
        row["value"] = rational(row["display"])
        row["error"] = rational(row["half_width"])
        identity = (row["meter"], row["t"])
        if (row["meter"] not in ids or identity in seen or row["error"] < 0
                or not days[0] <= row["t"] <= days[-1]):
            raise ValueError("invalid reading")
        seen.add(identity)
    for meter in ids:
        if (meter, days[0]) not in seen:
            raise ValueError("initial counter missing")
    routes = table("routes.csv")
    regions = set()
    for row in routes:
        row["lo"], row["hi"] = rational(row["start"]), rational(row["end"])
        if (row["branch"] not in meta["branches"] or not row["region"]
                or not days[0] <= row["lo"] < row["hi"] <= days[-1]):
            raise ValueError("invalid route")
        regions.add(row["region"])
    for branch in meta["branches"]:
        parts = sorted((r["lo"], r["hi"]) for r in routes if r["branch"] == branch)
        cursor = days[0]
        for lo, hi in parts:
            if lo != cursor:
                raise ValueError("route gap or overlap")
            cursor = hi
        if cursor != days[-1]:
            raise ValueError("incomplete route")
    return {"meta": meta, "days": days, "meters": meters, "readings": readings,
            "routes": routes, "regions": sorted(regions)}
