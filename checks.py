import json
import sys
from fractions import Fraction
import certificates as cert
from source import load


def must_reject(action):
    try:
        action()
    except ValueError:
        return
    raise AssertionError("invalid certificate accepted")


def verify():
    main = load("data/main")
    assert main["days"] == [Fraction(0), Fraction(24), Fraction(48)]
    assert len(main["readings"]) == 16 and len(main["routes"]) == 5
    a12 = next(r for r in main["readings"] if r["record"] == "a12")
    assert a12["value"] - a12["error"] == Fraction(271, 2)
    assert a12["value"] + a12["error"] == Fraction(273, 2)
    a_routes = [r for r in main["routes"] if r["branch"] == "A"]
    assert [(r["lo"], r["hi"], r["region"]) for r in a_routes] == [
        (Fraction(0), Fraction(30), "east"), (Fraction(30), Fraction(48), "west")]
    toy, bad = load("data/toy"), load("data/conflict")
    assert toy["readings"][1]["value"] == 2
    assert bad["readings"][1]["value"] == 3
    assert toy["regions"] == ["east", "south", "west"]
    rows = [{"id": "u", "a": {"x": "1"}, "b": "3/2"},
            {"id": "l", "a": {"x": "-1"}, "b": "0"}]
    assert cert.check_point(rows, {"x": "1/2"})
    assert cert.check_bound(rows, {"x": "2"}, "3", {"u": "2"})
    must_reject(lambda: cert.check_point(rows, {"x": "2"}))
    must_reject(lambda: cert.check_point(rows, {}))
    must_reject(lambda: cert.check_bound(rows, {"x": "1"}, "2", {"u": "1"}))
    must_reject(lambda: cert.check_bound(rows, {"x": "-1"}, "-3/2", {"u": "-1"}))
    impossible = [{"id": "p", "a": {"x": "1"}, "b": "0"},
                  {"id": "n", "a": {"x": "-1"}, "b": "-1"}]
    assert cert.check_conflict(impossible, {"p": "1", "n": "1"})
    must_reject(lambda: cert.check_conflict(rows, {"u": "1", "l": "1"}))
    must_reject(lambda: cert.check_conflict(impossible, {"unknown": "1"}))
    print(json.dumps({"main_readings": len(main["readings"]),
                      "a12_interval": ["271/2", "273/2"],
                      "switch_hour": "30", "generic_certificates": "checked"}))


def demo():
    pack = load("data/main")
    print(json.dumps({
        "readings": [{"record": r["record"], "meter": r["meter"], "hour": str(r["t"]),
                      "interval": [str(r["value"] - r["error"]),
                                   str(r["value"] + r["error"])]}
                     for r in pack["readings"]],
        "routes": [{"branch": r["branch"], "start": str(r["lo"]),
                    "end": str(r["hi"]), "region": r["region"],
                    "origin": r["origin"]} for r in pack["routes"]]
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    {"verify": verify, "demo": demo}[sys.argv[1]]()
