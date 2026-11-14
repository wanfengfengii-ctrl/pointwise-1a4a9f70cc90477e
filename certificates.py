"""Generic exact checks for finite rational linear inequalities."""
from fractions import Fraction


def q(value):
    if not isinstance(value, str):
        raise ValueError("certificate numbers must be strings")
    return Fraction(value)


def index(rows):
    result = {}
    for row in rows:
        if not isinstance(row["id"], str) or row["id"] in result:
            raise ValueError("duplicate or invalid inequality id")
        a = {key: q(value) for key, value in row["a"].items()}
        if any(not isinstance(key, str) for key in a):
            raise ValueError("invalid variable")
        result[row["id"]] = (a, q(row["b"]))
    return result


def check_point(rows, point):
    parsed = {key: q(value) for key, value in point.items()}
    for name, (a, b) in index(rows).items():
        if not set(a) <= set(parsed):
            raise ValueError("missing point coordinate")
        if sum((v * parsed[k] for k, v in a.items()), Fraction()) > b:
            raise ValueError("point violates " + name)
    return True


def combine(rows, weights):
    parsed = index(rows)
    if not set(weights) <= set(parsed):
        raise ValueError("unknown inequality")
    a, b = {}, Fraction()
    for name, value in weights.items():
        weight = q(value)
        if weight < 0:
            raise ValueError("negative multiplier")
        terms, rhs = parsed[name]
        b += weight * rhs
        for key, coefficient in terms.items():
            a[key] = a.get(key, Fraction()) + weight * coefficient
    return {key: value for key, value in a.items() if value}, b


def check_bound(rows, objective, upper, weights):
    a, b = combine(rows, weights)
    target = {key: q(value) for key, value in objective.items() if q(value)}
    if a != target or b != q(upper):
        raise ValueError("bound identity fails")
    return True


def check_conflict(rows, weights):
    a, b = combine(rows, weights)
    if a or b >= 0:
        raise ValueError("not a contradiction")
    return True
