#!/usr/bin/env python3
"""供水核账：封存两日片区用量的精确可达界、联合轨迹见证与有理证书。

用法:
  python3 analysis.py run MATERIAL_DIR OUTPUT_DIR     核定并写出 result.json 与 REPORT.md
  python3 analysis.py verify MATERIAL_DIR OUTPUT_DIR  从原材料重算并核验已有结果
  python3 analysis.py demo MATERIAL_DIR OUTPUT_DIR    重算并演示真实区间与一份端点证据

模型要点（详见 REPORT.md）：
  变量仅为各表期初累计真值 i_<表> 与有限时间段体积 v_<表>_<段>；
  断点取日界、全部读数时刻与全部换接时刻，段内路由恒定且无读数。
  行（a·x <= b）对应：读数量化闭区间、量程上限、累计非递减、段内干管=支路之和。
  同一读数在所有行中代表同一真值（期初+前缀段体积），跨日共享；
  换接只改变片区归属，不重置累计表。
  求解用 SciPy(HiGHS)，再由活跃集与对偶支撑精确重构有理顶点及 Farkas 乘子，
  最终证据一律经 certificates.py 有理核验；浮点重构失败时回退到精确单纯形。
"""
import json
import sys
import tempfile
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import certificates as cert
from source import load

ZERO = Fraction(0)
ONE = Fraction(1)


# ---------------------------------------------------------------- 模型

def build_model(pack):
    """从封存材料构造有限分段线性模型；返回 (model, segs)。segs 为 Fraction 时刻对。"""
    days = pack["days"]
    trunk = pack["meta"]["trunk"]
    branches = list(pack["meta"]["branches"])
    meter_order = [trunk] + branches
    cap = {m["meter"]: m["cap"] for m in pack["meters"]}
    mrec = {m["meter"]: m["record"] for m in pack["meters"]}

    pts = set(days)
    for r in pack["readings"]:
        pts.add(r["t"])
    for r in pack["routes"]:
        pts.add(r["lo"])
        pts.add(r["hi"])
    pts = sorted(pts)
    segs = [(pts[i], pts[i + 1]) for i in range(len(pts) - 1)]

    variables = {}
    for m in meter_order:
        variables["i_" + m] = "表%s在日界%s的期初累计真值(m3)" % (m, days[0])
    for k, (s, e) in enumerate(segs):
        for m in meter_order:
            variables["v_%s_%d" % (m, k)] = "表%s在时段%d [%s,%s) 的体积(m3)" % (m, k, s, e)

    rows = []

    def add(rid, a, b, sources):
        rows.append({"id": rid,
                     "a": {k: str(v) for k, v in a.items() if v},
                     "b": str(b),
                     "sources": list(sources)})

    # 读数量化闭区间：真值 = 期初 + 前缀段体积，同一读数跨日共享同一真值
    for r in sorted(pack["readings"],
                    key=lambda r: (meter_order.index(r["meter"]), r["t"], r["record"])):
        m, t = r["meter"], r["t"]
        a = {"i_" + m: ONE}
        for k, (s, e) in enumerate(segs):
            if e <= t:
                a["v_%s_%d" % (m, k)] = ONE
        src = ["readings.csv#" + r["record"]]
        add("rd_%s_hi" % r["record"], a, r["value"] + r["error"], src)
        add("rd_%s_lo" % r["record"], {k: -v for k, v in a.items()},
            -(r["value"] - r["error"]), src)

    # 量程上限（任意区间增量 <= 上限×时长，逐段表达即等价）与累计非递减
    for m in meter_order:
        for k, (s, e) in enumerate(segs):
            v = "v_%s_%d" % (m, k)
            src = ["meters.csv#" + mrec[m]]
            add("cap_%s_%d" % (m, k), {v: ONE}, cap[m] * (e - s), src)
            add("nn_%s_%d" % (m, k), {v: -ONE}, ZERO, src)

    # 每时刻干管流量等于各支路流量之和：对段积分得段体积守恒
    for k in range(len(segs)):
        a_le = {"v_%s_%d" % (trunk, k): ONE}
        a_ge = {"v_%s_%d" % (trunk, k): -ONE}
        for b in branches:
            a_le["v_%s_%d" % (b, k)] = -ONE
            a_ge["v_%s_%d" % (b, k)] = ONE
        add("bal_%d_le" % k, a_le, ZERO, ["meta.json"])
        add("bal_%d_ge" % k, a_ge, ZERO, ["meta.json"])

    model = {
        "days": [str(d) for d in days],
        "unit": pack["meta"].get("unit", "m3"),
        "trunk": trunk,
        "branches": branches,
        "regions": list(pack["regions"]),
        "variables": variables,
        "segments": [{"index": k, "start": str(s), "end": str(e)}
                     for k, (s, e) in enumerate(segs)],
        "rows": rows,
    }
    return model, segs


def build_objectives(pack, segs):
    """片区日用量 = 当日路由至该片区的支路段体积之和（线性系数映射）。"""
    days = pack["days"]
    branches = list(pack["meta"]["branches"])
    route_at = {}
    for b in branches:
        rts = sorted((r for r in pack["routes"] if r["branch"] == b),
                     key=lambda r: r["lo"])
        for k, (s, e) in enumerate(segs):
            hit = [r for r in rts if r["lo"] <= s and e <= r["hi"]]
            if len(hit) != 1:
                raise ValueError("route does not cover segment uniquely")
            route_at[(b, k)] = hit[0]["region"]
    objectives = {}
    for region in pack["regions"]:
        for d in (0, 1):
            obj = {}
            for b in branches:
                for k, (s, e) in enumerate(segs):
                    if days[d] <= s and e <= days[d + 1] and route_at[(b, k)] == region:
                        obj["v_%s_%d" % (b, k)] = ONE
            objectives[(region, d)] = obj
    return objectives


def parse_rows(rows_json):
    return [(r["id"], {k: Fraction(v) for k, v in r["a"].items()}, Fraction(r["b"]))
            for r in rows_json]


# ---------------------------------------------------------------- 精确线性代数

def gauss_square(mat, rhs):
    """精确解 n×n 线性方程组；奇异返回 None。"""
    n = len(mat)
    a = [list(mat[i]) + [rhs[i]] for i in range(n)]
    for c in range(n):
        piv = next((i for i in range(c, n) if a[i][c]), None)
        if piv is None:
            return None
        a[c], a[piv] = a[piv], a[c]
        pv = a[c][c]
        a[c] = [x / pv for x in a[c]]
        for i in range(n):
            if i != c and a[i][c]:
                f = a[i][c]
                a[i] = [x - f * y for x, y in zip(a[i], a[c])]
    return [a[i][n] for i in range(n)]


def gauss_under(mat, rhs):
    """精确求欠定方程组一组特解（主元列取最左可用列）；不相容返回 None。"""
    n_rows = len(mat)
    n_cols = len(mat[0]) if n_rows else 0
    a = [list(mat[i]) + [rhs[i]] for i in range(n_rows)]
    pivots = []
    r = 0
    for c in range(n_cols):
        piv = next((i for i in range(r, n_rows) if a[i][c]), None)
        if piv is None:
            continue
        a[r], a[piv] = a[piv], a[r]
        pv = a[r][c]
        a[r] = [x / pv for x in a[r]]
        for i in range(n_rows):
            if i != r and a[i][c]:
                f = a[i][c]
                a[i] = [x - f * y for x, y in zip(a[i], a[r])]
        pivots.append((c, r))
        r += 1
        if r == n_rows:
            break
    for i in range(r, n_rows):
        if a[i][n_cols]:
            return None
    sol = [ZERO] * n_cols
    for c, rr in pivots:
        sol[c] = a[rr][n_cols]
    return sol


# ---------------------------------------------------------------- 精确单纯形（回退与矛盾证书）

def simplex_exact(row_list, var_list, objective):
    """在 {x 自由 : a_i·x <= b_i} 上最大化 objective。

    返回 {"status":"optimal","value","point","weights"} 或
    {"status":"infeasible","farkas"}。全部 Fraction 精确运算，
    两阶段 + Bland 规则保证终止。weights/farkas 为行 id 到非负有理乘子。
    """
    m = len(row_list)
    nv = len(var_list)
    vidx = {v: j for j, v in enumerate(var_list)}
    flip = []
    amat = []
    bvec = []
    for (rid, a, bi) in row_list:
        row = [ZERO] * (2 * nv)
        for v, co in a.items():
            row[vidx[v]] = co
            row[nv + vidx[v]] = -co
        if bi < 0:
            row = [-z for z in row]
            bi = -bi
            flip.append(-1)
        else:
            flip.append(1)
        amat.append(row)
        bvec.append(bi)
    n_struct = 2 * nv
    art_of = {}
    ac = n_struct + m
    for i in range(m):
        if flip[i] < 0:
            art_of[i] = ac
            ac += 1
    n_cols = n_struct + m + len(art_of)
    tab = []
    for i in range(m):
        ti = [ZERO] * (n_cols + 1)
        for j in range(2 * nv):
            ti[j] = amat[i][j]
        if flip[i] > 0:
            ti[n_struct + i] = ONE
        else:
            ti[n_struct + i] = -ONE
            ti[art_of[i]] = ONE
        ti[n_cols] = bvec[i]
        tab.append(ti)
    basis = [(n_struct + i) if flip[i] > 0 else art_of[i] for i in range(m)]
    init_col = list(basis)  # 初始单位列；当前表中这些列即 B^{-1}
    art_cols = set(art_of.values())

    def pivot(r, j):
        pr = tab[r]
        pv = pr[j]
        if pv != 1:
            pr = [z / pv for z in pr]
            tab[r] = pr
        for i in range(m):
            if i != r:
                f = tab[i][j]
                if f:
                    tab[i] = [x - f * y for x, y in zip(tab[i], pr)]
        basis[r] = j

    def optimize(c, banned):
        while True:
            cb = [c[basis[r]] for r in range(m)]
            enter = -1
            for j in range(n_cols):
                if j in banned:
                    continue
                rc = c[j]
                for r in range(m):
                    if cb[r] and tab[r][j]:
                        rc -= cb[r] * tab[r][j]
                if rc > 0:
                    enter = j
                    break
            if enter < 0:
                val = sum((cb[r] * tab[r][n_cols] for r in range(m)), ZERO)
                return "optimal", val
            leave = -1
            best = None
            for i in range(m):
                tij = tab[i][enter]
                if tij > 0:
                    ratio = tab[i][n_cols] / tij
                    if (best is None or ratio < best
                            or (ratio == best and basis[i] < basis[leave])):
                        best = ratio
                        leave = i
            if leave < 0:
                return "unbounded", None
            pivot(leave, enter)

    # 阶段一：最大化 -Σ人工变量
    c1 = [ZERO] * n_cols
    for j in art_cols:
        c1[j] = -ONE
    st, val1 = optimize(c1, banned=set())
    assert st == "optimal"
    if val1 < 0:  # 不可行：取阶段一对偶为 Farkas 冲突乘子
        cb = [c1[basis[r]] for r in range(m)]
        farkas = {}
        for i in range(m):
            y = sum((cb[r] * tab[r][init_col[i]] for r in range(m)), ZERO)
            w = y if flip[i] > 0 else -y
            if w:
                assert w > 0
                farkas[row_list[i][0]] = w
        for v in var_list:
            s = sum((w * row_list[i][1].get(v, ZERO)
                     for i, (rid, _, _) in enumerate(row_list)
                     for w in [farkas.get(rid, ZERO)]), ZERO)
            assert s == 0
        rhs = sum((w * row_list[i][2]
                   for i, (rid, _, _) in enumerate(row_list)
                   for w in [farkas.get(rid, ZERO)]), ZERO)
        assert rhs < 0
        return {"status": "infeasible", "farkas": farkas}

    # 把零水平的人工变量换出基；全零行即冗余行，保留不影响后续
    for r in range(m):
        if basis[r] in art_cols:
            j = next((jj for jj in range(n_struct + m) if tab[r][jj]), None)
            if j is not None:
                pivot(r, j)
    # 阶段二
    c2 = [ZERO] * n_cols
    for v, co in objective.items():
        c2[vidx[v]] = co
        c2[nv + vidx[v]] = -co
    st, val = optimize(c2, banned=art_cols)
    if st != "optimal":
        raise ValueError("unbounded LP: model should be bounded by caps/readings")
    point = {}
    for v in var_list:
        xp = xq = ZERO
        for r in range(m):
            if basis[r] == vidx[v]:
                xp = tab[r][n_cols]
            elif basis[r] == nv + vidx[v]:
                xq = tab[r][n_cols]
        point[v] = xp - xq
    cb = [c2[basis[r]] for r in range(m)]
    weights = {}
    for i in range(m):
        y = sum((cb[r] * tab[r][init_col[i]] for r in range(m)), ZERO)
        w = y if flip[i] > 0 else -y
        if w:
            assert w > 0
            weights[row_list[i][0]] = w
    # 内部精确断言：原始可行、对偶恒等、强对偶
    for (rid, a, bi) in row_list:
        assert sum((co * point[v] for v, co in a.items()), ZERO) <= bi, rid
    for v in var_list:
        s = sum((w * row_list[i][1].get(v, ZERO)
                 for i, (rid, _, _) in enumerate(row_list)
                 for w in [weights.get(rid, ZERO)]), ZERO)
        assert s == objective.get(v, ZERO), v
    rhs = sum((w * row_list[i][2]
               for i, (rid, _, _) in enumerate(row_list)
               for w in [weights.get(rid, ZERO)]), ZERO)
    assert rhs == val
    assert sum((co * point[v] for v, co in objective.items()), ZERO) == val
    return {"status": "optimal", "value": val, "point": point, "weights": weights}


# ---------------------------------------------------------------- SciPy + 精确重构

class ReconstructError(Exception):
    pass


def scipy_optimize(row_list, var_list, objective):
    import numpy as np
    from scipy.optimize import linprog
    n = len(var_list)
    vidx = {v: j for j, v in enumerate(var_list)}
    a_ub = np.zeros((len(row_list), n))
    b_ub = np.zeros(len(row_list))
    for i, (rid, a, bi) in enumerate(row_list):
        for v, co in a.items():
            a_ub[i, vidx[v]] = float(co)
        b_ub[i] = float(bi)
    c = np.zeros(n)
    for v, co in objective.items():
        c[vidx[v]] = float(co)
    return linprog(-c, A_ub=a_ub, b_ub=b_ub,
                   bounds=[(None, None)] * n, method="highs")


def reconstruct(row_list, var_list, objective, res):
    """由 HiGHS 浮点解精确重构最优顶点与对偶乘子；任何精确检查失败都抛 ReconstructError。"""
    if res.status != 0 or res.ineqlin is None:
        raise ReconstructError("scipy status %s" % getattr(res, "status", "?"))
    n = len(var_list)
    vidx = {v: j for j, v in enumerate(var_list)}
    x = list(res.x)
    marg = list(res.ineqlin.marginals)  # min 问题乘子 <= 0；本问题权重 = -marg
    lhs = [sum(float(co) * x[vidx[v]] for v, co in a.items())
           for (_, a, _) in row_list]
    active = []
    for i, (rid, a, bi) in enumerate(row_list):
        if abs(lhs[i] - float(bi)) <= 1e-6 * max(1.0, abs(float(bi))):
            active.append(i)
    w_float = [-mf for mf in marg]
    active.sort(key=lambda i: (-w_float[i],
                               abs(lhs[i] - float(row_list[i][2]))))
    # 选 n 条精确线性无关的活跃行
    basis_rows = []
    selected = []
    for i in active:
        vec = [row_list[i][1].get(v, ZERO) for v in var_list]
        for pc, bv in basis_rows:
            if vec[pc]:
                f = vec[pc]
                vec = [a_ - f * b_ for a_, b_ in zip(vec, bv)]
        pc = next((j for j, z in enumerate(vec) if z), None)
        if pc is None:
            continue
        f = vec[pc]
        basis_rows.append((pc, [z / f for z in vec]))
        selected.append(i)
        if len(selected) == n:
            break
    if len(selected) < n:
        raise ReconstructError("active set rank %d < %d" % (len(selected), n))
    sol = gauss_square([[row_list[i][1].get(v, ZERO) for v in var_list]
                        for i in selected],
                       [row_list[i][2] for i in selected])
    if sol is None:
        raise ReconstructError("singular basis")
    point = {v: sol[j] for j, v in enumerate(var_list)}
    for (rid, a, bi) in row_list:
        if sum((co * point[v] for v, co in a.items()), ZERO) > bi:
            raise ReconstructError("primal infeasible at " + rid)
    beta = sum((co * point[v] for v, co in objective.items()), ZERO)
    if abs(float(beta) + float(res.fun)) > 1e-6 * max(1.0, abs(float(beta))):
        raise ReconstructError("value mismatch vs scipy")
    # 对偶：在活跃行（按浮点乘子降序）上精确解 A^T y = c
    cand = active
    meq = [[row_list[i][1].get(v, ZERO) for i in cand] for v in var_list]
    rhs = [objective.get(v, ZERO) for v in var_list]
    y = gauss_under(meq, rhs)
    if y is None:
        raise ReconstructError("dual system inconsistent")
    weights = {}
    for k, i in enumerate(cand):
        if y[k] < 0:
            raise ReconstructError("negative dual multiplier")
        if y[k]:
            weights[row_list[i][0]] = y[k]
    dual_val = sum((w * row_list[i][2]
                    for i, (rid, _, _) in enumerate(row_list)
                    for w in [weights.get(rid, ZERO)]), ZERO)
    if dual_val != beta:
        raise ReconstructError("duality gap")
    return {"status": "optimal", "value": beta, "point": point, "weights": weights}


def solve_bound(rows_json, var_list, objective):
    """求 max objective 的精确值、见证点与对偶乘子；不可行时给 Farkas 冲突乘子。"""
    row_list = parse_rows(rows_json)
    res = None
    try:
        res = scipy_optimize(row_list, var_list, objective)
    except Exception:
        res = None
    if res is not None and res.status == 2:  # HiGHS 判不可行：以精确单纯形结论为准
        return simplex_exact(row_list, var_list, objective)
    if res is not None and res.status == 0:
        try:
            return reconstruct(row_list, var_list, objective, res)
        except ReconstructError:
            pass
    return simplex_exact(row_list, var_list, objective)


# ---------------------------------------------------------------- 见证与核验

def make_witness(pack, segs, point, var_list):
    trunk = pack["meta"]["trunk"]
    meter_order = [trunk] + list(pack["meta"]["branches"])
    return {
        "point": {v: str(point[v]) for v in var_list},
        "initial": {m: str(point["i_" + m]) for m in meter_order},
        "segments": [
            {"start": str(s), "end": str(e),
             "volumes": {m: str(point["v_%s_%d" % (m, k)]) for m in meter_order}}
            for k, (s, e) in enumerate(segs)
        ],
    }


def usage_from_segments(pack, segments, region, day):
    """按原始 routes.csv 与轨迹段体积计算片区日用量（直接对原材料）。"""
    days = pack["days"]
    branches = list(pack["meta"]["branches"])
    total = ZERO
    for seg in segments:
        s, e = Fraction(seg["start"]), Fraction(seg["end"])
        if not (days[day] <= s and e <= days[day + 1]):
            continue
        for b in branches:
            hit = [r for r in pack["routes"]
                   if r["branch"] == b and r["lo"] <= s and e <= r["hi"]]
            if len(hit) != 1:
                raise ValueError("route coverage broken")
            if hit[0]["region"] == region:
                total += Fraction(seg["volumes"][b])
    return total


def check_trajectory(pack, model, witness, region=None, day=None, expect=None):
    """直接对原始材料核验完整物理轨迹：连续覆盖、量程、非递减、守恒、读数区间、
    point 与轨迹一致、达到声称端点。全部精确分数运算。"""
    days = pack["days"]
    trunk = pack["meta"]["trunk"]
    branches = list(pack["meta"]["branches"])
    meter_order = [trunk] + branches
    cap = {mm["meter"]: mm["cap"] for mm in pack["meters"]}
    segs = witness["segments"]
    if not segs:
        raise ValueError("empty trajectory")
    if Fraction(segs[0]["start"]) != days[0] or Fraction(segs[-1]["end"]) != days[-1]:
        raise ValueError("trajectory does not cover both days")
    for i in range(len(segs) - 1):
        if Fraction(segs[i]["end"]) != Fraction(segs[i + 1]["start"]):
            raise ValueError("trajectory not contiguous")
    point = {k: Fraction(v) for k, v in witness["point"].items()}
    if set(point) != set(model["variables"]):
        raise ValueError("point does not cover all model variables")
    for m in meter_order:
        if point["i_" + m] != Fraction(witness["initial"][m]):
            raise ValueError("initial inconsistent with point")
    for k, seg in enumerate(segs):
        s, e = Fraction(seg["start"]), Fraction(seg["end"])
        if not s < e:
            raise ValueError("degenerate segment")
        vols = {}
        for m in meter_order:
            v = Fraction(seg["volumes"][m])
            if point["v_%s_%d" % (m, k)] != v:
                raise ValueError("segment volume inconsistent with point")
            if v < 0:
                raise ValueError("counter decreased")
            if v > cap[m] * (e - s):
                raise ValueError("flow cap exceeded")
            vols[m] = v
        if vols[trunk] != sum((vols[b] for b in branches), ZERO):
            raise ValueError("trunk != sum of branches")
    for r in pack["readings"]:
        cv = Fraction(witness["initial"][r["meter"]])
        for seg in segs:
            if Fraction(seg["end"]) <= r["t"]:
                cv += Fraction(seg["volumes"][r["meter"]])
        if not (r["value"] - r["error"] <= cv <= r["value"] + r["error"]):
            raise ValueError("reading %s violated" % r["record"])
    cert.check_point(model["rows"], witness["point"])
    if region is not None:
        u = usage_from_segments(pack, segs, region, day)
        if u != Fraction(expect):
            raise ValueError("witness does not attain the bound")
    return True


def check_sources(pack, model):
    """模型每行的来源定位必须能在原材料中解析到真实记录。"""
    tables = {"readings.csv": {r["record"] for r in pack["readings"]},
              "meters.csv": {m["record"] for m in pack["meters"]},
              "routes.csv": {r["record"] for r in pack["routes"]}}
    for row in model["rows"]:
        if not row["sources"]:
            raise ValueError("row without sources: " + row["id"])
        for src in row["sources"]:
            if "#" in src:
                fname, rec = src.split("#", 1)
                if rec not in tables.get(fname, set()):
                    raise ValueError("unresolvable source " + src)
            elif src != "meta.json":
                raise ValueError("unknown source " + src)
    return True


# ---------------------------------------------------------------- 核定

def make_item(pack, model, segs, region, day, obj, lo, up):
    var_list = list(model["variables"])
    lower = -lo["value"]
    upper = up["value"]
    neg = {k: -v for k, v in obj.items()}
    item = {
        "region": region,
        "day": day,
        "lower": str(lower),
        "upper": str(upper),
        "lower_witness": make_witness(pack, segs, lo["point"], var_list),
        "upper_witness": make_witness(pack, segs, up["point"], var_list),
        "lower_proof": {"objective": {k: str(v) for k, v in neg.items()},
                        "upper": str(lo["value"]),
                        "weights": {k: str(v) for k, v in lo["weights"].items()}},
        "upper_proof": {"objective": {k: str(v) for k, v in obj.items()},
                        "upper": str(up["value"]),
                        "weights": {k: str(v) for k, v in up["weights"].items()}},
    }
    # 出厂前自检：证书与见证必须立即通过通用精确核验
    cert.check_point(model["rows"], item["lower_witness"]["point"])
    cert.check_point(model["rows"], item["upper_witness"]["point"])
    cert.check_bound(model["rows"], item["lower_proof"]["objective"],
                     item["lower_proof"]["upper"], item["lower_proof"]["weights"])
    cert.check_bound(model["rows"], item["upper_proof"]["objective"],
                     item["upper_proof"]["upper"], item["upper_proof"]["weights"])
    if sum((obj[v] * lo["point"][v] for v in obj), ZERO) != lower:
        raise AssertionError("lower witness off bound")
    if sum((obj[v] * up["point"][v] for v in obj), ZERO) != upper:
        raise AssertionError("upper witness off bound")
    return item


def assess(pack):
    """对全部片区×日求精确可达界。返回 (model, segs, objectives, items, farkas)。"""
    model, segs = build_model(pack)
    objectives = build_objectives(pack, segs)
    var_list = list(model["variables"])
    items = []
    for region in pack["regions"]:
        for d in (0, 1):
            obj = objectives[(region, d)]
            up = solve_bound(model["rows"], var_list, obj)
            if up["status"] == "infeasible":
                return model, segs, objectives, None, up["farkas"]
            lo = solve_bound(model["rows"], var_list, {k: -v for k, v in obj.items()})
            if lo["status"] == "infeasible":
                return model, segs, objectives, None, lo["farkas"]
            items.append(make_item(pack, model, segs, region, d, obj, lo, up))
    return model, segs, objectives, items, None


def tightened_bounds(pack, model, objectives, items):
    """对照：去掉守恒行后逐表独立的界；联合模型的界收紧之处。"""
    rows_nb = [r for r in model["rows"] if not r["id"].startswith("bal_")]
    var_list = list(model["variables"])
    out = []
    for item in items:
        region, day = item["region"], item["day"]
        obj = objectives[(region, day)]
        up = solve_bound(rows_nb, var_list, obj)
        lo = solve_bound(rows_nb, var_list, {k: -v for k, v in obj.items()})
        nl, nu = -lo["value"], up["value"]
        jl, ju = Fraction(item["lower"]), Fraction(item["upper"])
        if nl < jl or nu > ju:
            out.append({"region": region, "day": day,
                        "naive": [str(nl), str(nu)], "joint": [str(jl), str(ju)]})
    return out


def conflict_materials(model, farkas):
    srcs = set()
    by_id = {r["id"]: r for r in model["rows"]}
    for rid, w in farkas.items():
        if w > 0:
            srcs.update(by_id[rid]["sources"])
    return sorted(srcs)


# ---------------------------------------------------------------- 报告

def render_report(pack, model, segs, items, tightened):
    days = pack["days"]
    trunk = pack["meta"]["trunk"]
    branches = list(pack["meta"]["branches"])
    rdg = {r["record"]: r for r in pack["readings"]}
    mts = {m["record"]: m for m in pack["meters"]}
    rts = {r["record"]: r for r in pack["routes"]}
    ex_rd = next((r for r in pack["readings"] if r["error"] > 0), pack["readings"][0])
    sw = [r for r in pack["routes"] if r["lo"] > days[0]]

    L = []
    L.append("# 供水核账报告：两日片区用量核定")
    L.append("")
    L.append("材料目录封存了 %s 至 %s 的两个日窗口 [%s,%s)、[%s,%s)（单位 %s）。"
             "干管 %s，支管 %s；片区：%s。"
             % (days[0], days[-1], days[0], days[1], days[1], days[2],
                pack["meta"].get("unit", "m3"), trunk, "、".join(branches),
                "、".join(pack["regions"])))
    L.append("")
    L.append("## 核定区间（机器可核对）")
    L.append("")
    L.append("| 片区 | 日 | 下界 | 上界 |")
    L.append("| --- | --- | --- | --- |")
    for it in items:
        L.append("| %s | %d | %s | %s |" % (it["region"], it["day"],
                                             it["lower"], it["upper"]))
    L.append("")
    L.append("每个端点都是精确可达界：result.json 为其附一份覆盖全部表与两日的"
             "完整轨迹见证（initial + 连续 segments，段内按恒定流率解释）和一份"
             "有理数 Farkas 证书（proof.weights 为非负有理乘子，加权全部约束行后"
             "左侧恰等于该日片区用量或其相反数、右侧恰等于对应界）。"
             "端点不是逐表逐日局部极值的拼接，而是同一组联合可行轨迹。")
    L.append("")
    L.append("## 材料与模型")
    L.append("")
    L.append("- 读数量化：readings.csv 每行给示值 display 与半宽 half_width，真值属于"
             "闭区间 [display−half_width, display+half_width]；例如 readings.csv#%s "
             "为 [%s, %s]。同一记录的真值跨日共享，不为不同差分重复选择。"
             % (ex_rd["record"], ex_rd["value"] - ex_rd["error"],
                ex_rd["value"] + ex_rd["error"]))
    L.append("- 量程：meters.csv 给每表恒定流量上限，任意区间增量不超过上限×时长；"
             "例如 %s。" % "、".join("meters.csv#%s 上限 %s m3/h" % (k, mts[k]["cap"])
                                      for k in sorted(mts)))
    sw_txt = "；".join("routes.csv#%s：%s [%s,%s)→%s" % (k, rts[k]["branch"],
                       rts[k]["lo"], rts[k]["hi"], rts[k]["region"]) for k in sorted(rts))
    L.append("- 路由与换接：%s。换接时刻已知，累计表不随换接重置。" % sw_txt)
    L.append("- 变量：i_<表> 为日界 %s 的期初累计真值；v_<表>_<段> 为该表在对应时段的"
             "体积。断点取日界、全部读数时刻与全部换接时刻，共 %d 段：%s。"
             % (days[0], len(segs),
                "、".join("[%s,%s)" % (s, e) for s, e in segs)))
    L.append("- 约束行（result.json model.rows，每行附 sources 材料定位）："
             "rd_* 读数量化闭区间（计量规则：量化闭区间）；cap_* 量程上限"
             "（任意区间增量≤上限×时长，逐段表达即等价）；nn_* 累计非递减；"
             "bal_* 段内干管体积=各支路体积之和（每时刻干管流量=支路流量之和，"
             "对段积分）。模型精确保留材料允许的读数与片区日用量联合可行集合，"
             "未加入插值、整数或误差独立性假设。")
    L.append("")
    L.append("## 有限分段为何不遗漏任何可行的日量与读数组合")
    L.append("")
    L.append("段内无读数且路由恒定。(1) 任一真实物理轨迹（连续、非递减、任意区间"
             "增量≤上限×时长、每时刻干管流量=支路流量之和）在各段上的积分体积必满足"
             "全部行：读数时刻都是断点，故读数行精确表达量化闭区间；任意区间规则"
             "含段区间，故有量程行；非递减给出 nn 行；流量守恒对段积分给出 bal 行。"
             "片区日用量是段体积的线性组合（日界与换接都是断点）。")
    L.append("(2) 反之，任一行可行点把每段体积按恒定流率 体积/段长 摊到段内："
             "流率≤上限，故任意区间增量≤上限×时长；累计折线连续非递减；"
             "bal 行使守恒在每时刻成立；读数行使各断点累计真值落在量化闭区间内。"
             "于是行可行集与连续轨迹在（期初值、段体积）上的投影完全一致，"
             "日用量与读数真值的联合可达集合不变；对线性目标求得的最值即连续问题的"
             "精确可达界，Farkas 证书对整个连续可行集有效。")
    L.append("")
    L.append("## 联合约束如何收紧片区范围")
    L.append("")
    if tightened:
        L.append("对照“逐表独立”（去掉 bal_* 守恒行，各表仅用自身读数与量程）的界，"
                 "联合模型因干管读数（如 %s）与每时刻守恒、以及同一读数跨日共享"
                 "而收紧以下片区日量：" % "、".join(
                     "readings.csv#%s" % r["record"] for r in pack["readings"]
                     if r["meter"] == trunk))
        L.append("")
        L.append("| 片区 | 日 | 逐表独立界 | 联合核定界 |")
        L.append("| --- | --- | --- | --- |")
        for t in tightened:
            L.append("| %s | %d | [%s, %s] | [%s, %s] |"
                     % (t["region"], t["day"], t["naive"][0], t["naive"][1],
                        t["joint"][0], t["joint"][1]))
        L.append("")
        L.append("收紧机理：干管各段体积被其读数夹住后，守恒把它分摊到三支路；"
                 "某支路被自身读数与量程限制时，其余支路的段体积随之被夹紧；"
                 "跨日共享的读数（如 %s）再把相邻时段的可行分配联动起来。"
                 % "、".join("readings.csv#%s" % r["record"] for r in pack["readings"]
                             if days[0] < r["t"] < days[-1]))
    else:
        L.append("本批材料下联合约束未产生超出逐表独立界的额外收紧。")
    L.append("")
    L.append("## 为何不能唯一核定")
    L.append("")
    if any(r["error"] > 0 for r in pack["readings"]):
        L.append("- 量化误差无概率分布：readings.csv 各行 half_width>0（如 readings.csv#%s），"
                 "真值只需落在闭区间内，不能按示值当精确读数，也不输出置信区间。"
                 % ex_rd["record"])
    else:
        L.append("- 本批读数半宽为零，但读数只覆盖部分时刻；真值在缺测时段内仍自由，"
                 "不输出置信区间。")
    miss = [m for m in branches + [trunk]
            if {r["t"] for r in pack["readings"] if r["meter"] == m}
            != set(days[1:-1]) | {days[0], days[-1]}]
    L.append("- 缺测小时不是零：源日志只保留 readings.csv 出现的时刻（如表 %s 在部分"
             "日界/换接时刻无读数），段内分配自由，不得以插值核定缺测量。" % "、".join(miss))
    if sw:
        L.append("- 换接重分配：%s 等换接改变支路去向而不重置累计表，片区归属的体积"
                 "在换接前后重新划分，单表累计差分不能唯一确定片区日量。"
                 % "、".join("routes.csv#%s" % r["record"] for r in sw))
    else:
        L.append("- 本批无换接；但守恒只约束干管=支路之和，不约束支路之间的分配，"
                 "单表累计差分不能唯一确定片区日量。")
    L.append("- 守恒只约束干管=支路之和，不约束支路之间的分配；同一读数跨日共享"
             "只能联合收紧，不能逐表逐日独立确定。")
    L.append("")
    L.append("## 证据与复核")
    L.append("")
    L.append("- 每个端点：完整轨迹（全部表、两日连续覆盖）+ 有理证书；"
             "certificates.check_point / check_bound 精确核验；"
             "verify 还会绕过模型行，直接对原始 readings/meters/routes 复核轨迹物理性"
             "（连续覆盖、量程、非递减、守恒、读数区间、达到端点）。")
    L.append("- 小例互核：data/toy 期望 east、west 每日 [0,1]、south 每日 [0,0]；"
             "data/conflict 每表单独可解释而联合矛盾；另自建“换接与读数同刻”与"
             "“全部日量唯一”两个微型对照，用完整枚举整数轨迹与求解结果互核。"
             "微型例中连续可行域的极值由整数轨迹覆盖，是因为有理证书已对整个连续"
             "可行集有效、而枚举的整数轨迹达到了同一界；这不声称一般材料只有整数解。")
    L.append("- 重跑：python3 analysis.py run data/main outputs/main；"
             "核验：python3 analysis.py verify data/main outputs/main；"
             "演示：python3 analysis.py demo data/main outputs/main。")
    L.append("")
    return "\n".join(L)


def render_conflict_report(pack, model, farkas, materials):
    days = pack["days"]
    by_id = {r["id"]: r for r in model["rows"]}
    loc = {}
    for r in pack["readings"]:
        loc["readings.csv#" + r["record"]] = "表%s于时刻%s示值%s±%s，真值∈[%s,%s]" % (
            r["meter"], r["t"], r["value"], r["error"],
            r["value"] - r["error"], r["value"] + r["error"])
    for m in pack["meters"]:
        loc["meters.csv#" + m["record"]] = "表%s流量上限 %s m3/h" % (m["meter"], m["cap"])
    loc["meta.json"] = "干管 %s 与支管 %s 的拓扑（每时刻干管流量=支路流量之和）" % (
        pack["meta"]["trunk"], "、".join(pack["meta"]["branches"]))
    L = []
    L.append("# 供水核账报告：材料联合矛盾")
    L.append("")
    L.append("两日窗口 [%s,%s)、[%s,%s) 的封存材料不存在任何一组同时满足全部计量规则"
             "的轨迹：虽然逐表单独看都可解释，但联合（读数量化闭区间 + 量程 + 非递减 + "
             "干管=支路守恒）后无可行解。" % (days[0], days[1], days[1], days[2]))
    L.append("")
    L.append("冲突由 Farkas 乘子证明：result.json conflict.weights 给每条约束行一个"
             "非负有理乘子，加权后所有变量系数恰为零、右侧严格小于零"
             "（certificates.check_conflict 通过）。正权重行涉及的材料如下。")
    L.append("")
    L.append("## 冲突记录（机器可核对）")
    L.append("")
    L.append("| 材料定位 |")
    L.append("| --- |")
    for s in materials:
        L.append("| %s |" % s)
    L.append("")
    L.append("## 冲突说明")
    L.append("")
    for s in materials:
        L.append("- %s：%s" % (s, loc.get(s, "（见原表）")))
    L.append("")
    L.append("正权重约束行：%s。"
             % "、".join("%s(×%s)" % (rid, w) for rid, w in sorted(farkas.items())))
    L.append("")
    L.append("核验：python3 analysis.py verify MATERIAL_DIR OUTPUT_DIR；"
             "演示：python3 analysis.py demo MATERIAL_DIR OUTPUT_DIR。")
    L.append("")
    return "\n".join(L)


# ---------------------------------------------------------------- 读写结果

def no_json_float(x):
    if isinstance(x, float):
        raise ValueError("JSON float forbidden")
    if isinstance(x, dict):
        for v in x.values():
            no_json_float(v)
    elif isinstance(x, (list, tuple)):
        for v in x:
            no_json_float(v)


def write_outputs(out_dir, result, report):
    no_json_float(result)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    (out / "result.json").write_text(json.dumps(result, ensure_ascii=False, indent=1),
                                     encoding="utf-8")
    (out / "REPORT.md").write_text(report, encoding="utf-8")


def run(material_dir, output_dir):
    try:
        pack = load(material_dir)
    except Exception as exc:
        print("读取材料失败: %s" % exc, file=sys.stderr)
        return 1
    model, segs, objectives, items, farkas = assess(pack)
    if items is None:
        materials = conflict_materials(model, farkas)
        result = {"status": "inconsistent", "model": model,
                  "conflict": {"weights": {k: str(v) for k, v in farkas.items()},
                               "materials": materials}}
        cert.check_conflict(model["rows"], result["conflict"]["weights"])
        write_outputs(output_dir, result,
                      render_conflict_report(pack, model, farkas, materials))
    else:
        tightened = tightened_bounds(pack, model, objectives, items)
        result = {"status": "feasible", "model": model, "bounds": items}
        write_outputs(output_dir, result,
                      render_report(pack, model, segs, items, tightened))
    return 0


# ---------------------------------------------------------------- 整数穷举（微型互核）

def enumerate_integer(pack, segs, objectives):
    """完整枚举整数段体积轨迹（要求零半宽、整数材料），返回 {(region,day):(min,max), ...}
    与可行轨迹条数。"""
    days = pack["days"]
    trunk = pack["meta"]["trunk"]
    branches = list(pack["meta"]["branches"])
    meter_order = [trunk] + branches
    cap = {m["meter"]: m["cap"] for m in pack["meters"]}
    for r in pack["readings"]:
        if r["error"] != 0 or r["value"].denominator != 1:
            raise ValueError("enumeration needs exact integer readings")
    init = {}
    for m in meter_order:
        r0 = [r for r in pack["readings"] if r["meter"] == m and r["t"] == days[0]]
        if len(r0) != 1:
            raise ValueError("initial reading missing")
        init[m] = r0[0]["value"]
    keys = [(m, k) for k in range(len(segs)) for m in meter_order]
    grids = []
    for (m, k) in keys:
        s, e = segs[k]
        hi = cap[m] * (e - s)
        if hi.denominator != 1 or hi > 64:
            raise ValueError("enumeration grid too large")
        grids.append(range(int(hi) + 1))
    rdg = {}
    for r in pack["readings"]:
        rdg.setdefault(r["meter"], []).append(r)
    best = {}
    count = 0
    import itertools
    for assign in itertools.product(*grids):
        vol = dict(zip(keys, assign))
        ok = True
        for k in range(len(segs)):
            if vol[(trunk, k)] != sum(vol[(b, k)] for b in branches):
                ok = False
                break
        if not ok:
            continue
        for m, rs in rdg.items():
            for r in rs:
                cv = init[m] + sum(vol[(m, k)] for k, (s, e) in enumerate(segs)
                                   if e <= r["t"])
                if cv != r["value"]:
                    ok = False
                    break
            if not ok:
                break
        if not ok:
            continue
        count += 1
        segs_json = [{"start": str(s), "end": str(e),
                      "volumes": {m: str(vol[(m, k)]) for m in meter_order}}
                     for k, (s, e) in enumerate(segs)]
        for (region, d) in objectives:
            u = usage_from_segments(pack, segs_json, region, d)
            key = (region, d)
            if key not in best:
                best[key] = [u, u]
            else:
                best[key][0] = min(best[key][0], u)
                best[key][1] = max(best[key][1], u)
    return best, count


# ---------------------------------------------------------------- 微型对照材料

def write_material(root, meta, meters, readings, routes):
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    (root / "meta.json").write_text(json.dumps(meta), encoding="utf-8")

    def wcsv(name, header, rows):
        lines = [",".join(header)]
        for r in rows:
            lines.append(",".join(str(x) for x in r))
        (root / name).write_text("\n".join(lines) + "\n", encoding="utf-8")

    wcsv("meters.csv", ["record", "meter", "max_m3ph", "origin"], meters)
    wcsv("readings.csv", ["record", "meter", "hour", "display", "half_width", "origin"],
         readings)
    wcsv("routes.csv", ["record", "branch", "start", "end", "region", "origin"], routes)


def micro_switch_reading(root):
    """换接与读数同刻：换接时刻 1 同时是干管读数时刻。"""
    write_material(
        root,
        {"days": ["0", "1", "2"], "trunk": "T", "branches": ["A", "B"], "unit": "m3"},
        [("mT", "T", "4", "micro1:capT"), ("mA", "A", "2", "micro1:capA"),
         ("mB", "B", "2", "micro1:capB")],
        [("t0", "T", "0", "0", "0", "micro1:T0"), ("t1", "T", "1", "2", "0", "micro1:T1"),
         ("t2", "T", "2", "3", "0", "micro1:T2"),
         ("a0", "A", "0", "0", "0", "micro1:A0"), ("a2", "A", "2", "2", "0", "micro1:A2"),
         ("b0", "B", "0", "0", "0", "micro1:B0"), ("b2", "B", "2", "1", "0", "micro1:B2")],
        [("rA0", "A", "0", "1", "east", "micro1:v0"), ("rA1", "A", "1", "2", "west", "micro1:v1"),
         ("rB0", "B", "0", "1", "west", "micro1:v0"), ("rB1", "B", "1", "2", "east", "micro1:v1")])


def micro_unique(root):
    """全部日量唯一：每个片区每日用量下界=上界。"""
    write_material(
        root,
        {"days": ["0", "1", "2"], "trunk": "T", "branches": ["A", "B", "C"], "unit": "m3"},
        [("mT", "T", "2", "micro2:capT"), ("mA", "A", "1", "micro2:capA"),
         ("mB", "B", "1", "micro2:capB"), ("mC", "C", "0", "micro2:capC")],
        [("t0", "T", "0", "0", "0", "micro2:T0"), ("t1", "T", "1", "1", "0", "micro2:T1"),
         ("t2", "T", "2", "2", "0", "micro2:T2"),
         ("a0", "A", "0", "0", "0", "micro2:A0"), ("a1", "A", "1", "1", "0", "micro2:A1"),
         ("a2", "A", "2", "1", "0", "micro2:A2"),
         ("b0", "B", "0", "0", "0", "micro2:B0"), ("b1", "B", "1", "0", "0", "micro2:B1"),
         ("b2", "B", "2", "1", "0", "micro2:B2"),
         ("c0", "C", "0", "0", "0", "micro2:C0"), ("c1", "C", "1", "0", "0", "micro2:C1"),
         ("c2", "C", "2", "0", "0", "micro2:C2")],
        [("rA0", "A", "0", "1", "east", "micro2:v0"), ("rA1", "A", "1", "2", "west", "micro2:v1"),
         ("rB0", "B", "0", "1", "west", "micro2:v0"), ("rB1", "B", "1", "2", "east", "micro2:v1"),
         ("rC", "C", "0", "2", "south", "micro2:fx")])


# ---------------------------------------------------------------- 核验

def check_item_full(pack, model, objectives, item):
    """核验一个界条目：目标一致、证书精确、轨迹物理可行且达到端点。"""
    region, day = item["region"], item["day"]
    obj = objectives[(region, day)]
    neg = {k: -v for k, v in obj.items()}
    if item["upper_proof"]["objective"] != {k: str(v) for k, v in obj.items()}:
        raise ValueError("upper objective mismatch")
    if item["lower_proof"]["objective"] != {k: str(v) for k, v in neg.items()}:
        raise ValueError("lower objective mismatch")
    if Fraction(item["lower_proof"]["upper"]) != -Fraction(item["lower"]):
        raise ValueError("lower proof upper != -lower")
    if Fraction(item["upper_proof"]["upper"]) != Fraction(item["upper"]):
        raise ValueError("upper proof upper != upper")
    cert.check_bound(model["rows"], item["lower_proof"]["objective"],
                     item["lower_proof"]["upper"], item["lower_proof"]["weights"])
    cert.check_bound(model["rows"], item["upper_proof"]["objective"],
                     item["upper_proof"]["upper"], item["upper_proof"]["weights"])
    check_trajectory(pack, model, item["lower_witness"], region, day, item["lower"])
    check_trajectory(pack, model, item["upper_witness"], region, day, item["upper"])
    return True


def parse_md_table(text, heading):
    lines = text.splitlines()
    out = []
    in_sec = False
    for ln in lines:
        if ln.startswith("## "):
            in_sec = ln.startswith(heading)
            continue
        if in_sec and ln.strip().startswith("|"):
            cells = [c.strip() for c in ln.strip().strip("|").split("|")]
            if not cells or cells[0] in ("片区", "材料定位"):
                continue
            if set(cells[0]) <= set("-: "):
                continue
            out.append(cells)
    return out


def verify_feasible_dir(pack, result, report_text, expect_bounds=None):
    """重算并核验一个 feasible 结果；返回条目数。"""
    model, segs, objectives, items, farkas = assess(pack)
    if items is None:
        raise ValueError("recompute says inconsistent, stored says feasible")
    if result.get("status") != "feasible":
        raise ValueError("stored status not feasible")
    if result.get("model") != model:
        raise ValueError("stored model rows differ from raw-material rebuild")
    check_sources(pack, result["model"])
    got = {(it["region"], it["day"]): (it["lower"], it["upper"]) for it in result["bounds"]}
    for it in items:
        key = (it["region"], it["day"])
        if got.get(key) != (it["lower"], it["upper"]):
            raise ValueError("bound mismatch at %s" % (key,))
    if expect_bounds is not None:
        for key, val in expect_bounds.items():
            if got.get(key) != val:
                raise ValueError("expected bound %s at %s, got %s" % (val, key, got.get(key)))
    stored = {(it["region"], it["day"]): it for it in result["bounds"]}
    for it in items:
        check_item_full(pack, model, objectives, stored[(it["region"], it["day"])])
    tbl = parse_md_table(report_text, "## 核定区间")
    if len(tbl) != len(items):
        raise ValueError("report bounds table size mismatch")
    for row in tbl:
        region, day, lo, hi = row[0], int(row[1]), row[2], row[3]
        key = (region, day)
        if got.get(key) != (lo, hi):
            raise ValueError("report table mismatch at %s" % (key,))
    return len(items)


def verify_inconsistent(pack, model_rows, weights):
    """核验矛盾：证书精确、重算不可行、逐表单独可解释。"""
    cert.check_conflict(model_rows, weights)
    row_list = parse_rows(model_rows)
    var_list = sorted({v for (_, a, _) in row_list for v in a})
    st = simplex_exact(row_list, var_list, {})
    if st["status"] != "infeasible":
        raise ValueError("recompute says feasible, stored says inconsistent")
    # 逐表单独均可解释
    for m in [pack["meta"]["trunk"]] + list(pack["meta"]["branches"]):
        sub = [r for r in row_list
               if r[0].startswith("rd_") and any(
                   rd["record"] == r[0][3:-3] and rd["meter"] == m
                   for rd in pack["readings"])
               or r[0].startswith("cap_%s_" % m) or r[0].startswith("nn_%s_" % m)]
        st_m = simplex_exact(sub, var_list, {})
        if st_m["status"] != "optimal":
            raise ValueError("meter %s alone not explainable" % m)
    return True


def verify(material_dir, output_dir):
    try:
        pack = load(material_dir)
    except Exception as exc:
        print("读取材料失败: %s" % exc, file=sys.stderr)
        return 1
    out = Path(output_dir)
    try:
        result = json.loads((out / "result.json").read_text(encoding="utf-8"))
        report_text = (out / "REPORT.md").read_text(encoding="utf-8")
    except Exception as exc:
        print("读取已有结果失败: %s" % exc, file=sys.stderr)
        return 1
    no_json_float(result)
    summary = {"main": None, "toy": None, "conflict": None, "micro": []}

    if result.get("status") == "feasible":
        n = verify_feasible_dir(pack, result, report_text)
        summary["main"] = {"status": "feasible", "bounds_checked": n}
    elif result.get("status") == "inconsistent":
        model, segs = build_model(pack)
        if result.get("model") != model:
            raise ValueError("stored model rows differ from raw-material rebuild")
        check_sources(pack, result["model"])
        verify_inconsistent(pack, result["model"]["rows"],
                            result["conflict"]["weights"])
        mats = conflict_materials(model, {k: Fraction(v) for k, v in
                                          result["conflict"]["weights"].items()})
        if sorted(result["conflict"]["materials"]) != mats:
            raise ValueError("conflict materials mismatch")
        tbl = parse_md_table(report_text, "## 冲突记录")
        if sorted(row[0] for row in tbl) != mats:
            raise ValueError("report conflict table mismatch")
        summary["main"] = {"status": "inconsistent", "materials": len(mats)}
    else:
        raise ValueError("unknown status in result.json")

    # 契约小例：data/toy 期望区间；data/conflict 联合矛盾
    base = Path(material_dir).parent
    toy = load(base / "toy")
    t_model, t_segs, t_obj, t_items, t_farkas = assess(toy)
    if t_items is None:
        raise ValueError("data/toy should be feasible")
    expect = {}
    for region in ("east", "west"):
        for d in (0, 1):
            expect[(region, d)] = ("0", "1")
    for d in (0, 1):
        expect[("south", d)] = ("0", "0")
    got = {(it["region"], it["day"]): (it["lower"], it["upper"]) for it in t_items}
    if got != expect:
        raise ValueError("data/toy bounds mismatch: %s" % got)
    for it in t_items:
        check_item_full(toy, t_model, t_obj, it)
    enum, cnt = enumerate_integer(toy, t_segs, t_obj)
    for key, (lo, hi) in enum.items():
        if (str(lo), str(hi)) != expect[key]:
            raise ValueError("toy enumeration mismatch at %s" % (key,))
    summary["toy"] = {"bounds": len(t_items), "integer_trajectories": cnt}

    bad = load(base / "conflict")
    b_model, b_segs, b_obj, b_items, b_farkas = assess(bad)
    if b_items is not None:
        raise ValueError("data/conflict should be inconsistent")
    cert.check_conflict(b_model["rows"], {k: str(v) for k, v in b_farkas.items()})
    verify_inconsistent(bad, b_model["rows"], {k: str(v) for k, v in b_farkas.items()})
    summary["conflict"] = {"farkas_rows": len(b_farkas),
                           "materials": conflict_materials(b_model, b_farkas)}

    # 自建微型对照：换接与读数同刻；全部日量唯一。整数穷举与求解互核。
    tmp = Path(tempfile.mkdtemp(prefix="micro_"))
    micro_switch_reading(tmp / "switch_reading")
    micro_unique(tmp / "unique")
    mp = load(tmp / "switch_reading")
    m_model, m_segs, m_obj, m_items, m_farkas = assess(mp)
    if m_items is None:
        raise ValueError("micro switch_reading should be feasible")
    sw_hour = [r for r in mp["routes"] if r["lo"] == Fraction(1)]
    rd_hour = [r for r in mp["readings"] if r["t"] == Fraction(1)]
    if not sw_hour or not rd_hour:
        raise ValueError("micro theme broken: switch and reading must share an hour")
    m_got = {(it["region"], it["day"]): (it["lower"], it["upper"]) for it in m_items}
    m_expect = {("east", 0): ("1", "2"), ("east", 1): ("0", "1"),
                ("west", 0): ("0", "1"), ("west", 1): ("0", "1")}
    if m_got != m_expect:
        raise ValueError("micro switch_reading bounds mismatch: %s" % m_got)
    for it in m_items:
        check_item_full(mp, m_model, m_obj, it)
    enum, cnt = enumerate_integer(mp, m_segs, m_obj)
    for key, (lo, hi) in enum.items():
        if (str(lo), str(hi)) != m_expect[key]:
            raise ValueError("micro enumeration mismatch at %s" % (key,))
    summary["micro"].append({"case": "switch_and_reading_same_hour",
                             "integer_trajectories": cnt,
                             "bounds": {"%s/%d" % k: v for k, v in m_got.items()}})

    up = load(tmp / "unique")
    u_model, u_segs, u_obj, u_items, u_farkas = assess(up)
    if u_items is None:
        raise ValueError("micro unique should be feasible")
    for it in u_items:
        if it["lower"] != it["upper"]:
            raise ValueError("micro unique: bound not unique")
        check_item_full(up, u_model, u_obj, it)
    enum, cnt = enumerate_integer(up, u_segs, u_obj)
    for key, (lo, hi) in enum.items():
        if lo != hi:
            raise ValueError("micro unique enumeration not unique")
    summary["micro"].append({"case": "all_daily_usage_unique",
                             "integer_trajectories": cnt,
                             "unique_bounds": len(u_items)})

    print(json.dumps({"verify": "ok", **summary}, ensure_ascii=False, indent=1))
    print("微型例连续可行域的极值由整数轨迹覆盖：有理证书已对整个连续可行集有效，"
          "而完整枚举的整数轨迹达到同一界；这不声称一般材料只有整数解。")
    return 0


# ---------------------------------------------------------------- 演示

def demo(material_dir, output_dir):
    try:
        pack = load(material_dir)
    except Exception as exc:
        print("读取材料失败: %s" % exc, file=sys.stderr)
        return 1
    model, segs, objectives, items, farkas = assess(pack)
    if items is None:
        materials = conflict_materials(model, farkas)
        result = {"status": "inconsistent", "model": model,
                  "conflict": {"weights": {k: str(v) for k, v in farkas.items()},
                               "materials": materials}}
        cert.check_conflict(model["rows"], result["conflict"]["weights"])
        write_outputs(output_dir, result,
                      render_conflict_report(pack, model, farkas, materials))
        print("材料联合矛盾（每表单独均可解释）。Farkas 冲突乘子：")
        print(json.dumps(result["conflict"]["weights"], ensure_ascii=False, indent=1))
        print("涉及材料: %s" % "、".join(materials))
        print("certificates.check_conflict: True")
        return 0
    tightened = tightened_bounds(pack, model, objectives, items)
    result = {"status": "feasible", "model": model, "bounds": items}
    write_outputs(output_dir, result,
                  render_report(pack, model, segs, items, tightened))
    print("真实核定区间（m3）：")
    for it in items:
        print("  %s 日%d: [%s, %s]" % (it["region"], it["day"], it["lower"], it["upper"]))
    it = items[0]
    w = it["upper_witness"]
    print("\n端点轨迹演示：%s 日%d 上界 %s 的完整联合见证" % (it["region"], it["day"],
                                                             it["upper"]))
    print("  期初累计: %s" % json.dumps(w["initial"], ensure_ascii=False))
    for seg in w["segments"]:
        print("  段 [%s,%s): %s" % (seg["start"], seg["end"],
                                    json.dumps(seg["volumes"], ensure_ascii=False)))
    ok1 = cert.check_point(model["rows"], w["point"])
    ok2 = cert.check_bound(model["rows"], it["upper_proof"]["objective"],
                           it["upper_proof"]["upper"], it["upper_proof"]["weights"])
    ok3 = check_trajectory(pack, model, w, it["region"], it["day"], it["upper"])
    print("证书与轨迹核验: check_point=%s check_bound=%s 物理轨迹核验=%s"
          % (ok1, ok2, ok3))
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 4 or sys.argv[1] not in ("run", "verify", "demo"):
        print(__doc__)
        sys.exit(2)
    cmd, mat, outd = sys.argv[1], sys.argv[2], sys.argv[3]
    try:
        code = {"run": run, "verify": verify, "demo": demo}[cmd](mat, outd)
    except Exception as exc:
        import traceback
        traceback.print_exc()
        print("失败: %s" % exc, file=sys.stderr)
        sys.exit(1)
    sys.exit(code)
