"""Hash every model output, so a refactor can be shown to move nothing.

Writes one hash per output (per dataclass field, per trace key) to a JSON file. Comparing two
files lists outputs that changed and outputs that exist on one side only (a deleted trace or
field shows up there, not as a change). Writes nothing to the repo.

usage (from the repo root, inside t700-dev):
    PYTHONPATH=src python3 tools/model_digest.py <out.json>
    python3 tools/model_digest.py --compare <before.json> <after.json>
"""

import dataclasses
import enum
import hashlib
import inspect
import json
import sys

import numpy as np


def compare(a_path, b_path):
    a, b = json.load(open(a_path)), json.load(open(b_path))
    changed = sorted(k for k in a.keys() & b.keys() if a[k] != b[k])
    only_a, only_b = sorted(a.keys() - b.keys()), sorted(b.keys() - a.keys())
    print(f"{len(a)} before, {len(b)} after, {len(a.keys() & b.keys())} shared")
    print(f"CHANGED ({len(changed)}):", *changed[:60], sep="\n  ")
    print(
        f"ONLY BEFORE ({len(only_a)}):", *sorted({k.rsplit(" ", 1)[-1] for k in only_a}), sep="\n  "
    )
    print(
        f"ONLY AFTER ({len(only_b)}):", *sorted({k.rsplit(" ", 1)[-1] for k in only_b}), sep="\n  "
    )
    print("VERDICT:", "IDENTICAL on every shared output" if not changed else "OUTPUTS MOVED")
    return 1 if changed else 0


if sys.argv[1] == "--compare":
    sys.exit(compare(sys.argv[2], sys.argv[3]))

from t700 import constants as c  # noqa: E402
from t700 import maps, realtime, trim  # noqa: E402
from t700.control import loop, schedules  # noqa: E402
from t700.engine import Ambient, frame  # noqa: E402
from t700.linear import DOF, extract  # noqa: E402

OUT: dict[str, str] = {}


def enc(v, depth=0):
    """A canonical, exact text form of a value."""
    if depth > 6:
        return "<deep>"
    if isinstance(v, np.ndarray):
        if v.dtype.kind in "fc":
            return f"{v.dtype}{v.shape}:" + ",".join(
                (x.hex() if isinstance(x, float) else f"{x.real.hex()}j{x.imag.hex()}")
                for x in v.ravel().tolist()
            )
        return f"{v.dtype}{v.shape}:{v.ravel().tolist()}"
    if isinstance(v, (bool, np.bool_)):
        return repr(bool(v))
    if isinstance(v, (float, np.floating)):
        return float(v).hex()
    if isinstance(v, (int, np.integer)):
        return repr(int(v))
    if isinstance(v, complex):
        return f"{v.real.hex()}j{v.imag.hex()}"
    if isinstance(v, enum.Enum):
        return repr(v.value)
    if isinstance(v, str) or v is None:
        return repr(v)
    if isinstance(v, dict):
        return (
            "{"
            + ",".join(
                f"{k!r}:{enc(x, depth + 1)}"
                for k, x in sorted(v.items(), key=lambda kv: repr(kv[0]))
            )
            + "}"
        )
    if isinstance(v, (list, tuple)):
        return "[" + ",".join(enc(x, depth + 1) for x in v) + "]"
    if callable(v):
        return "<callable>"
    return type(v).__name__


def put(key, v):
    OUT[key] = hashlib.sha256(enc(v).encode()).hexdigest()[:20]


def fields(obj):
    if dataclasses.is_dataclass(obj):
        return {f.name: getattr(obj, f.name) for f in dataclasses.fields(obj)}
    if hasattr(obj, "_asdict"):
        return obj._asdict()
    d = dict(getattr(obj, "__dict__", {}))
    for s in getattr(type(obj), "__slots__", ()):
        if hasattr(obj, s):
            d[s] = getattr(obj, s)
    return d


def put_obj(prefix, obj):
    """One hash per field, so a deleted field reads as 'only before', not as a change."""
    fs = fields(obj)
    if not fs:
        put(prefix, obj)
    for name, v in fs.items():
        if dataclasses.is_dataclass(v) or hasattr(v, "_asdict"):
            put_obj(f"{prefix}.{name}", v)
        else:
            put(f"{prefix} {name}", v)


AMB = Ambient(14.696, 518.67)
HOT = Ambient(14.0, 560.0)
WF_PPH = {1: 476.3, 2: 349.3, 3: 267.7}
DQ = {1: 0.019471, 2: 0.015869, 3: 0.012950}


def pps(pph):
    return pph / 3600.0


# ---- every table, after physical conditioning
for mod in (maps, schedules):
    for name, fn in sorted(vars(mod).items()):
        if not callable(fn) or name.startswith("_") or inspect.isclass(fn):
            continue
        try:
            sig = inspect.signature(fn)
        except (TypeError, ValueError):
            continue
        if any(
            p.default is p.empty and p.kind in (p.POSITIONAL_ONLY, p.POSITIONAL_OR_KEYWORD)
            for p in sig.parameters.values()
        ):
            continue
        try:
            obj = fn()
        except Exception as e:  # noqa: BLE001
            put(f"table {mod.__name__}.{name} raised", type(e).__name__)
            continue
        if not isinstance(obj, (maps.Curve, maps.SpeedMap)):
            continue
        put_obj(f"table {mod.__name__}.{name}", obj)
        for probe in (0.0, 0.5, 1.0, 50.0, 90.0, 100.0, 1000.0):
            try:
                put(f"table {mod.__name__}.{name} at {probe}", obj(probe))
            except Exception:  # noqa: BLE001
                pass
put("maps.CONDITIONING", maps.CONDITIONING)

# ---- trims: the printed points at both speeds, a hot day, and the whole continuation sweep
for t, wf in WF_PPH.items():
    for np_rpm in (20895.0, c.NP_DES):
        r = trim.solve(pps(wf), np_rpm, AMB)
        put_obj(f"trim t{t} np{np_rpm}", r)
    put_obj(f"trim t{t} hot", trim.solve(pps(wf), c.NP_DES, HOT))
for i, r in enumerate(trim.sweep(np.arange(125.0, 801.0, 25.0), c.NP_DES, AMB)):
    put_obj(f"sweep {125 + 25 * i}", r)

# ---- linear models: every variant, both steps, with and without the rotor load
for t, wf in WF_PPH.items():
    r = trim.solve(pps(wf), 20895.0, AMB)
    for dof in DOF:
        for rel in (0.02, 1e-5):
            m = extract(r, pps(wf), dof, AMB, j_load=c.J_LOAD_UH60A, dq_req_dnp=DQ[t], rel_step=rel)
            put_obj(f"linear t{t} {dof.value} rel{rel}", m)
            put(f"linear t{t} {dof.value} rel{rel} eig", np.sort_complex(np.linalg.eigvals(m.A)))
        put_obj(f"linear t{t} {dof.value} bare", extract(r, pps(wf), dof, AMB))

# ---- open-loop transients: Figures 9 and 10 and their switches
wf0 = pps(400.0)
r0 = trim.solve(wf0, c.NP_DES, AMB)
f0 = frame(r0.state, wf0, AMB)
put_obj("frame400", f0)
for fig, hi, t_step in ((9, 775.0, 0.5232), (10, 125.0, 0.5217)):
    for hs in (True, False):
        for lag in (True, False):
            for multi in (True, False):
                st = realtime.from_trim(r0, lag_whole_flow=lag)
                put_obj(f"rtseed fig{fig} lag{lag}", st)
                tr = realtime.run(
                    st,
                    lambda t, h=hi / 3600.0, s=t_step: wf0 if t < s else h,
                    AMB,
                    duration_s=3.0,
                    dt=0.007,
                    multirate=multi,
                    q_req_ftlbf=f0.q_pt_ftlbf,
                    integrate_np=False,
                    heat_sink=hs,
                    lag_whole_flow=lag,
                )
                for k, v in tr.items():
                    put(f"rt fig{fig} hs{hs} lag{lag} multi{multi} {k}", v)
# NP integrated, rotor inertia and a load torque, default frame
st = realtime.from_trim(r0)
tr = realtime.run(
    st,
    lambda t: wf0 if t < 0.3 else pps(500.0),
    AMB,
    duration_s=2.0,
    q_req_ftlbf=f0.q_pt_ftlbf,
    j_load=c.J_LOAD_UH60A,
)
for k, v in tr.items():
    put(f"rt np-integrated {k}", v)

# ---- closed loop: settle at each trim, then a collective step and a load step
for t, wf in WF_PPH.items():
    r = trim.solve(pps(wf), c.NP_DES, AMB)
    f = frame(r.state, pps(wf), AMB)
    xcpc = {1: 52.75, 2: 35.87, 3: 23.01}[t]
    pilot = loop.Pilot(xcpc_pct=xcpc, pas_deg=100.0, pcprf_pct=20895.0 * 100.0 / c.NP_DES)
    s = loop.seed(r, pilot, AMB)
    put_obj(f"loop t{t} seed", s)
    q0 = f.q_pt_ftlbf
    for n in range(1500):
        if n == 600:
            pilot = loop.Pilot(xcpc_pct=xcpc + 8.0, pas_deg=110.0, pcprf_pct=100.0)
        scale = 1.15 if n >= 1000 else 1.0
        s, e, h, fr = loop.step(
            s,
            pilot,
            AMB,
            dt=0.007,
            load=lambda v, q=q0, k=scale: k * q + 0.2 * (v - 20895.0),
            j_load=c.J_LOAD_UH60A,
            heat_sink=True,
        )
        if n % 50 == 49:
            put_obj(f"loop t{t} n{n} state", s)
            put_obj(f"loop t{t} n{n} ecu", e)
            put_obj(f"loop t{t} n{n} hmu", h)
            put_obj(f"loop t{t} n{n} frame", fr)

json.dump(OUT, open(sys.argv[1], "w"), indent=0, sort_keys=True)
print(sys.argv[1], len(OUT), "outputs hashed")
