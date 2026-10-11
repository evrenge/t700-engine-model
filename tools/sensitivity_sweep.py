"""What is each map's read error worth? Push +1 sigma through the model and measure.

Every digitized CSV carries a per-point 1 sigma in its own units. This shifts one map by
exactly that, one axis at a time, resolves all three trims and both open-loop transients,
and reports what moved. It is the only statement this project makes about what the
digitizing uncertainty is *worth*, as opposed to what it is.

## Why it runs in a throwaway worktree, one subprocess per case

Neither half of that is caution. `t700.maps` caches every loaded `Curve` and freezes its
arrays -- deliberately, after `maps.f2().y *= 1.02` silently corrupted a table for a whole
process -- so a perturbation applied in memory either raises or, worse, does nothing while
appearing to work. A first attempt at this sweep reported every sensitivity as exactly
0.000 % for that reason. And `DATA_ROOT` is resolved at import, so pointing the model at a
different data directory means a different checkout, not a different variable.

So: `git worktree add` a scratch copy of HEAD, perturb a CSV in it, run one subprocess,
restore the file, repeat. The real tree is never written to.

## Reading the output

The steady column is the worst shaft-power change over the three trims. The transient column
is the largest change in any whole-curve panel rms, in points of that panel's excursion --
and it is the column that stops the table lying about dynamic-only quantities. `f_hs` sets a
heat-sink time constant that does no work at equilibrium, so it costs exactly 0.000 % at
every trim however wrong it is, and shows up only in the transient.

    python3 tools/sensitivity_sweep.py                 # every map, both axes
    python3 tools/sensitivity_sweep.py --keep          # leave the worktree for inspection
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent

CASES = [
    ("f2", "maps/f2_compressor_temperature.csv", "sx", "sy"),
    ("f3", "maps/f3_seal_bleed_fraction.csv", "sigma_x", "sigma_y"),
    ("f4", "maps/f4_pt_balance_bleed_fraction.csv", "sigma_x", "sigma_y"),
    ("f5", "maps/f5_tip_leak_cooling_bleed_fraction.csv", "sigma_x", "sigma_y"),
    ("f6", "maps/f6_combustor_efficiency.csv", "sigma_x", "sigma_y"),
    ("f8", "maps/f8_pt_energy.csv", "sigma_x", "sigma_y"),
    ("f9", "maps/f9_pt_mass_flow.csv", "sx", "sy"),
    ("f10", "maps/f10_exhaust_pressure_loss.csv", "sigma_x", "sigma_y"),
    ("f_hs", "maps/fhs_heat_sink_constant.csv", "sigma_x", "sigma_y"),
]

PROBE = """
import json, sys
sys.path.insert(0, "src"); sys.path.insert(0, "validation")
from t700 import trim
from t700.units import wf_pps_from_pph
from t700.engine import Ambient
import report as R
AMB = Ambient(14.696, 518.67)
out = {"trims": {}}
for k, w in (("hover", 476.3), ("level", 349.3), ("descent", 267.7)):
    r = trim.solve(wf_pps_from_pph(w), 20895.0, AMB)
    out["trims"][k] = {"shp": float(r.shp), "ngc": float(r.frame.ngc_pct)}
out["curve"] = {k: float(v) for k, v in R.whole_curve().items() if isinstance(v, (int, float))}
b1 = R.table_b1()
out["b1"] = float(b1["rms"])
print(json.dumps(out))
"""


def _run(sandbox: Path):
    p = subprocess.run(
        [sys.executable, "sweep_probe.py"], cwd=sandbox, capture_output=True, text=True
    )
    if p.returncode != 0:
        tail = p.stderr.strip().splitlines()[-1] if p.stderr.strip() else "no output"
        return None, tail
    return json.loads(p.stdout), None


def _perturb(sandbox: Path, rel: str, col: str, sig_col: str) -> bool:
    """Rewrite one column shifted by +1 of its own per-point sigma. False if it cannot be."""
    full = sandbox / "data" / rel
    lines = full.read_text().splitlines()
    head = next(i for i, ln in enumerate(lines) if ln.strip() and not ln.startswith("#"))
    names = [c.strip() for c in lines[head].split(",")]
    if col not in names or sig_col not in names:
        return False
    ci, si = names.index(col), names.index(sig_col)
    body = [i for i in range(head + 1, len(lines)) if lines[i].strip()]
    vals = np.array([float(lines[i].split(",")[ci]) for i in body])
    sig = np.array([float(lines[i].split(",")[si]) for i in body])
    new = vals + sig
    if col == "x" and not np.all(np.diff(new) > 0):
        return False  # a shifted abscissa that is no longer increasing is not a map
    for k, i in enumerate(body):
        cells = [c.strip() for c in lines[i].split(",")]
        cells[ci] = f"{new[k]:.10g}"
        lines[i] = ",".join(cells)
    full.write_text("\n".join(lines) + "\n")
    return True


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--keep", action="store_true", help="leave the scratch worktree behind")
    args = ap.parse_args(argv)

    tmp = Path(tempfile.mkdtemp(prefix="t700-sweep-"))
    sandbox = tmp / "wt"
    subprocess.run(
        ["git", "worktree", "add", "-f", str(sandbox), "HEAD"],
        cwd=ROOT,
        check=True,
        capture_output=True,
    )
    (sandbox / "sweep_probe.py").write_text(PROBE)
    try:
        base, err = _run(sandbox)
        if base is None:
            print(f"the unperturbed model does not run in the worktree: {err}")
            return 1
        print("+1 SIGMA OF EACH MAP'S OWN READ ERROR, PUSHED THROUGH THE MODEL")
        shp0 = ", ".join(f"{k} {v['shp']:.2f}" for k, v in base["trims"].items())
        print(f"baseline SHP  {shp0}")
        print(
            f"\n{'map':>6} {'axis':>5} {'hover':>9} {'level':>9} {'descent':>9} "
            f"{'worst':>8} | {'worst transient panel':>28} | {'B.1 rms':>8}"
        )
        print("-" * 96)
        for name, rel, sx, sy in CASES:
            orig = (sandbox / "data" / rel).read_text()
            for axis, col, sig_col in (("x", "x", sx), ("y", "y", sy)):
                if not _perturb(sandbox, rel, col, sig_col):
                    print(
                        f"{name:>6} {axis:>5}   not perturbable (column missing, or the "
                        f"shift breaks monotonicity)"
                    )
                    (sandbox / "data" / rel).write_text(orig)
                    continue
                got, err = _run(sandbox)
                (sandbox / "data" / rel).write_text(orig)
                if got is None:
                    print(f"{name:>6} {axis:>5}   model failed: {err}")
                    continue
                shp = [
                    (got["trims"][k]["shp"] / base["trims"][k]["shp"] - 1) * 100
                    for k in ("hover", "level", "descent")
                ]
                dcur = {k: got["curve"][k] - base["curve"][k] for k in base["curve"]}
                ck = max(dcur, key=lambda k: abs(dcur[k]))
                print(
                    f"{name:>6} {axis:>5} "
                    + " ".join(f"{v:+9.3f}" for v in shp)
                    + f" {max(abs(v) for v in shp):8.3f} | {ck:>18} {dcur[ck]:+8.4f} "
                    f"| {got['b1'] - base['b1']:+8.4f}"
                )
    finally:
        if not args.keep:
            subprocess.run(
                ["git", "worktree", "remove", "--force", str(sandbox)],
                cwd=ROOT,
                capture_output=True,
            )
        else:
            print(f"\nworktree left at {sandbox}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
