"""C-3: 파라메트릭 스터디 — 형상 비 (H/p, t/p, R_c/t) 와 쉐브론 각 β 가 등가 강성·등가 두께·응력 집중·단위 응답에 미치는 영향.

기준: 05 §5 가정 기본값 (원호+접선 p=9, H=3.2, t=0.6, R_c=R_v=1.5 mm; 316L E=193 GPa, ν=0.3; β=60°).
한 번에 한 변수만 바꾼다 (one-at-a-time). 출력:
    docs/param_study/sweep_<var>.csv     각 케이스의 형상량·채택식 A,D·등가 두께·K_t·단위 응답
    docs/fig/param_<var>.png             경향 그림 (소형 다중 패널, 단일 축, 계열 ≤3, 직접 라벨)
    docs/11_parametric_study.md 에 요약 표·해석 (수동 작성, 이 스크립트 출력 인용)
실행: python calc/parametric.py
"""
from __future__ import annotations

import csv
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from geometry import arc_tangent                       # noqa: E402
from stiffness import Sheet, adopted                    # noqa: E402
from transform import zone_stiffness, resultants_to_local  # noqa: E402
from equivalent_plate import equivalent_thicknesses     # noqa: E402
from stress_recovery import recovery_from_resultants, kt_x_tension  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_CSV = os.path.join(HERE, "..", "docs", "param_study")
OUT_FIG = os.path.join(HERE, "..", "docs", "fig")
M = Sheet(193000.0, 0.3)
BASE = dict(p=9.0, H=3.2, t=0.6, R_c=1.5, beta=60.0)

UNIT_CASES = (("N_x", [1, 0, 0], [0, 0, 0]), ("N_y", [0, 1, 0], [0, 0, 0]), ("N_xy", [0, 0, 1], [0, 0, 0]),
              ("M_x", [0, 0, 0], [1, 0, 0]), ("M_y", [0, 0, 0], [0, 1, 0]), ("M_xy", [0, 0, 0], [0, 0, 1]))


def evaluate(p: float, H: float, t: float, R_c: float, beta: float) -> dict:
    """한 케이스의 모든 출력량."""
    g = arc_tangent(p, H, t, R_c=R_c)
    K = adopted(g, M)
    th = equivalent_thicknesses(K, g)
    row = dict(p=p, H=H, t=t, R_c=R_c, beta=beta, H_over_p=H / p, t_over_p=t / p, Rc_over_t=R_c / t, f_over_t=H / 2 / t,
               alpha_deg=g.extra["alpha_deg"], T_L=g.extra["T_L"], lam=g.l / g.c, J1=g.avg_inv_sqrt_a, J2=g.eps2_avg_phi2_sqrt_a,
               t_over_Rmin=t / g.extra["R_min"])
    row.update({k: K[k] for k in ("A11", "A12", "A22", "A66", "D11", "D12", "D22", "D66")})
    row.update(A22_over_A11=K["A22"] / K["A11"], A66_over_A11=K["A66"] / K["A11"], D22_over_D11=K["D22"] / K["D11"],
               t_b=th["t_b"], t_b_over_t=th["t_b"] / t, t_s=th["t_s"], Kt=kt_x_tension(g.f, g.h))
    # 등가 판 강성을 원판 평판(두께 t)에 대해 정규화: 능선 방향 A22/(E t) 등
    row.update(A11_over_Et=K["A11"] / (M.E * t), A22_over_Et=K["A22"] / (M.E * t), A66_over_Gt=K["A66"] / (M.G * t),
               D11_over_Dflat=K["D11"] / (M.E * t**3 / 12 / (1 - M.nu**2)), D22_over_Dflat=K["D22"] / (M.E * t**3 / 12 / (1 - M.nu**2)),
               D66_over_Gt3=K["D66"] / (M.G * t**3 / 12))
    for name, N, Mm in UNIT_CASES:                        # 국부축 단위 응답 [MPa per unit]
        r = recovery_from_resultants(g, M, N, Mm, n=400)
        row[f"vm_per_{name}"] = r.max_vm()[0]
    # 판축 응답: 우측 영역에 판축 단위 합력 (N_xp, N_yp, M_xp, M_yp) → 국부 변환 → 최대 vM
    for name, N, Mm in (("N_xp", [1, 0, 0], [0, 0, 0]), ("N_yp", [0, 1, 0], [0, 0, 0]), ("M_xp", [0, 0, 0], [1, 0, 0]), ("M_yp", [0, 0, 0], [0, 1, 0])):
        Nl, Ml = resultants_to_local(np.array(N, float), beta, "right"), resultants_to_local(np.array(Mm, float), beta, "right")
        row[f"vm_per_{name}"] = recovery_from_resultants(g, M, Nl, Ml, n=400).max_vm()[0]
    Z = zone_stiffness(K, beta, "right")
    row.update({f"Z{k}": Z[k] for k in ("A11", "A12", "A22", "A66", "A16", "A26", "D11", "D12", "D22", "D66", "D16", "D26")})
    row.update(coupA16=Z["A16"] / math.sqrt(Z["A11"] * Z["A66"]), coupD16=Z["D16"] / math.sqrt(Z["D11"] * Z["D66"]),
               ZA22_over_ZA11=Z["A22"] / Z["A11"], ZD22_over_ZD11=Z["D22"] / Z["D11"])
    return row


def sweep(var: str, values) -> list[dict]:
    rows = []
    for v in values:
        kw = dict(BASE); kw[var] = v
        if var == "H":                      # H 가 작아지면 원호가 겹치지 않도록 R_c 를 H 에 비례해 줄임 (기본 R_c/H = 0.469 유지)
            kw["R_c"] = min(BASE["R_c"], 0.469 * v)
        try:
            rows.append(evaluate(kw["p"], kw["H"], kw["t"], kw["R_c"], kw["beta"]))
        except ValueError as e:
            print(f"  skip {var}={v}: {e}")
    return rows


def write_csv(name: str, rows: list[dict]) -> None:
    os.makedirs(OUT_CSV, exist_ok=True)
    with open(os.path.join(OUT_CSV, f"sweep_{name}.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)


# ----------------------------------------------------------------------------- 그림
C = ["#2a78d6", "#eb6834", "#1baf7a"]          # 검증된 범주형 3색 (dataviz 참조 팔레트, 고정 순서)
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e6e6e3"


def _style(ax, xlabel, ylabel, title, log=False):
    import matplotlib.pyplot as plt  # noqa
    ax.set_title(title, loc="left", fontsize=10.5, color=INK, fontweight="bold", pad=8)
    ax.set_xlabel(xlabel, color=INK2, fontsize=9); ax.set_ylabel(ylabel, color=INK2, fontsize=9)
    ax.grid(True, axis="y", color=GRID, lw=0.8); ax.set_axisbelow(True)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"): ax.spines[sp].set_color("#bdbdb8")
    ax.tick_params(colors=INK2, labelsize=8.5)
    if log: ax.set_yscale("log")


def _lines(ax, x, series, base_x=None, log=False):
    """series: [(label, y)], ≤3. 선 끝 직접 라벨(겹침 회피), 기준값 위치 마커."""
    ends = []
    for i, (lab, y) in enumerate(series):
        ax.plot(x, y, color=C[i], lw=2, solid_capstyle="round")
        ends.append((lab, float(y[-1])))
        if base_x is not None:
            yb = np.interp(base_x, x, y)
            ax.plot(base_x, yb, "o", color=C[i], ms=6, mec="white", mew=1.5)
    if base_x is not None:
        ax.axvline(base_x, color="#bdbdb8", lw=0.8, ls=":")
    ax.margins(x=0.02); ax.set_xlim(x[0], x[-1] + (x[-1] - x[0]) * 0.24)
    if log: ax.set_yscale("log")
    # 라벨 겹침 회피: 표시 좌표(픽셀)에서 최소 간격 확보
    ax.figure.canvas.draw()
    tr, inv = ax.transData, ax.transData.inverted()
    pix = [tr.transform((x[-1], v))[1] for _, v in ends]
    order = np.argsort(pix); adj = list(pix); gap = 11.0
    for k in range(1, len(order)):
        lo, hi = order[k - 1], order[k]
        if adj[hi] - adj[lo] < gap: adj[hi] = adj[lo] + gap
    for (lab, v), py in zip(ends, adj):
        yl = inv.transform((tr.transform((x[-1], v))[0], py))[1]
        ax.annotate(lab, (x[-1], yl), (5, 0), textcoords="offset points", va="center", fontsize=8.5, color=INK)


def fig_sweep(name: str, rows: list[dict], xkey: str, xlabel: str, base_x: float, panels: list) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9})
    x = np.array([r[xkey] for r in rows])
    n = len(panels)
    fig, axs = plt.subplots(1, n, figsize=(4.1 * n, 3.6))
    axs = np.atleast_1d(axs)
    for ax, (title, ylabel, series, log) in zip(axs, panels):
        _style(ax, xlabel, ylabel, title, log)
        _lines(ax, x, [(lab, np.array([r[k] for r in rows])) for lab, k in series], base_x, log)
    fig.suptitle(f"Parametric study — {name} sweep (others at default; marker = default value)", x=0.01, ha="left", fontsize=10, color=INK2, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    os.makedirs(OUT_FIG, exist_ok=True)
    fig.savefig(os.path.join(OUT_FIG, f"param_{name}.png"), dpi=140); plt.close(fig)


def main() -> None:
    sweeps = {
        "Hp": ("H", np.round(np.linspace(1.4, 5.4, 21), 3), "H_over_p", "H/p  (depth / pitch)", BASE["H"] / BASE["p"]),
        "tp": ("t", np.round(np.linspace(0.3, 1.2, 19), 3), "t_over_p", "t/p  (sheet thickness / pitch)", BASE["t"] / BASE["p"]),
        "Rt": ("R_c", np.round(np.linspace(0.6, 2.6, 21), 3), "Rc_over_t", "R_c/t  (ridge radius / thickness)", BASE["R_c"] / BASE["t"]),
        "beta": ("beta", np.round(np.linspace(20, 80, 25), 3), "beta", "β  [deg]  (ridge to flow axis)", BASE["beta"]),
    }
    results = {}
    for name, (var, vals, xkey, xlabel, bx) in sweeps.items():
        rows = sweep(var, vals); results[name] = rows; write_csv(name, rows)
        print(f"{name}: {len(rows)} cases")
    # ---- 그림 정의 (패널: 제목, y라벨, [(라벨, 키)], 로그)
    fig_sweep("Hp", results["Hp"], "H_over_p", "H/p", BASE["H"] / BASE["p"], [
        ("Stiffness normalised by flat sheet", "ratio (log)", [("A22 / E t", "A22_over_Et"), ("A11 / E t", "A11_over_Et"), ("D22 / D_flat", "D22_over_Dflat")], True),
        ("Orthotropy and equivalent thickness", "ratio (log)", [("A22 / A11", "A22_over_A11"), ("t_b / t", "t_b_over_t")], True),
        ("Ridge stress concentration K_t = 1+6f/t", "K_t", [("K_t", "Kt")], False),
        ("Max von Mises per unit resultant", "MPa per (N/mm) or (N·mm/mm)", [("N_x", "vm_per_N_x"), ("M_y", "vm_per_M_y"), ("M_xy", "vm_per_M_xy")], True),
    ])
    fig_sweep("tp", results["tp"], "t_over_p", "t/p", BASE["t"] / BASE["p"], [
        ("Stiffness normalised by flat sheet", "ratio (log)", [("A22 / E t", "A22_over_Et"), ("A11 / E t", "A11_over_Et"), ("D22 / D_flat", "D22_over_Dflat")], True),
        ("Equivalent thickness", "mm", [("t_b", "t_b"), ("t_s = hλ", "t_s"), ("t", "t")], False),
        ("Ridge stress concentration K_t", "K_t", [("K_t", "Kt")], False),
        ("Max von Mises per unit resultant", "MPa per unit (log)", [("N_x", "vm_per_N_x"), ("M_y", "vm_per_M_y"), ("M_xy", "vm_per_M_xy")], True),
    ])
    fig_sweep("Rt", results["Rt"], "Rc_over_t", "R_c/t", BASE["R_c"] / BASE["t"], [
        ("Section shape", "value", [("α [deg]", "alpha_deg"), ("T_L [mm]", "T_L")], False),
        ("Stiffness normalised by flat sheet", "ratio (log)", [("A11 / E t", "A11_over_Et"), ("D22 / D_flat", "D22_over_Dflat"), ("λ = l/c", "lam")], True),
        ("Shell-theory validity  t / R_min", "t / R_min", [("t/R_min", "t_over_Rmin")], False),
        ("Max von Mises per unit resultant", "MPa per unit (log)", [("N_x", "vm_per_N_x"), ("M_y", "vm_per_M_y"), ("M_xy", "vm_per_M_xy")], True),
    ])
    fig_sweep("beta", results["beta"], "beta", "β [deg]", BASE["beta"], [
        ("Right-zone plate-axis membrane stiffness", "N/mm (log)", [("Ā11 (x_p)", "ZA11"), ("Ā22 (y_p)", "ZA22"), ("Ā66", "ZA66")], True),
        ("Right-zone plate-axis bending stiffness", "N·mm (log)", [("D̄11 (x_p)", "ZD11"), ("D̄22 (y_p)", "ZD22"), ("D̄66", "ZD66")], True),
        ("Coupling ratio  X16 / √(X11·X66)  (<1)", "ratio", [("A16", "coupA16"), ("D16", "coupD16")], False),
        ("Max von Mises per unit plate-axis load", "MPa per unit (log)", [("N_yp", "vm_per_N_yp"), ("M_yp", "vm_per_M_yp"), ("M_xp", "vm_per_M_xp")], True),
    ])
    print("figures written to", OUT_FIG)
    # 요약 출력
    def pick(rows, key, xkey, xs):
        return [np.interp(x, [r[xkey] for r in rows], [r[key] for r in rows]) for x in xs]
    print("\n[H/p] 0.2 / 0.356(base) / 0.5 :")
    for key in ("A22_over_A11", "Kt", "t_b_over_t", "vm_per_N_x", "vm_per_M_y", "vm_per_M_xy"):
        print(f"  {key:14s}", [f"{v:.3g}" for v in pick(results["Hp"], key, "H_over_p", [0.2, BASE["H"] / BASE["p"], 0.5])])
    print("[t/p] 0.044 / 0.0667(base) / 0.1 :")
    for key in ("A22_over_A11", "Kt", "t_b", "vm_per_N_x", "vm_per_M_y", "vm_per_M_xy"):
        print(f"  {key:14s}", [f"{v:.3g}" for v in pick(results["tp"], key, "t_over_p", [0.044, BASE["t"] / BASE["p"], 0.1])])
    print("[R_c/t] 1.5 / 2.5(base) / 4 :")
    for key in ("alpha_deg", "A11_over_Et", "D22_over_Dflat", "t_over_Rmin", "vm_per_N_x", "vm_per_M_xy"):
        print(f"  {key:14s}", [f"{v:.3g}" for v in pick(results["Rt"], key, "Rc_over_t", [1.5, 2.5, 4.0])])
    print("[β] 30 / 45 / 60(base) / 75 :")
    for key in ("ZA11", "ZA22", "ZA66", "ZD11", "ZD22", "ZD66", "coupA16", "coupD16", "vm_per_N_yp", "vm_per_M_yp", "vm_per_M_xp"):
        print(f"  {key:14s}", [f"{v:.4g}" for v in pick(results["beta"], key, "beta", [30, 45, 60, 75])])


if __name__ == "__main__":
    main()
