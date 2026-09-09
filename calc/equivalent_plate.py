"""등가 판 두께 정의 및 ANSYS 등가 판 입력(preintegrated general shell section) 생성 — B-4 구현.

docs/09_equivalent_thickness.md 의 정의를 구현한다.

등가 두께 정의 (국부축 채택식 A, D 로부터):
    t_b  = sqrt(12 D11 / A11)   굽힘-막 정합 두께 (11 방향).  채택식에서는 12 방향, 22 방향도 같은 값 → 단일 t_b
    t_s  = sqrt(12 D66 / A66)   전단-비틀림 정합 두께  (= h λ, t_b 와 크게 다름 → 단일 균질 판으로는 A66/D66 동시 만족 불가)
    t_m  = h λ                  질량(면적) 등가 두께
FE 입력 방식:
    (i)  preintegrated general shell section (SECTYPE,,GENS + SSPA/SSPB/SSPD/SSPE/SSPM) — A, D 를 성분 그대로 입력. 채택.
    (ii) 균질 직교이방성 단일층 (t_b + 공학상수) — A66 또는 D66 중 하나만 맞음. 국부축 비교·간이 검토용.
ANSYS 명령 인자 순서 (ansys-mapdl-core 0.74 명령 레퍼런스 docstring 으로 확인, 2026-09-09):
    SSPA, A11, A21, A31, A22, A32, A33 [,T]      SSPB, B11, B21, B31, B22, B32, B33 [,T, B12, B13, B23]
    SSPD, D11, D21, D31, D22, D32, D33 [,T]      SSPE, E11, E21, E22 [,T]      SSPM, DENS [,T]
    인덱스 3 = xy(전단/비틀림). 일반화 변형률은 {ε_x, ε_y, γ_xy, κ_x, κ_y, κ_xy(공학 비틀림)} 으로 가정 →
    E 단계에서 단일 요소 비틀림 시험으로 M_xy/κ_xy 규약을 확인한다 (docs/09 §5).
"""
from __future__ import annotations

import math
from typing import Dict, Optional

import numpy as np

from geometry import CorrugationGeometry, arc_tangent
from stiffness import Sheet, adopted
from transform import zone_stiffness


# ----------------------------------------------------------------------------- 등가 두께
def equivalent_thicknesses(K: Dict[str, float], g: CorrugationGeometry) -> Dict[str, float]:
    lam = g.l / g.c
    out = dict(
        t_b11=math.sqrt(12 * K["D11"] / K["A11"]),
        t_b22=math.sqrt(12 * K["D22"] / K["A22"]),
        t_b12=math.sqrt(12 * K["D12"] / K["A12"]) if K["A12"] > 0 else float("nan"),
        t_s=math.sqrt(12 * K["D66"] / K["A66"]),
        t_m=g.h * lam,
        t_sheet=g.h,
    )
    out["t_b"] = out["t_b11"]
    out["t_b_spread"] = max(out["t_b11"], out["t_b22"], out["t_b12"]) / min(out["t_b11"], out["t_b22"], out["t_b12"]) - 1
    return out


def homogeneous_orthotropic(K: Dict[str, float], t_e: float, shear_from: str = "D") -> Dict[str, float]:
    """방식 (ii): 두께 t_e 의 균질 직교이방성 판.  11/12/22 는 A 기준(= D 기준, t_e=t_b 일 때),
    전단 G 는 'A' → A66/t_e, 'D' → 12 D66/t_e³ 중 택일.  반환에 각 방식의 잔차 포함."""
    det = K["A11"] * K["A22"] - K["A12"] ** 2
    Ex, Ey = det / (K["A22"] * t_e), det / (K["A11"] * t_e)
    nu_xy, nu_yx = K["A12"] / K["A22"], K["A12"] / K["A11"]
    G_A, G_D = K["A66"] / t_e, 12 * K["D66"] / t_e**3
    G = G_D if shear_from == "D" else G_A
    # 재현 검사
    s = t_e**3 / 12
    resid = dict(
        D11=(Ex * s / (1 - nu_xy * nu_yx)) / K["D11"] - 1,
        D22=(Ey * s / (1 - nu_xy * nu_yx)) / K["D22"] - 1,
        A66=(G * t_e) / K["A66"] - 1,
        D66=(G * s) / K["D66"] - 1,
    )
    return dict(t_e=t_e, E_x=Ex, E_y=Ey, nu_xy=nu_xy, nu_yx=nu_yx, G_xy=G, G_from_A=G_A, G_from_D=G_D, resid=resid)


# ----------------------------------------------------------------------------- 횡전단·질량 (Mindlin 쉘 입력 보조량)
def transverse_shear_estimate(g: CorrugationGeometry, m: Sheet, factor: float = 1.0) -> Dict[str, float]:
    """SSPE 용 횡전단 강성 1차 추정 [N/mm].  Kirchhoff 등가모델에는 없는 량이므로 원판 전개면적 기준
    k·G·h·λ (k=5/6) 를 두 방향에 같이 두고, E 단계에서 ×0.1 / ×10 감도로 영향이 없음을 확인한다."""
    E = factor * (5 / 6) * m.G * g.h * (g.l / g.c)
    return dict(E11=E, E21=0.0, E22=E)


def mass_per_area(rho: float, g: CorrugationGeometry) -> float:
    """SSPM 입력: 단위 두께 가정 밀도 = 단위 투영면적당 질량 = ρ h λ.  (ρ tonne/mm³ → tonne/mm²)"""
    return rho * g.h * g.l / g.c


# ----------------------------------------------------------------------------- APDL 스니펫
def _fmt(v: float) -> str:
    return f"{v:.6e}"


def gens_snippet(K: Dict[str, float], secid: int, name: str, E_ts: Dict[str, float], dens_area: float,
                 comment: str = "") -> str:
    """하나의 등가 판 영역에 대한 preintegrated general shell section APDL 블록.
    K: 판축(또는 국부축) 강성 dict (A11,A12,A22,A66[,A16,A26], D..).  인덱스 3 = xy."""
    a16, a26 = K.get("A16", 0.0), K.get("A26", 0.0)
    d16, d26 = K.get("D16", 0.0), K.get("D26", 0.0)
    lines = [
        f"! --- {comment}" if comment else "!",
        f"SECTYPE,{secid},GENS,,{name}          ! preintegrated general shell section",
        f"SSPA,{_fmt(K['A11'])},{_fmt(K['A12'])},{_fmt(a16)},{_fmt(K['A22'])},{_fmt(a26)},{_fmt(K['A66'])}   ! A11,A21,A31,A22,A32,A33 [N/mm]",
        f"SSPB,0,0,0,0,0,0                      ! B = 0 (대칭 단면)",
        f"SSPD,{_fmt(K['D11'])},{_fmt(K['D12'])},{_fmt(d16)},{_fmt(K['D22'])},{_fmt(d26)},{_fmt(K['D66'])}   ! D11,D21,D31,D22,D32,D33 [N·mm]",
        f"SSPE,{_fmt(E_ts['E11'])},{_fmt(E_ts['E21'])},{_fmt(E_ts['E22'])}   ! 횡전단 E11,E21,E22 [N/mm] — 1차 추정, 감도 확인 필요",
        f"SSPM,{_fmt(dens_area)}                 ! 단위두께 가정 밀도 = ρ h λ [tonne/mm²]",
    ]
    return "\n".join(lines)


def chevron_gens_deck(g: CorrugationGeometry, m: Sheet, rho: float, beta_deg: float,
                      mode: str = "plate_axes", secid_right: int = 11, secid_left: int = 12) -> str:
    """쉐브론 판 두 영역의 GENS 절 정의 블록.
    mode='plate_axes': 판축으로 회전한 K(16/26 포함) 를 넣고 요소좌표계는 판축(x_p,y_p)에 맞춤 (ESYS 회전 불필요).
    mode='local_axes': 국부축 K(직교이방)를 넣고 요소좌표계를 능선 방향으로 ∓β 회전 (ESYS 로컬 CS 필요)."""
    K = adopted(g, m)
    E_ts, dens = transverse_shear_estimate(g, m), mass_per_area(rho, g)
    hdr = [
        "!=== SPHX equivalent plate sections (generated by calc/equivalent_plate.py) — units mm-N-s-MPa-tonne ===",
        f"! profile={g.profile} p={g.p} H={g.H} t={g.t}  E={m.E} nu={m.nu} rho={rho}  beta={beta_deg} deg  mode={mode}",
        f"! equivalent thicknesses: t_b={equivalent_thicknesses(K, g)['t_b']:.3f} mm (bending-membrane), t_s={equivalent_thicknesses(K, g)['t_s']:.3f} mm (shear-twist), t_m={equivalent_thicknesses(K, g)['t_m']:.3f} mm (mass)",
        "! generalized strain order assumed {eps_x, eps_y, gam_xy, kap_x, kap_y, kap_xy(engineering twist)} — verify with single-element twist test (docs/09 §5)",
    ]
    if mode == "plate_axes":
        KR, KL = zone_stiffness(K, beta_deg, "right"), zone_stiffness(K, beta_deg, "left")
        body = [gens_snippet(KR, secid_right, "SPHX_R", E_ts, dens, f"right zone (+beta): theta = -{beta_deg} deg, includes A16/A26/D16/D26"),
                gens_snippet(KL, secid_left, "SPHX_L", E_ts, dens, f"left zone (-beta): theta = +{beta_deg} deg, A16/A26/D16/D26 sign-reversed"),
                "! element coordinate system for both zones = plate axes (x_p, y_p); apex line must coincide with element edges"]
    else:
        body = [gens_snippet(K, secid_right, "SPHX_LOC", E_ts, dens, "local-axis (ridge) stiffness, orthotropic; same section for both zones"),
                f"! use ESYS with a local CS rotated by -{beta_deg} deg (right zone) and +{beta_deg} deg (left zone) about z_p so that element x = across-ridge direction"]
    return "\n".join(hdr + body)


# ----------------------------------------------------------------------------- 자체 검증
def _selfcheck() -> None:
    g = arc_tangent(9.0, 3.2, 0.6, R_c=1.5); m = Sheet(193000.0, 0.3); K = adopted(g, m)
    th = equivalent_thicknesses(K, g)
    print("등가 두께 (가정 기본값):", {k: round(v, 4) for k, v in th.items()})
    assert th["t_b_spread"] < 1e-9, "채택식에서는 11/12/22 의 D/A 비가 동일해야 함"
    hoA, hoD = homogeneous_orthotropic(K, th["t_b"], "A"), homogeneous_orthotropic(K, th["t_b"], "D")
    print("균질 직교이방 단일층 (t_e = t_b):")
    print(f"  E_x={hoA['E_x']:.1f} MPa  E_y={hoA['E_y']:.1f} MPa  nu_xy={hoA['nu_xy']:.4f}  G(A기준)={hoA['G_from_A']:.1f}  G(D기준)={hoA['G_from_D']:.1f} MPa  → 비 {hoA['G_from_A']/hoA['G_from_D']:.1f}")
    print(f"  잔차 G=A기준: {hoA['resid']}\n  잔차 G=D기준: {hoD['resid']}")
    assert abs(hoA["resid"]["D11"]) < 1e-9 and abs(hoA["resid"]["D22"]) < 1e-9
    print("\n" + chevron_gens_deck(g, m, 7.9e-9, 60.0))
    print("\nselfcheck OK")


if __name__ == "__main__":
    _selfcheck()
