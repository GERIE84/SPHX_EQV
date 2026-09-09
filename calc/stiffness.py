"""일방향 주름판(단위셀, 국부축 x⊥능선 / y∥능선)의 등가 직교이방성 강성 — B-2 구현.

docs/06_equivalent_stiffness.md 의 수식을 구현한다.  입력은 geometry.CorrugationGeometry (c, l, I1, I2, h)
와 원판 물성 E, nu.  출력은 등가 판의 A(면내), D(굽힘) 행렬 성분 (B 연성은 대칭 단면에서 0).

기호:  lam = l/c (전개비),  J1 = I1/(2c) = <1/√a>,  J2 = I2/(2c) = ε²<φ²√a>  [length²],  h = t,  G = E/(2(1+ν))

모델:
    xia      : Xia, Friswell & Saavedra Flores (2012)          — RVE, Dirichlet BC (D22 상한)
    lang_su  : Lang & Su (2022)                                 — RVE + Kirchhoff 단순화
    ye_thin  : Ye, Berdichevsky & Yu (2014) Eq.(29)             — VAM, 얇은판·대칭 선두항
    ye_full  : Ye (2014) Eq.(19)                                — VAM, 대칭 단면 전체식 (A66, D66 곡률 보정 포함)
    adopted  : 채택식 = ye_full 에서 A66, D66 의 (h κ)² 보정만 생략 (docs/06 §3)
구성관계 관례:  {N_x,N_y,N_xy} = A {ε_x, ε_y, γ_xy},  {M_x,M_y,M_xy} = D {κ_x, κ_y, 2κ_xy}  (공학 전단변형률/비틀림)
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, Iterable, Optional, Tuple

import numpy as np

from geometry import CorrugationGeometry, arc_tangent, sinusoidal, trapezoidal, round_corrugation

KEYS = ("A11", "A12", "A22", "A66", "D11", "D12", "D22", "D66")


@dataclass
class Sheet:
    """원판(등방성) 물성."""
    E: float
    nu: float

    @property
    def G(self) -> float:
        return self.E / (2 * (1 + self.nu))


def _base(g: CorrugationGeometry, m: Sheet):
    return g.h, g.l / g.c, g.I1 / (2 * g.c), g.I2 / (2 * g.c), m.E, m.nu, m.G


# ----------------------------------------------------------------------------- 모델별 등가 강성
def xia(g: CorrugationGeometry, m: Sheet) -> Dict[str, float]:
    h, lam, J1, J2, E, nu, G = _base(g, m)
    A11 = E * h**3 / ((1 - nu**2) * (12 * J2 + h**2 * J1))
    D11 = E * h**3 / (12 * (1 - nu**2) * lam)
    return dict(A11=A11, A12=nu * A11, A22=nu**2 * A11 + E * h * lam, A66=G * h / lam,
                D11=D11, D12=nu * D11, D22=(E * h * J2 + E * h**3 * J1 / 12) / (1 - nu**2), D66=G * h**3 * lam / 12)


def lang_su(g: CorrugationGeometry, m: Sheet) -> Dict[str, float]:
    h, lam, J1, J2, E, nu, G = _base(g, m)
    A11 = E * h**3 / (12 * J2 * (1 - nu**2) + h**2 * J1)
    D11 = E * h**3 / (12 * (1 - nu**2) * lam)
    return dict(A11=A11, A12=nu * A11, A22=nu**2 * A11 + (1 - nu**2) * E * h * lam, A66=G * h / lam,
                D11=D11, D12=nu * D11, D22=E * h * J2 + E * h**3 * J1 / 12, D66=G * h**3 * lam / 12)


def ye_thin(g: CorrugationGeometry, m: Sheet) -> Dict[str, float]:
    h, lam, J1, J2, E, nu, G = _base(g, m)
    A11 = E * h**3 / (12 * (1 - nu**2) * J2)
    D11 = E * h**3 / (12 * (1 - nu**2) * lam)
    return dict(A11=A11, A12=nu * A11, A22=E * h * lam, A66=G * h / lam,
                D11=D11, D12=nu * D11, D22=E * h * J2, D66=G * h**3 * lam / 12)


def curvature_segments(g: CorrugationGeometry, n: int = 2000) -> Iterable[Tuple[float, float]]:
    """한 주기 중앙면을 (호 길이 ds, 곡률 κ) 조각으로 반환.  Ye 전체식의 A66, D66 곡률 보정용."""
    if g.profile in ("arc_tangent", "round"):
        if g.profile == "round":
            Rc = Rv = g.extra["R"]; alpha = math.pi / 2; TL = 2 * g.extra["L"]
        else:
            Rc, Rv, alpha, TL = g.extra["R_c"], g.extra["R_v"], math.radians(g.extra["alpha_deg"]), g.extra["T_L"]
        # 한 주기 = 2 × (산 원호 + 접선 + 골 원호)
        return [(2 * Rc * alpha, 1 / Rc), (2 * TL, 0.0), (2 * Rv * alpha, 1 / Rv)]
    if g.profile == "trapezoidal":
        return [(2 * g.l, 0.0)]                      # 날카로운 모서리의 곡률은 무시 (문서 §3.4 참조)
    if g.profile == "sinusoidal":
        f, c = g.f, g.c
        x = np.linspace(0, 2 * c, n + 1); xm = 0.5 * (x[1:] + x[:-1])
        zp = -f * math.pi / c * np.sin(math.pi * xm / c)
        zpp = -f * (math.pi / c) ** 2 * np.cos(math.pi * xm / c)
        sa = np.sqrt(1 + zp**2)
        return list(zip(sa * np.diff(x), zpp / sa**3))
    raise ValueError(g.profile)


def ye_full(g: CorrugationGeometry, m: Sheet) -> Dict[str, float]:
    """Ye (2014) Eq.(19), 대칭 단면.  <φ𝒜> = <φ²√a> 를 사용 (부분적분, 대칭)."""
    h, lam, J1, J2, E, nu, G = _base(g, m)
    eps = 2 * g.c
    C = -(12 * J2 / h**2 + J1)                                   # 𝒞
    A11 = (E / (1 - nu**2)) * (12 * J2 / h + h * J1) / C**2      # = E h³ / ((1-ν²)(12 J2 + h² J1))  (Xia 와 동일)
    D11 = E * h**3 / (12 * (1 - nu**2) * lam)
    # A66 = μ h α1 ,  α1 = 1 / < √a / (1 + h²κ²/48) >  ;  D66 = μh/4 < √a h²/3 − (1/√a)(h⁴κ²/144)/(1+h²κ²/48) >
    inv_a1 = 0.0; d66_int = 0.0
    for ds, kap in curvature_segments(g):
        corr = 1 + h**2 * kap**2 / 48
        inv_a1 += ds / corr
        d66_int += ds * (h**2 / 3 - (h**4 * kap**2 / 144) / corr)
    alpha1 = eps / inv_a1
    A66 = G * h * alpha1
    D66 = G * h / 4 * d66_int / eps
    return dict(A11=A11, A12=nu * A11, A22=E * h * lam + nu**2 * A11, A66=A66,
                D11=D11, D12=nu * D11, D22=E * h * J2 + E * h**3 * J1 / 12 + nu**2 * D11, D66=D66)


def adopted(g: CorrugationGeometry, m: Sheet) -> Dict[str, float]:
    """채택식 (docs/06 §3): Ye 전체식에서 A66, D66 의 곡률 보정만 생략."""
    h, lam, J1, J2, E, nu, G = _base(g, m)
    A11 = E * h**3 / ((1 - nu**2) * (12 * J2 + h**2 * J1))
    D11 = E * h**3 / (12 * (1 - nu**2) * lam)
    return dict(A11=A11, A12=nu * A11, A22=E * h * lam + nu**2 * A11, A66=G * h / lam,
                D11=D11, D12=nu * D11, D22=E * h * J2 + E * h**3 * J1 / 12 + nu**2 * D11, D66=G * h**3 * lam / 12)


MODELS = dict(xia=xia, lang_su=lang_su, ye_thin=ye_thin, ye_full=ye_full, adopted=adopted)


# ----------------------------------------------------------------------------- 등가 공학상수 변환 (B-4 예비)
def engineering_constants(K: Dict[str, float], t_e: float, basis: str = "membrane") -> Dict[str, float]:
    """등가 두께 t_e 를 가정했을 때의 직교이방성 공학상수.
    basis='membrane': A 행렬로부터,  basis='bending': D 행렬로부터 (Lang & Su Eq. 16 형태).
    두 결과가 다르면 단일 두께의 균질 직교이방성 판으로는 A, D 를 동시에 만족시킬 수 없다는 뜻 (B-4 에서 다룸)."""
    if basis == "membrane":
        a11, a12, a22, a66, s = K["A11"], K["A12"], K["A22"], K["A66"], t_e
    else:
        a11, a12, a22, a66, s = K["D11"], K["D12"], K["D22"], K["D66"], t_e**3 / 12
    det = a11 * a22 - a12**2
    return dict(E_x=det / (a22 * s), E_y=det / (a11 * s), nu_xy=a12 / a22, nu_yx=a12 / a11, G_xy=a66 / s)


def rho_equivalent(rho: float, g: CorrugationGeometry) -> float:
    """등가 밀도 (질량 보존, 두께 h 유지 시): ρ_eq = ρ l / c  (Lang & Su Eq. 17)."""
    return rho * g.l / g.c


# ----------------------------------------------------------------------------- 표 출력
def table(g: CorrugationGeometry, m: Sheet, models=("xia", "lang_su", "ye_thin", "ye_full", "adopted"), ref: Optional[str] = "adopted") -> str:
    res = {k: MODELS[k](g, m) for k in models}
    w = 12
    head = f"{'':6}" + "".join(f"{k:>{w}}" for k in models) + (f"   {'ratio to ' + ref:>18}" if ref else "")
    lines = [head]
    for key in KEYS:
        row = f"{key:6}" + "".join(f"{res[k][key]:{w}.5g}" for k in models)
        if ref:
            row += "   " + " ".join(f"{res[k][key] / res[ref][key]:6.3f}" for k in models if k != ref)
        lines.append(row)
    return "\n".join(lines)


# ----------------------------------------------------------------------------- 자체 검증
def _selfcheck() -> None:
    ok = True

    def chk(name, got, ref, tol=3e-4):
        nonlocal ok
        rel = abs(got - ref) / abs(ref)
        ok &= rel < tol
        print(f"  [{'OK ' if rel < tol else 'BAD'}] {name}: {got:.6g} (ref {ref:.6g}, rel {rel:.1e})")

    print("1) 사인형 (Lang & Su Table 4/5, Ye Table 1): E=30 GPa, ν=0.2, c=0.32, f=0.11, h=0.005 m")
    g = sinusoidal(p=0.64, H=0.22, t=0.005); m = Sheet(30e9, 0.2)
    refs = dict(
        xia=dict(A11=47613, A12=9523, A22=187.08e6, A66=50.113e6, D11=261.004, D12=52.20, D22=1068260, D66=162.39),
        lang_su=dict(A11=47612, A12=9522, A22=179.60e6, A66=50.113e6, D11=261.004, D12=52.20, D22=1025530, D66=162.39),
        ye_full=dict(A11=47613, A12=9523, A22=1.8708e8, A66=5.0113e7, D11=261.004, D12=52.20, D22=1025540, D66=162.39),
    )
    for name, r in refs.items():
        K = MODELS[name](g, m)
        for k in KEYS:
            chk(f"{name}.{k}", K[k], r[k], 5e-4)

    print("2) 사다리꼴 (Lang & Su Table 6 / Xia Table 5): E=21 GPa, ν=0.3, c=0.0508, f=0.0127, h=0.00635, α=45°")
    g = trapezoidal(p=0.1016, H=0.0254, t=0.00635, alpha_deg=45); m = Sheet(21e9, 0.3)
    refs = dict(
        xia=dict(A11=4.052e6, A12=1.216e6, A22=161.332e6, A66=42.489e6, D11=407.917, D12=122.375, D22=17809, D66=208.032),
        lang_su=dict(A11=4.042e6, A12=1.213e6, A22=146.844e6, A66=42.489e6, D11=407.917, D12=122.375, D22=16206, D66=208.032),
    )
    for name, r in refs.items():
        K = MODELS[name](g, m)
        for k in KEYS:
            chk(f"{name}.{k}", K[k], r[k], 5e-4)
    # MSG-TW 기준값(Table 6): D22=16242, A22=161.332e6 → 채택식과 비교
    K = adopted(g, m)
    chk("adopted.D22 vs MSG-TW 16242", K["D22"], 16242, 5e-3)
    chk("adopted.A22 vs MSG-TW 161.332e6", K["A22"], 161.332e6, 5e-3)

    print("3) 평판 극한 (H→0): 모든 모델이 등방 판 강성으로 환원되는가 (사인형 H=1e-6)")
    g = sinusoidal(p=10.0, H=1e-6, t=1.0); m = Sheet(200e3, 0.3)
    Eh, D = m.E * 1.0 / (1 - m.nu**2), m.E / (12 * (1 - m.nu**2))
    for name in ("xia", "ye_full", "adopted"):
        K = MODELS[name](g, m)
        chk(f"{name}.A11→Eh/(1-ν²)", K["A11"], Eh, 1e-3); chk(f"{name}.A22→Eh/(1-ν²)", K["A22"], Eh, 1e-3)
        chk(f"{name}.D11→D", K["D11"], D, 1e-3); chk(f"{name}.D22→D", K["D22"], D, 1e-3); chk(f"{name}.D66→Gh³/12", K["D66"], m.G / 12, 1e-3)
    K = lang_su(g, m)
    print(f"  (참고) lang_su 평판 극한: A11/(Eh/(1-ν²)) = {K['A11']/Eh:.4f}, D22/D = {K['D22']/D:.4f}  → 평판으로 환원되지 않음 (Ye 의 지적과 일치)")
    print("ALL OK" if ok else "SOME CHECKS FAILED")


if __name__ == "__main__":
    _selfcheck()
    print("\n=== 가정 기본값 (docs/05 §5): arc_tangent p=9, H=3.2, t=0.6, R=1.5 mm, 316L E=193 GPa, ν=0.3  [N, mm] ===")
    g = arc_tangent(9.0, 3.2, 0.6, R_c=1.5); m = Sheet(193000.0, 0.3)
    print(g.summary()); print()
    print(table(g, m))
    K = adopted(g, m)
    print("\n채택식 → 공학상수 (t_e = t 로 가정한 경우):")
    for basis in ("membrane", "bending"):
        print(f"  {basis:9}", {k: round(v, 4 if k.startswith('nu') else 1) for k, v in engineering_constants(K, g.t, basis).items()})
    print(f"\n직교이방성 비:  A22/A11 = {K['A22']/K['A11']:.1f},  D22/D11 = {K['D22']/K['D11']:.1f},  A66/A11 = {K['A66']/K['A11']:.2f}")
