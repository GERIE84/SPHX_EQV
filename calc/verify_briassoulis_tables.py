"""Briassoulis (1986) Table 1–3 재계산 및 독립 FE 벤치마크 비교.

D. Briassoulis, Equivalent orthotropic properties of corrugated sheets, Comput. Struct. 23(2) (1986) 129–138.
§4.1 형상: 원호+접선(arc-and-tangent) 표준 주름  c=2 in, f=0.21875 in, l=2.046 in, α=17.54°, R=2.0208 in, t=0.25 in
재료: E=30e6 psi, μ=0.3.  FE: 9절점 Lagrangian 쉘, 한 주기(2c) × L=3 in 에 일정 변형률/곡률 상태 부과.

관례 주의: 그의 비틀림 강성 B_xy 는 M_xy = B_xy · w_,xy 로 정의 → Ye/Xia/본 프로젝트의 D66 (M_xy = D66 · 2κ_xy) 의 2배.
"""
import math
import sys

from geometry import arc_tangent
from stiffness import KEYS, MODELS, Sheet

E, mu, t, f, c = 30e6, 0.3, 0.25, 0.21875, 2.0
g = arc_tangent(p=2 * c, H=2 * f, t=t, R_c=2.0208, R_v=2.0208)
m = Sheet(E, mu)
l = g.l

# 논문 수치 (kips/in → lb/in ; kips·in/in → lb·in/in).  D66 열은 그의 B_xy 값 그대로.
FE = dict(A11=1440e3, A12=439e3, A22=7709e3, A66=2800e3, D11=41.5e3, D12=12.5e3, D22=216.5e3, Bxy=30.0e3)
OLD = dict(A11=1794e3, A12=550e3, A22=7673e3, A66=2820e3, D11=42.0e3, D12=0.0, D22=187.6e3, Bxy=30.7e3)
PRESENT = dict(A11=1421e3, A12=434e3, A22=7673e3, A66=2885e3, D11=42.0e3, D12=12.4e3, D22=222.4e3, Bxy=30.0e3)


def briassoulis_present(E, mu, t, f, c, l):
    """Table 3 'Present' (Eq. 8, 9, 10, 11, 13, 14)."""
    A11 = E * t / (1 + (f / t) ** 2 * 6 * (1 - mu**2) * (l**2 / c**2 - (l / (2 * math.pi * c)) * math.sin(2 * math.pi * l / c)))
    D11 = E * t**3 / (12 * (1 - mu**2)) * c / l
    return dict(A11=A11, A12=mu * A11, A22=E * t * l / c, A66=E * t / (2 * (1 + mu)),
                D11=D11, D12=mu * D11, D22=E * t**3 / (12 * (1 - mu**2)) + E * t * f**2 / 2,
                Bxy=E * t**3 / (12 * (1 + mu)))


def stress_concentration_x(f, t):
    """Appendix B Eq. (B2): 일정 ε_x 상태에서 능선(z=f) 내측 섬유 최대응력 / 등가 판 공칭응력 N_x/t."""
    return 1 + 6 * f / t


if __name__ == "__main__":
    print(g.summary())
    print(f"  paper: l=2.046 in, alpha=17.54 deg  (계산 l={l:.4f}, alpha={g.extra['alpha_deg']:.2f} → 공표값 자체 반올림 0.6%)")
    print(f"  shallow: f/c={f/c:.3f}, l/c={l/c:.4f}, t/R={t/2.0208:.3f}, f/t={f/t:.3f}")

    print("\n[1] Briassoulis 'Present' 식 재계산 vs Table 3 수치 (kips/in)")
    b = briassoulis_present(E, mu, t, f, c, l)
    for k in list(KEYS[:7]) + ["Bxy"]:
        print(f"  {k:4} calc={b[k]/1e3:9.1f}  paper={PRESENT[k]/1e3:9.1f}")

    print("\n[2] 본 프로젝트 모델 vs Briassoulis FE (Table 1, 2)  — kips/in, 괄호 = 모델/FE")
    names = list(MODELS)
    print(f"  {'':5}{'FE':>9}" + "".join(f"{n:>11}" for n in names))
    for key in KEYS:
        row = {n: MODELS[n](g, m)[key] for n in names}
        fe = FE["Bxy"] if key == "D66" else FE[key]
        if key == "D66":
            row = {n: 2 * v for n, v in row.items()}   # → B_xy 관례로 변환
        print(f"  {key + ('*' if key == 'D66' else ''):5}{fe/1e3:9.1f}" + "".join(f"{row[n]/1e3:7.1f}({row[n]/fe:5.3f})" for n in names))
    print("  * D66 행은 2·D66 = B_xy 로 환산해 비교")

    print("\n[3] 비틀림 관례 확인")
    D66_ye_conv = E * t**3 / (24 * (1 + mu))
    print(f"  Briassoulis B_xy = Et³/(12(1+μ)) = {b['Bxy']/1e3:.2f} kips·in/in  = 2 × {D66_ye_conv/1e3:.2f} (= Ye 관례 D66 = Et³/(24(1+ν)))")
    # Lang & Su Table 5 사인형 판 (E=30 GPa, ν=0.2, h=0.005 m): 260.42 vs Ye 130.21
    Eh3 = 30e9 * 0.005**3
    print(f"  사인형 검증판: Et³/(12(1+ν)) = {Eh3/(12*1.2):.2f} (Lang & Su Table 5 'Briassoulis' 260.42)  vs  Et³/(24(1+ν)) = {Eh3/(24*1.2):.2f} (Ye 전사 130.21)")
    print("  → Lang & Su 는 B_xy 관례 값을 D66 관례로 변환하지 않고 옮겼음. Ye 전사가 옳다.")

    print("\n[4] 응력집중 (Appendix B, Eq. B2): K_t = 1 + 6 f/t")
    print(f"  Briassoulis 판 f/t={f/t:.3f} → K_t={stress_concentration_x(f, t):.2f};   PHE 기본값 f=1.6, t=0.6 → K_t={stress_concentration_x(1.6, 0.6):.1f}")
