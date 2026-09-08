"""Ye, Berdichevsky & Yu (2014) 및 고전 등가 강성식 재계산 스크립트.

docs/02_literature.md §3에 전사한 수식이 Ye (2014) Table 1(사인형) 및
Lang & Su (2022) Table 5의 Briassoulis 열을 재현하는지 확인한다.

Ye 기호: ε 주기(=2c), T 주름 높이(중앙면 기준, =f), S 한 주기 전개길이(=2l),
         X = x/ε ∈ [-1/2, 1/2], x3 = ε φ(X), φ = dφ/dX, a = 1+φ², ⟨f⟩ = ∫_{-1/2}^{1/2} f dX
대응:   ⟨√a⟩ = S/ε = l/c,  ε²⟨φ²√a⟩ = I2/(2c),  ⟨1/√a⟩ = I1/(2c)
"""
import numpy as np
from scipy.integrate import quad


def cell_averages(phi, dphi, d2phi, eps, h):
    """한 주기 평균량. phi, dphi(=φ), d2phi(=φ') 는 X의 함수."""
    a = lambda X: 1 + dphi(X) ** 2
    sa = lambda X: np.sqrt(a(X))
    avg = lambda g: quad(g, -0.5, 0.5, limit=200)[0]
    out = {}
    out["sqrt_a"] = avg(sa)                                    # ⟨√a⟩
    out["inv_sqrt_a"] = avg(lambda X: 1 / sa(X))               # ⟨1/√a⟩
    out["phi2_sqrt_a"] = avg(lambda X: phi(X) ** 2 * sa(X))    # ⟨φ²√a⟩
    out["sqrt_a_phi"] = avg(lambda X: sa(X) * phi(X))          # ⟨√a φ⟩  (대칭이면 0)
    B = out["sqrt_a_phi"] / out["sqrt_a"]                      # Eq. (12) ℬ
    # 𝒜(X) = -∫_0^X √a φ dY + ℬ ∫_0^X √a dY   — Eq. (13)
    def calA(X):
        s = 1.0 if X >= 0 else -1.0
        i1 = quad(lambda Y: sa(Y) * phi(Y), 0, X)[0]
        i2 = quad(sa, 0, X)[0]
        return -i1 + B * i2
    out["phiA"] = avg(lambda X: dphi(X) * calA(X))             # ⟨φ𝒜⟩
    out["B"] = B
    # α1, D66 관련 (Eq. 12, 19)
    corr = lambda X: 1 + d2phi(X) ** 2 * h**2 / (48 * eps**2 * a(X) ** 3)
    out["alpha1"] = 1 / avg(lambda X: sa(X) / corr(X))
    out["D66_int"] = avg(lambda X: sa(X) * h**2 / 3
                         - (1 / sa(X)) * (h**4 * d2phi(X) ** 2 / (144 * eps**2 * a(X) ** 2)) / corr(X))
    return out


def ye_2014_symmetric(E, nu, eps, h, av):
    """Ye et al. (2014) Eq. (19) — 대칭 주름 (B=0)."""
    mu = E / (2 * (1 + nu))
    C = -12 * av["phiA"] * eps**2 / h**2 - av["inv_sqrt_a"]     # 𝒞, Eq. (12)
    A11 = (E / (1 - nu**2)) * 12 * eps**2 * av["phiA"] / (h * C**2) \
        + (E * h / (1 - nu**2)) * av["inv_sqrt_a"] / C**2
    D11 = E * h**3 / (12 * (1 - nu**2) * av["sqrt_a"])
    return dict(
        A11=A11,
        A12=nu * A11,
        A22=E * h * av["sqrt_a"] + nu**2 * A11,
        A66=mu * h * av["alpha1"],
        D11=D11,
        D12=nu * D11,
        D22=E * h * eps**2 * av["phi2_sqrt_a"] + E * h**3 / 12 * av["inv_sqrt_a"] + nu**2 * D11,
        D66=mu * h / 4 * av["D66_int"],
    )


def ye_2014_thin_symmetric(E, nu, eps, h, av):
    """Ye et al. (2014) Eq. (29) — h/ε ≪ 1 선두항, 대칭 주름."""
    mu = E / (2 * (1 + nu))
    D11 = E * h**3 / (12 * (1 - nu**2) * av["sqrt_a"])
    A11 = E * h**3 / (12 * (1 - nu**2) * eps**2 * av["phi2_sqrt_a"])
    return dict(A11=A11, A12=nu * A11, A22=E * h * av["sqrt_a"], A66=mu * h / av["sqrt_a"],
                D11=D11, D12=nu * D11, D22=E * h * eps**2 * av["phi2_sqrt_a"], D66=mu * h**3 / 12 * av["sqrt_a"])


def seydel_1931_bending(E, nu, eps, S, Iy, h):
    """Ye Eq. (5). Iy = h ε² ⟨φ²√a⟩ (Ye Eq. 4)."""
    return dict(D11=eps / S * E * h**3 / (12 * (1 - nu**2)), D12=0.0, D22=E * Iy,
                D66=S / eps * E * h**3 / (24 * (1 + nu)))


def classical_extension_1960s(E, nu, eps, S, T, h):
    """Ye Eq. (8). T = 주름 높이; Ye는 T² := ⟨x3²√a⟩/2 로 두면 정확하다고 지적."""
    A11 = E * h**3 / (6 * (1 - nu**2) * T**2)
    return dict(A11=A11, A12=nu * A11, A22=S / eps * E * h, A66=eps / S * E * h / (2 * (1 + nu)))


def briassoulis_1986(E, nu, eps, S, T, h):
    """Ye Eq. (6), (9) — 사인형 x3 = T sin(2πx/ε) 가정."""
    D11 = eps / S * E * h**3 / (12 * (1 - nu**2))
    A11 = E * h**3 / (h**2 + 6 * (1 - nu**2) * T**2 * (S**2 / eps**2 - S / (2 * np.pi * eps) * np.sin(2 * np.pi * S / eps)))
    return dict(A11=A11, A12=nu * A11, A22=S / eps * E * h, A66=E * h / (2 * (1 + nu)),
                D11=D11, D12=nu * D11, D22=E * h * T**2 / 2 + E * h**3 / (12 * (1 - nu**2)),
                D66=E * h**3 / (24 * (1 + nu)))


def show(title, calc, ref):
    print(f"\n=== {title} ===")
    for k in calc:
        r = ref.get(k)
        rs = f"{r:12.6g}" if r is not None else f"{'-':>12}"
        print(f"  {k:4} calc={calc[k]:12.6g}  paper={rs}")


if __name__ == "__main__":
    # Ye (2014) §7.1 사인형: ε=0.64, T=0.11, h=0.005, E=30 GPa, ν=0.2  (= Lang & Su Table 4 사인형)
    E, nu, eps, T, h = 30e9, 0.2, 0.64, 0.11, 0.005
    phi = lambda X: T / eps * np.sin(2 * np.pi * X)
    dphi = lambda X: T / eps * 2 * np.pi * np.cos(2 * np.pi * X)
    d2phi = lambda X: -T / eps * (2 * np.pi) ** 2 * np.sin(2 * np.pi * X)
    av = cell_averages(phi, dphi, d2phi, eps, h)
    S = eps * av["sqrt_a"]
    Iy = h * eps**2 * av["phi2_sqrt_a"]
    print(f"S={S:.5f} (=2l)  Iy={Iy:.6e}  ⟨√a⟩={av['sqrt_a']:.5f}  ⟨1/√a⟩={av['inv_sqrt_a']:.5f}  ⟨φ𝒜⟩={av['phiA']:.6e}  ℬ={av['B']:.2e}")

    show("Ye (2014) Eq.(19) full symmetric  [Ye Table 1 'Present']",
         ye_2014_symmetric(E, nu, eps, h, av),
         dict(A11=47613, A12=9523, A22=1.8708e8, A66=5.0113e7, D11=261.004, D12=52.20, D22=1025540, D66=162.39))
    show("Ye (2014) Eq.(29) thin-plate leading order  [= Lang & Su Table 2 'Ye symmetric' column]",
         ye_2014_thin_symmetric(E, nu, eps, h, av), {})
    ref_classical = dict(A11=53805, A12=10761, A22=1.8708e8, A66=5.0113e7, D11=261.004, D12=0, D22=1025270, D66=162.39)
    show("Seydel (1931) Eq.(5) + 1960s extension Eq.(8)  [Ye Table 1 'Eqs.(5),(8)']",
         {**classical_extension_1960s(E, nu, eps, S, T, h), **seydel_1931_bending(E, nu, eps, S, Iy, h)}, ref_classical)
    show("Briassoulis (1986) Eq.(6),(9)  [Lang & Su Table 5 'Briassoulis']",
         briassoulis_1986(E, nu, eps, S, T, h),
         dict(A11=39639, A12=7928, A22=187.08e6, A66=62.5e6, D11=261.004, D12=52.20, D22=907830, D66=260.42))
