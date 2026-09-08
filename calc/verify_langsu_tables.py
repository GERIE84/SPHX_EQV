"""Lang & Su (2022) Table 5 / Table 6 재계산 스크립트.

docs/02_literature.md §1.5에 전사한 등가 강성식(Xia et al. 2012, Lang & Su 2022)이
논문 Table 5(사인형), Table 6(사다리꼴)의 값을 재현하는지 확인한다.

기호: c 반주기, f 반높이, h 원판 두께, l 반주기 전개길이,
      I1 = ∫(dx/ds)^2 ds, I2 = ∫ z^2 ds (한 주기 2l 구간).
좌표: x = 파형 진행 방향(능선에 수직, 유연), y = 능선 방향(강).
"""
import numpy as np
from scipy.integrate import quad


def xia_2012(E, nu, c, l, I1, I2, h):
    """Xia, Friswell, Saavedra Flores (2012) — Lang & Su Eq. (7)."""
    k = 2 * c * E * h**3 / (I1 * h**2 + 12 * I2)
    return dict(
        A11=k / (1 - nu**2),
        A12=nu * k / (1 - nu**2),
        A22=nu**2 * k / (1 - nu**2) + l * E * h / c,
        A66=c * E * h / (2 * l * (1 + nu)),
        D11=c * E * h**3 / (12 * l * (1 - nu**2)),
        D12=nu * c * E * h**3 / (12 * l * (1 - nu**2)),
        D22=(12 * I2 * E * h + I1 * E * h**3) / (24 * c * (1 - nu**2)),
        D66=l * E * h**3 / (24 * c * (1 + nu)),
    )


def lang_su_2022(E, nu, c, l, I1, I2, h):
    """Lang & Su (2022) — Table 2 'Equivalent model derived in this study'."""
    k = 2 * c * E * h**3 / (I1 * h**2 + 12 * I2 * (1 - nu**2))
    return dict(
        A11=k,
        A12=nu * k,
        A22=nu**2 * k + (1 - nu**2) * l * E * h / c,
        A66=c * E * h / (2 * l * (1 + nu)),
        D11=c * E * h**3 / (12 * l * (1 - nu**2)),
        D12=nu * c * E * h**3 / (12 * l * (1 - nu**2)),
        D22=(12 * I2 * E * h + I1 * E * h**3) / (24 * c),
        D66=l * E * h**3 / (24 * c * (1 + nu)),
    )


def sinusoidal_geometry(c, f):
    """Eq. (14): z = f sin(pi x / c)."""
    zp = lambda x: f * np.pi / c * np.cos(np.pi * x / c)
    ds = lambda x: np.sqrt(1 + zp(x) ** 2)
    l = quad(ds, 0, c)[0]
    I1 = quad(lambda x: 1 / ds(x), 0, 2 * c)[0]  # cos^2(theta) ds = dx/ds dx
    I2 = quad(lambda x: (f * np.sin(np.pi * x / c)) ** 2 * ds(x), 0, 2 * c)[0]
    return l, I1, I2


def trapezoidal_geometry(c, f, alpha):
    """Eq. (15): 경사각 alpha [rad]."""
    l = 2 * f / np.sin(alpha) + c - 2 * f / np.tan(alpha)
    I1 = 4 * f * np.cos(alpha) ** 2 / np.sin(alpha) + 2 * c - 4 * f / np.tan(alpha)
    I2 = 4 * f**3 / (3 * np.sin(alpha)) + 2 * f**2 * (c - 2 * f / np.tan(alpha))
    return l, I1, I2


def report(title, calc_xia, calc_ls, ref_xia, ref_ls):
    print(f"\n=== {title} ===")
    print(f"{'':6} {'Xia calc':>12} {'Xia paper':>12} {'LS calc':>12} {'LS paper':>12}")
    for k in calc_xia:
        print(f"{k:6} {calc_xia[k]:12.6g} {ref_xia[k]:12.6g} {calc_ls[k]:12.6g} {ref_ls[k]:12.6g}")


if __name__ == "__main__":
    # Table 4 — 사인형
    E, nu, c, f, h = 30e9, 0.2, 0.32, 0.11, 0.005
    l, I1, I2 = sinusoidal_geometry(c, f)
    print(f"sinusoidal: l={l:.5f} I1={I1:.5f} I2={I2:.6f}")
    report(
        "Table 5 sinusoidal",
        xia_2012(E, nu, c, l, I1, I2, h),
        lang_su_2022(E, nu, c, l, I1, I2, h),
        dict(A11=47613, A12=9523, A22=187.08e6, A66=50.113e6, D11=261.004, D12=52.20, D22=1068260, D66=162.39),
        dict(A11=47612, A12=9522, A22=179.60e6, A66=50.113e6, D11=261.004, D12=52.20, D22=1025530, D66=162.39),
    )

    # Table 4 — 사다리꼴
    E, nu, c, f, h, a = 21e9, 0.3, 0.0508, 0.0127, 0.00635, np.radians(45)
    l, I1, I2 = trapezoidal_geometry(c, f, a)
    print(f"\ntrapezoidal: l={l:.5f} I1={I1:.5f} I2={I2:.6e}")
    report(
        "Table 6 trapezoidal",
        xia_2012(E, nu, c, l, I1, I2, h),
        lang_su_2022(E, nu, c, l, I1, I2, h),
        dict(A11=4.052e6, A12=1.216e6, A22=161.332e6, A66=42.489e6, D11=407.917, D12=122.375, D22=17809, D66=208.032),
        dict(A11=4.042e6, A12=1.213e6, A22=146.844e6, A66=42.489e6, D11=407.917, D12=122.375, D22=16206, D66=208.032),
    )
