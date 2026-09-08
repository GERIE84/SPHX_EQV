"""Xia, Friswell & Saavedra Flores (2012) Table 2–5 재계산 스크립트.

docs/02_literature.md §3에 전사한 수식을 원문 Table 5(사다리꼴, E=21 GPa, ν=0.3,
c=0.0508, f=0.0127, t=0.00635, α=45°)의 수치와 대조한다.

Xia 기호: c 반주기, l 반주기 전개길이, I1=∫(dx/ds)²ds, I2=∫z²ds (한 주기 2l),
원판(sheet) 강성 A11..D66 은 일반 직교이방성 입력 가능 (여기서는 등방성).
"""
import numpy as np
from scipy.integrate import quad


def isotropic_sheet(E, nu, t):
    """등방성 박판의 국부 강성 (Xia Eq. 7 계수)."""
    k, d = E * t / (1 - nu**2), E * t**3 / (12 * (1 - nu**2))
    return dict(A11=k, A12=nu * k, A22=k, A66=E * t / (2 * (1 + nu)),
                D11=d, D12=nu * d, D22=d, D66=E * t**3 / (24 * (1 + nu)))


def xia_2012(sheet, c, l, I1, I2):
    """Xia Table 2 — 임의 단면, 임의(직교이방성) 원판."""
    s = sheet
    A11 = 2 * c / (I1 / s["A11"] + I2 / s["D11"])          # Eq. (21)
    A12 = s["A12"] / s["A11"] * A11                          # Eq. (23)
    D11 = c / l * s["D11"]                                   # Eq. (36)
    return dict(
        A11=A11,
        A12=A12,
        A22=s["A12"] / s["A11"] * A12 + l / c * (s["A11"] * s["A22"] - s["A12"] ** 2) / s["A11"],  # Eq. (28)
        A66=c / l * s["A66"],
        D11=D11,
        D12=s["D12"] / s["D11"] * D11,
        D22=(I2 * s["A22"] + I1 * s["D22"]) / (2 * c),
        D66=l / c * s["D66"],
    )


def samanta_mukhopadhyay_as_printed(E, nu, t, c, l, I2):
    """Xia Table 3 'Samanta and Mukhopadhyay (1999)' 열, 인쇄된 수식 그대로."""
    return dict(
        A11=2 * c / I2 * E * t**3 / 12,
        A22=l / c * E * t,
        A66=c / l * E * t / (2 * (1 + nu)),
        D11=c / l * E * t**3 / 12,
        D12=0.0,
        D22=E * t / (2 * c) * I2,
        D66=l / c * E * t**3 / (6 * (1 + nu)),
    ) | {"A12": nu * (2 * c / I2 * E * t**3 / 12)}


def trapezoid(c, f, alpha):
    l = 2 * f / np.sin(alpha) + c - 2 * f / np.tan(alpha)
    I1 = 4 * f * np.cos(alpha) ** 2 / np.sin(alpha) + 2 * c - 4 * f / np.tan(alpha)
    I2 = 4 * f**3 / (3 * np.sin(alpha)) + 2 * f**2 * (c - 2 * f / np.tan(alpha))
    return l, I1, I2


def round_corrugation_closed(R, L):
    """Xia Table 4: c=2R, l=πR+2L, I1=πR, I2=4L³/3+2πL²R+8LR²+πR³."""
    return 2 * R, np.pi * R + 2 * L, np.pi * R, 4 * L**3 / 3 + 2 * np.pi * L**2 * R + 8 * L * R**2 + np.pi * R**3


def round_corrugation_numeric(R, L):
    """반주기 = 수직직선 L + 반원 R + 수직직선 L (폭 2R); 전주기는 상하 대칭."""
    # 반주기 상부: I1(cos²θ ds), I2(z² ds)
    I1_half = 0.0 + quad(lambda p: np.sin(p) ** 2 * R, 0, np.pi)[0]           # 직선부 cosθ=0
    I2_half = 2 * quad(lambda z: z**2, 0, L)[0] + quad(lambda p: (L + R * np.sin(p)) ** 2 * R, 0, np.pi)[0]
    return 2 * I1_half, 2 * I2_half


def show(title, calc, ref):
    print(f"\n=== {title} ===")
    for k in ["A11", "A12", "A22", "A66", "D11", "D12", "D22", "D66"]:
        r = ref.get(k)
        ratio = f"{calc[k]/r:8.4f}" if r else "       -"
        print(f"  {k:4} calc={calc[k]:12.6g}  paper={r if r is not None else '-':>12}  ratio={ratio}")


if __name__ == "__main__":
    E, nu, t, c, f, a = 21e9, 0.3, 0.00635, 0.0508, 0.0127, np.radians(45)
    l, I1, I2 = trapezoid(c, f, a)
    print(f"trapezoid: l={l:.5f} I1={I1:.5f} I2={I2:.6e}")

    ref_xia = dict(A11=4.052e6, A12=1.216e6, A22=161.332e6, A66=42.489e6, D11=407.917, D12=122.375, D22=17809, D66=208.032)
    ref_fem = dict(A11=4.051e6, A12=1.215e6, A22=163.910e6, A66=42.797e6, D11=406.512, D12=121.954, D22=17805, D66=207.96)
    ref_sm = dict(A11=4.150e6, A12=1.245e6, A22=176.888e6, A66=42.489e6, D11=371.205, D12=0, D22=31647, D66=208.032)

    show("Xia Table 2 (general form, isotropic sheet)  vs Table 5 'Proposed'", xia_2012(isotropic_sheet(E, nu, t), c, l, I1, I2), ref_xia)
    show("Xia Table 2  vs Table 5 'FEM' (ANSYS SHELL63 unit cell, generalized-strain BC)", xia_2012(isotropic_sheet(E, nu, t), c, l, I1, I2), ref_fem)
    show("Samanta & Mukhopadhyay — Xia Table 3 formulas AS PRINTED  vs Table 5 'S&M' values",
         samanta_mukhopadhyay_as_printed(E, nu, t, c, l, I2), ref_sm)

    # Ye (2014) Eq.(10)이 전사한 Xia A22 vs Xia 원문 Table 2 A22
    S, eps = 2 * l, 2 * c
    ye_transcribed = nu**2 * ref_xia["A11"] + S / eps * E * t * (1 / (1 - nu**2) - (1 - nu**2) / (4 * (1 + nu) ** 2))
    print(f"\nYe(2014) Eq.(10) transcription of Xia A22 = {ye_transcribed:.6g}  (Xia original Table 2 → {ref_xia['A22']:.6g})")

    # 원호 주름 I1, I2 폐형식 vs 수치적분
    R, L = 0.003, 0.003
    c_r, l_r, I1c, I2c = round_corrugation_closed(R, L)
    I1n, I2n = round_corrugation_numeric(R, L)
    print(f"\nround corrugation R=L=3mm: c={c_r:.4f} l={l_r:.6f}  I1 closed={I1c:.6e} numeric={I1n:.6e}  I2 closed={I2c:.6e} numeric={I2n:.6e}")
