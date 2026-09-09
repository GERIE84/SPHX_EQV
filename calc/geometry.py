"""주름 단면 형상 라이브러리 — 등가모델 입력용 파생 형상량 계산.

docs/05_geometry_parameters.md 의 정의를 구현한다.

1차 입력 (docs/04_plate_dimension_form.md):
    p   피치 (능선 수직 방향, 산-산)          [mm]
    H   깊이 (중앙면 산-골)                    [mm]
    t   판 두께                                [mm]
    단면형식별 추가: R_c, R_v (원호 반경), alpha (플랭크 각) 또는 T_L (접선 길이), w_c/w_v (사다리꼴 평탄폭)

파생 형상량 (문헌 기호):
    c = p/2, f = H/2, h = t
    l   반주기 전개길이                      = S/2 (Ye), 반주기 호 길이
    I1  = ∫_0^{2l} (dx/ds)^2 ds              (한 주기)
    I2  = ∫_0^{2l} z^2 ds                    (한 주기, 중앙면 기준, <z>=0)
    Ye 셀 평균: <sqrt a> = l/c, <1/sqrt a> = I1/(2c), eps^2 <phi^2 sqrt a> = I2/(2c)  (eps = 2c)

모든 단면은 산(x=0, z=+H/2)에서 시작해 골(x=p/2, z=-H/2)로 내려가는 반주기를 정의하고,
전주기는 x=p/2 에 대한 거울 대칭으로 만든다. 따라서 z 의 주기 평균은 0 이다.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable, Dict, Optional

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq


# ----------------------------------------------------------------------------- 결과 컨테이너
@dataclass
class CorrugationGeometry:
    profile: str
    p: float
    H: float
    t: float
    l: float          # 반주기 전개길이
    I1: float         # ∫(dx/ds)^2 ds  (한 주기)
    I2: float         # ∫ z^2 ds       (한 주기)
    extra: Dict[str, float] = field(default_factory=dict)

    # --- 문헌 기호 ---
    @property
    def c(self) -> float:
        return self.p / 2

    @property
    def f(self) -> float:
        return self.H / 2

    @property
    def h(self) -> float:
        return self.t

    # --- Ye (2014) 셀 평균 ---
    @property
    def avg_sqrt_a(self) -> float:          # <√a> = S/ε = l/c
        return self.l / self.c

    @property
    def avg_inv_sqrt_a(self) -> float:      # <1/√a> = I1/(2c)
        return self.I1 / (2 * self.c)

    @property
    def eps2_avg_phi2_sqrt_a(self) -> float:  # ε²<φ²√a> = I2/(2c)  [length^2]
        return self.I2 / (2 * self.c)

    # --- 단위폭 단면 특성 (Lang & Su Eq. 28-30) ---
    @property
    def A_u(self) -> float:                 # 단위폭당 단면적 = h l / c
        return self.t * self.l / self.c

    @property
    def I_u(self) -> float:                 # 단위폭당 2차 모멘트 = h I2 / (2c)
        return self.t * self.I2 / (2 * self.c)

    @property
    def W_u(self) -> float:                 # 단위폭당 단면계수 = I_u / (f + h/2)
        return self.I_u / (self.f + self.t / 2)

    # --- 무차원 군 ---
    def dimensionless(self) -> Dict[str, float]:
        d = dict(H_over_p=self.H / self.p, t_over_p=self.t / self.p, t_over_H=self.t / self.H,
                 l_over_c=self.l / self.c)
        if "R_min" in self.extra:
            d["t_over_Rmin"] = self.t / self.extra["R_min"]
        return d

    def summary(self) -> str:
        lines = [f"profile={self.profile}  p={self.p:g}  H={self.H:g}  t={self.t:g}",
                 f"  c={self.c:g}  f={self.f:g}  l={self.l:.6g}  I1={self.I1:.6g}  I2={self.I2:.6g}",
                 f"  <sqrt a>={self.avg_sqrt_a:.6g}  <1/sqrt a>={self.avg_inv_sqrt_a:.6g}  I2/(2c)={self.eps2_avg_phi2_sqrt_a:.6g}",
                 f"  A_u={self.A_u:.6g}  I_u={self.I_u:.6g}  W_u={self.W_u:.6g}"]
        if self.extra:
            lines.append("  extra: " + ", ".join(f"{k}={v:.6g}" for k, v in self.extra.items()))
        return "\n".join(lines)


# ----------------------------------------------------------------------------- 공통: 수치 적분 (임의 z(x))
def integrals_from_zx(z: Callable[[float], float], dz: Callable[[float], float], p: float) -> tuple[float, float, float]:
    """z(x), z'(x) 가 주어진 임의 단면의 l, I1, I2 (한 주기 [0, p])."""
    ds = lambda x: math.sqrt(1 + dz(x) ** 2)
    l = quad(ds, 0, p / 2, limit=200)[0]
    I1 = quad(lambda x: 1 / ds(x), 0, p, limit=200)[0]           # (dx/ds)^2 ds = dx/ds dx = dx/√a
    I2 = quad(lambda x: z(x) ** 2 * ds(x), 0, p, limit=200)[0]
    return l, I1, I2


# ----------------------------------------------------------------------------- 단면 1: 사인형
def sinusoidal(p: float, H: float, t: float) -> CorrugationGeometry:
    """z = f sin(πx/c) 와 동일 (산 위치만 이동). 폐형식 없음 → 수치적분."""
    f, c = H / 2, p / 2
    z = lambda x: f * math.cos(math.pi * x / c)
    dz = lambda x: -f * math.pi / c * math.sin(math.pi * x / c)
    l, I1, I2 = integrals_from_zx(z, dz, p)
    R_min = c**2 / (math.pi**2 * f)  # 산·골에서의 곡률반경 1/|z''|
    return CorrugationGeometry("sinusoidal", p, H, t, l, I1, I2, extra=dict(R_min=R_min))


# ----------------------------------------------------------------------------- 단면 2: 원호 + 접선 (arc-and-tangent)
def _arc_tangent_solve(p: float, H: float, R_c: float, R_v: float) -> tuple[float, float]:
    """R_c, R_v 가 주어질 때 플랭크 각 alpha [rad], 접선 길이 T_L 을 구한다.
    수평: (R_c+R_v) sin α + T_L cos α = p/2 ;  수직: (R_c+R_v)(1-cos α) + T_L sin α = H
    """
    Rs = R_c + R_v
    # 특수해: 반원 + 수직 직선 (round) → α = 90°, T_L = H - Rs
    if abs(p / 2 - Rs) <= 1e-9 * max(p, 1.0) and H >= Rs - 1e-12:
        return math.pi / 2, H - Rs
    # T_L ≥ 0 이 되는 α 의 상한: 수평식 → sin α ≤ p/(2Rs), 수직식 → 1 - cos α ≤ H/Rs
    a_h = math.asin(min(1.0, p / (2 * Rs)))
    a_v = math.acos(max(-1.0, 1 - H / Rs))
    a_hi = min(a_h, a_v, math.pi / 2) - 1e-9
    T_h = lambda a: (p / 2 - Rs * math.sin(a)) / math.cos(a)
    T_v = lambda a: (H - Rs * (1 - math.cos(a))) / math.sin(a)
    g = lambda a: T_h(a) - T_v(a)          # g(0+) → -∞ ; 유효 근에서 부호 변화
    if g(a_hi) <= 0:
        raise ValueError("arc-tangent: 주어진 p, H, R_c, R_v 로 접선 구간(T_L ≥ 0)을 만들 수 없음 — 원호 반경이 너무 큼")
    alpha = brentq(g, 1e-6, a_hi, xtol=1e-14)
    T_L = T_h(alpha) if math.cos(alpha) >= math.sin(alpha) else T_v(alpha)
    if T_L < -1e-9:
        raise ValueError("arc-tangent: T_L < 0")
    return alpha, max(T_L, 0.0)


def _arc_tangent_from_alpha(p: float, H: float, alpha: float) -> tuple[float, float]:
    """R_c = R_v = R 가정, alpha 주어질 때 R, T_L (선형 2식)."""
    s, cs = math.sin(alpha), math.cos(alpha)
    # 2R s + T cs = p/2 ; 2R(1-cs) + T s = H
    A = np.array([[2 * s, cs], [2 * (1 - cs), s]])
    R, T = np.linalg.solve(A, np.array([p / 2, H]))
    if R <= 0 or T < -1e-9:
        raise ValueError("arc-tangent(alpha): 해가 물리적이지 않음 (R<=0 또는 T_L<0)")
    return R, max(T, 0.0)


def arc_tangent(p: float, H: float, t: float, R_c: Optional[float] = None, R_v: Optional[float] = None,
                alpha_deg: Optional[float] = None) -> CorrugationGeometry:
    """원호+접선 단면. 입력 모드: (R_c, R_v) 또는 (alpha_deg, R_c=R_v).  폐형식 l, I1, I2 (Lang & Su Eq. 27 일반화)."""
    if alpha_deg is not None and R_c is None:
        alpha = math.radians(alpha_deg)
        R, T_L = _arc_tangent_from_alpha(p, H, alpha)
        R_c = R_v = R
    else:
        if R_v is None:
            R_v = R_c
        alpha, T_L = _arc_tangent_solve(p, H, R_c, R_v)
    l = (R_c + R_v) * alpha + T_L
    q = alpha / 2 + math.sin(2 * alpha) / 4               # ∫_0^α cos²θ dθ
    I1 = 2 * ((R_c + R_v) * q + T_L * math.cos(alpha) ** 2)
    a_c, a_v = H / 2 - R_c, H / 2 - R_v                   # 원호 중심의 |z|
    I2_crest = R_c * (a_c**2 * alpha + 2 * a_c * R_c * math.sin(alpha) + R_c**2 * q)
    I2_valley = R_v * (a_v**2 * alpha + 2 * a_v * R_v * math.sin(alpha) + R_v**2 * q)
    z1 = H / 2 - R_c * (1 - math.cos(alpha))              # 플랭크 시작 z
    z2 = z1 - T_L * math.sin(alpha)                       # 플랭크 끝 z (= -H/2 + R_v(1-cos α))
    I2_flank = (z1**3 - z2**3) / (3 * math.sin(alpha)) if T_L > 0 else 0.0
    I2 = 2 * (I2_crest + I2_flank + I2_valley)
    return CorrugationGeometry("arc_tangent", p, H, t, l, I1, I2,
                               extra=dict(R_c=R_c, R_v=R_v, alpha_deg=math.degrees(alpha), T_L=T_L, R_min=min(R_c, R_v)))


def arc_tangent_profile_xz(g: CorrugationGeometry, n: int = 400) -> tuple[np.ndarray, np.ndarray]:
    """작도·FE 형상 생성용 중앙면 좌표 (한 주기)."""
    R_c, R_v, alpha, T_L, H, p = g.extra["R_c"], g.extra["R_v"], math.radians(g.extra["alpha_deg"]), g.extra["T_L"], g.H, g.p
    th = np.linspace(0, alpha, n // 4)
    xa, za = R_c * np.sin(th), H / 2 - R_c + R_c * np.cos(th)
    s = np.linspace(0, T_L, max(n // 4, 2))
    xf, zf = xa[-1] + s * np.cos(alpha), za[-1] - s * np.sin(alpha)
    th2 = np.linspace(alpha, 0, n // 4)
    xv, zv = p / 2 - R_v * np.sin(th2), -H / 2 + R_v - R_v * np.cos(th2)
    xh = np.concatenate([xa, xf[1:], xv[1:]]); zh = np.concatenate([za, zf[1:], zv[1:]])
    return np.concatenate([xh, p - xh[::-1][1:]]), np.concatenate([zh, zh[::-1][1:]])


# ----------------------------------------------------------------------------- 단면 3: 사다리꼴 (모서리 날카로움)
def trapezoidal(p: float, H: float, t: float, alpha_deg: float) -> CorrugationGeometry:
    """Lang & Su Eq. (15) / Xia Table 3.  산·골 평탄폭은 w = c - H/tan α (자동)."""
    c, f, a = p / 2, H / 2, math.radians(alpha_deg)
    w_flat = c - 2 * f / math.tan(a)
    if w_flat < -1e-9:
        raise ValueError("trapezoidal: 평탄폭 < 0 (알파가 너무 작음)")
    l = 2 * f / math.sin(a) + w_flat
    I1 = 4 * f * math.cos(a) ** 2 / math.sin(a) + 2 * w_flat
    I2 = 4 * f**3 / (3 * math.sin(a)) + 2 * f**2 * w_flat
    return CorrugationGeometry("trapezoidal", p, H, t, l, I1, I2, extra=dict(alpha_deg=alpha_deg, w_flat=w_flat))


# ----------------------------------------------------------------------------- 단면 4: 원호(round) — 반원 + 수직 직선 (Xia Table 4)
def round_corrugation(R: float, L: float, t: float) -> CorrugationGeometry:
    """Xia (2012) Table 4: c=2R, l=πR+2L, I1=πR, I2=4L³/3+2πL²R+8LR²+πR³.  p=4R, H=2R+2L."""
    p, H = 4 * R, 2 * R + 2 * L
    l = math.pi * R + 2 * L
    I1 = math.pi * R
    I2 = 4 * L**3 / 3 + 2 * math.pi * L**2 * R + 8 * L * R**2 + math.pi * R**3
    return CorrugationGeometry("round", p, H, t, l, I1, I2, extra=dict(R=R, L=L, R_min=R))


# ----------------------------------------------------------------------------- 팩토리
def from_dict(d: dict) -> CorrugationGeometry:
    """docs/04 양식(B절) 값을 담은 dict 로부터 형상 객체 생성.  키: profile, p, H, t, R_c, R_v, alpha_deg, R, L"""
    prof = d["profile"]
    if prof == "sinusoidal":
        return sinusoidal(d["p"], d["H"], d["t"])
    if prof == "arc_tangent":
        return arc_tangent(d["p"], d["H"], d["t"], R_c=d.get("R_c"), R_v=d.get("R_v"), alpha_deg=d.get("alpha_deg"))
    if prof == "trapezoidal":
        return trapezoidal(d["p"], d["H"], d["t"], d["alpha_deg"])
    if prof == "round":
        return round_corrugation(d["R"], d["L"], d["t"])
    raise ValueError(f"unknown profile {prof}")


# ----------------------------------------------------------------------------- 자체 검증
def _selfcheck() -> None:
    ok = True

    def chk(name, got, ref, tol=2e-4):
        nonlocal ok
        rel = abs(got - ref) / abs(ref)
        flag = "OK " if rel < tol else "BAD"
        ok &= rel < tol
        print(f"  [{flag}] {name}: {got:.6g} (ref {ref:.6g}, rel {rel:.1e})")

    print("1) 사인형 — Lang & Su Table 4 (c=0.32, f=0.11 m) → verify_langsu_tables.py 값")
    g = sinusoidal(p=0.64, H=0.22, t=0.005)
    chk("l", g.l, 0.39910); chk("I1", g.I1, 0.52265); chk("I2", g.I2, 0.004374, 3e-4)

    print("2) 사다리꼴 — Lang & Su Table 4 (c=0.0508, f=0.0127, α=45°)")
    g = trapezoidal(p=0.1016, H=0.0254, t=0.00635, alpha_deg=45)
    chk("l", g.l, 0.06132); chk("I1", g.I1, 0.08672); chk("I2", g.I2, 1.2056e-5)

    print("3) 원호(round) — Xia Table 4 폐형식 vs arc_tangent(α=90°) 일반식")
    R, L = 3.0, 3.0
    g1 = round_corrugation(R, L, t=0.2)
    g2 = arc_tangent(p=4 * R, H=2 * R + 2 * L, t=0.2, R_c=R, R_v=R)
    chk("alpha_deg", g2.extra["alpha_deg"], 90.0); chk("T_L", g2.extra["T_L"], 2 * L)
    chk("l", g2.l, g1.l); chk("I1", g2.I1, g1.I1); chk("I2", g2.I2, g1.I2)

    print("4) 원호-접선 — Lang & Su Table 10 (200×55×3 암거판: r=53, T_L=32.3, α0=45.187°, c=100, f=27.5)")
    g = arc_tangent(p=200.0, H=55.0, t=3.0, R_c=53 + 1.5, R_v=53 + 1.5)
    chk("alpha_deg", g.extra["alpha_deg"], 45.187, 2e-3)
    chk("T_L (공표값 세트 자체가 반올림되어 수평식 잔차 0.1 mm → 0.4% 허용)", g.extra["T_L"], 32.3, 6e-3)
    # Lang & Su Eq.(27) 와 비교
    R, a, TL = 54.5, math.radians(g.extra["alpha_deg"]), g.extra["T_L"]
    chk("l vs Eq.27", g.l, TL + 2 * R * a); chk("I1 vs Eq.27", g.I1, R * (2 * a + math.sin(2 * a)) + 2 * TL * math.cos(a) ** 2)

    print("5) 원호-접선 폐형식 I2 vs 수치적분 (비대칭 R_c≠R_v 포함)")
    g = arc_tangent(p=10.0, H=3.0, t=0.5, R_c=1.2, R_v=1.8)
    x, z = arc_tangent_profile_xz(g, n=4000)
    ds = np.hypot(np.diff(x), np.diff(z)); zm = 0.5 * (z[1:] + z[:-1])
    chk("l (numeric)", np.sum(ds) / 2, g.l, 1e-3); chk("I2 (numeric)", np.sum(zm**2 * ds), g.I2, 2e-3)
    cs = np.diff(x) / ds
    chk("I1 (numeric)", np.sum(cs**2 * ds), g.I1, 2e-3)
    print("ALL OK" if ok else "SOME CHECKS FAILED")


if __name__ == "__main__":
    _selfcheck()
    print("\n예시 (가정값, PHE 통상 범위): arc_tangent p=9, H=3.2, t=0.6, R_c=R_v=1.5 mm")
    print(arc_tangent(9.0, 3.2, 0.6, R_c=1.5).summary())
