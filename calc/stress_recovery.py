"""등가 판 합력/변형률 → 주름 원판 국부 응력 복원 (B-4).

두 계열을 제공한다.
[1] Briassoulis (1986) Appendix B 능선 폐형식 (교차검증용, 단축 케이스)  — kt_x_tension, sigma_x_ridge, sigma_y_curvature
[2] Ye, Berdichevsky & Yu (2014) VAM 복원식 Eq. (14)–(17), (20) (채택, 6개 거시 변형률 동시)  — ye_recovery

국부축: x ⊥ 능선(파형 진행), y ∥ 능선.  기호는 docs/05, 06 (f = H/2, h = t, λ = l/c, J1, J2).
Ye 기호 → 물리량:  ε = p (주기), X = x/p,  a = 1 + z'²,  √a = ds/dx,  φ' = ε z'',  κ_shell = z''/a^{3/2} (단면 곡률)
    φ'² h²/(48 ε² a³) = h² z''²/(48 a³) = h² κ_shell²/48        h² φ' κ_xy/(12 ε a) = h² z'' κ_xy/(12 a)        φ'/(2εa) = z''/(2a)

복원식 (대칭 단면, ℬ = α2 = 0). 거시 변형률 m = {ε_xx, ε_yy, γ_xy(=2ε_xy), κ_xx, κ_yy, κ_xy}  (κ_αβ = +w,αβ):
    𝒞 = −(12 J2/h² + J1),  c1 = −(ε_xx + ν ε_yy)/𝒞,  c4 = (κ_xx + ν κ_yy)/λ,  α1 = 1/⟨√a/(1+h²κ²/48)⟩,  c2 = α1 γ_xy
    γ11 = c1 √a − ν a (ε_yy + z κ_yy)                      2γ12 = (√a c2 − h² z'' κ_xy/(12a)) / (1 + h²κ²/48)
    γ22 = ε_yy + z κ_yy
    ρ11 = a (12 c1 z/h² + c4) − ν √a κ_yy                  2ρ12 = (−2√a κ_xy + z'' c2/(2a)) / (1 + h²κ²/48)
    ρ22 = κ_yy/√a
  (ρ11 의 −ν√a κ_yy 와 ρ22 의 + 부호는 논문 Eq. (73), (91), (93) 유도와 평판 극한·에너지 일치 검사로 확정. 인쇄본 Eq. (14)
   의 부호와 다르므로 docs/09 §4 에 근거를 기록함.)
물리 성분 (공변 → 물리): ε_s = γ11/a, ε_y = γ22, γ_sy = 2γ12/√a, κ_s = ρ11/a, κ_y = ρ22, 2κ_sy = 2ρ12/√a
표면 변형률: e(ζ) = ε − ζ κ  (ζ: 법선 n 방향 두께 좌표, n 은 산에서 +z).  평면응력 → σ_s, σ_y, τ_sy, von Mises.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, Tuple

import numpy as np

from geometry import CorrugationGeometry, arc_tangent, sinusoidal, trapezoidal
from stiffness import Sheet, adopted, ye_full


# ============================================================================= [1] Briassoulis 능선 폐형식
@dataclass
class RidgeStress:
    inner: float
    outer: float
    nominal: float


def kt_x_tension(f: float, t: float) -> float:
    """능선 내측 섬유 응력집중계수 (Briassoulis Eq. B2): K_t = 1 + 6 f/t."""
    return 1.0 + 6.0 * f / t


def sigma_x_ridge(N_x: float, f: float, t: float) -> RidgeStress:
    s0 = N_x / t
    return RidgeStress(inner=s0 * (1 + 6 * f / t), outer=s0 * (1 - 6 * f / t), nominal=s0)


def sigma_y_curvature(M_y: float, f: float, t: float, mu: float, z: float) -> tuple[float, float]:
    """Briassoulis Eq. (B5): 능선 방향 굽힘 M_y* 상태, 높이 z 의 (외측, 내측) 섬유 응력."""
    k = 6 * M_y / t**2 / ((f / t) ** 2 * 6 * (1 - mu**2) + 1)
    a = (z / t) * 2 * (1 - mu**2)
    return (k * (-a - 1), k * (-a + 1))


def uniform_moment_and_axial(M_y: float, f: float, t: float) -> tuple[float, float]:
    A, I = t, t**3 / 12
    return M_y / (1 + A * f**2 / (2 * I)), M_y * f / (I / A - f**2 / 2)


# ============================================================================= [2] Ye (2014) VAM 복원
def profile_samples(g: CorrugationGeometry, n: int = 800) -> Dict[str, np.ndarray]:
    """한 주기 중앙면을 호 길이 균등 분할한 표본: x, z, zp(=dz/dx), kappa(단면 곡률, 산에서 음), ds, dx, sqrt_a, a."""
    p, f = g.p, g.f
    if g.profile in ("arc_tangent", "round"):
        if g.profile == "round":
            Rc = Rv = g.extra["R"]; alpha = math.pi / 2; TL = 2 * g.extra["L"]
        else:
            Rc, Rv, alpha, TL = g.extra["R_c"], g.extra["R_v"], math.radians(g.extra["alpha_deg"]), g.extra["T_L"]
        segs = []  # (length, fn(s_local)->(x,z,theta,kappa)) for the half period crest→valley
        segs.append((Rc * alpha, lambda s: (Rc * math.sin(s / Rc), f - Rc + Rc * math.cos(s / Rc), -s / Rc, -1 / Rc)))
        x1, z1 = Rc * math.sin(alpha), f - Rc + Rc * math.cos(alpha)
        segs.append((TL, lambda s: (x1 + s * math.cos(alpha), z1 - s * math.sin(alpha), -alpha, 0.0)))
        x2, z2 = x1 + TL * math.cos(alpha), z1 - TL * math.sin(alpha)
        segs.append((Rv * alpha, lambda s: (p / 2 - Rv * math.sin(alpha - s / Rv), -f + Rv - Rv * math.cos(alpha - s / Rv), -(alpha - s / Rv), +1 / Rv)))
        L = sum(s[0] for s in segs); nh = max(n // 2, 50)
        s_mid = (np.arange(nh) + 0.5) * L / nh
        xs, zs, ths, kps = [], [], [], []
        for s in s_mid:
            acc = 0.0
            for length, fn in segs:
                if s <= acc + length or (length, fn) is segs[-1]:
                    x, z, th, kp = fn(min(max(s - acc, 0.0), length)); break
                acc += length
            xs.append(x); zs.append(z); ths.append(th); kps.append(kp)
        xs, zs, ths, kps = map(np.array, (xs, zs, ths, kps))
        ds = np.full(nh, L / nh)
        # 전주기: x=p/2 거울 → z 대칭, 기울기 반대, 곡률 대칭
        x = np.concatenate([xs, p - xs[::-1]]); z = np.concatenate([zs, zs[::-1]])
        zp = np.concatenate([np.tan(ths), -np.tan(ths[::-1])]); kap = np.concatenate([kps, kps[::-1]])
        ds = np.concatenate([ds, ds[::-1]])
    elif g.profile == "sinusoidal":
        c = g.c
        xe = np.linspace(0, p, n + 1); x = 0.5 * (xe[1:] + xe[:-1])
        z = f * np.cos(np.pi * x / c); zp = -f * np.pi / c * np.sin(np.pi * x / c); zpp = -f * (np.pi / c) ** 2 * np.cos(np.pi * x / c)
        sa = np.sqrt(1 + zp**2); ds = sa * np.diff(xe); kap = zpp / sa**3
    elif g.profile == "trapezoidal":
        c, a = g.c, math.radians(g.extra["alpha_deg"]); w = g.extra["w_flat"]
        # 산 평탄(w/2) - 경사 - 골 평탄(w/2) 반주기, 모서리 곡률 무시
        pts = [(0, f), (w / 2, f), (w / 2 + 2 * f / math.tan(a), -f), (c, -f)]
        xs, zs, zps, dss = [], [], [], []
        for (xa, za), (xb, zb) in zip(pts[:-1], pts[1:]):
            L = math.hypot(xb - xa, zb - za); k = max(int(n / 6), 20)
            for i in range(k):
                t_ = (i + 0.5) / k
                xs.append(xa + t_ * (xb - xa)); zs.append(za + t_ * (zb - za)); zps.append((zb - za) / (xb - xa) if xb != xa else 0.0); dss.append(L / k)
        xs, zs, zps, dss = map(np.array, (xs, zs, zps, dss))
        x = np.concatenate([xs, p - xs[::-1]]); z = np.concatenate([zs, zs[::-1]]); zp = np.concatenate([zps, -zps[::-1]])
        ds = np.concatenate([dss, dss[::-1]]); kap = np.zeros_like(x)
    else:
        raise ValueError(g.profile)
    sa = np.sqrt(1 + zp**2)
    return dict(x=x, z=z, zp=zp, kappa=kap, ds=ds, dx=ds / sa, sqrt_a=sa, a=sa**2)


@dataclass
class RecoveryResult:
    x: np.ndarray; z: np.ndarray
    eps_s: np.ndarray; eps_y: np.ndarray; gam_sy: np.ndarray        # 물리 막변형률
    kap_s: np.ndarray; kap_y: np.ndarray; kap_sy2: np.ndarray        # 물리 곡률 (kap_sy2 = 2κ_sy)
    sig_s: np.ndarray; sig_y: np.ndarray; tau: np.ndarray            # (2, n): [0]=외측(ζ=+h/2, 법선 n 쪽), [1]=내측
    vm: np.ndarray                                                   # (2, n) von Mises
    N_s: np.ndarray; N_y: np.ndarray; M_s: np.ndarray; M_y: np.ndarray  # 국부 합력 (단위 호 길이당)

    def max_vm(self) -> Tuple[float, int, int]:
        i = np.unravel_index(np.argmax(self.vm), self.vm.shape)
        return float(self.vm[i]), int(i[0]), int(i[1])


def ye_recovery(g: CorrugationGeometry, m: Sheet, macro: Dict[str, float], n: int = 800) -> RecoveryResult:
    """VAM 복원 (대칭 단면). macro 키: exx, eyy, gxy, kxx, kyy, kxy (없으면 0)."""
    E, nu, G, h = m.E, m.nu, m.G, g.h
    exx, eyy, gxy = macro.get("exx", 0.0), macro.get("eyy", 0.0), macro.get("gxy", 0.0)
    kxx, kyy, kxy = macro.get("kxx", 0.0), macro.get("kyy", 0.0), macro.get("kxy", 0.0)
    S = profile_samples(g, n)
    z, zp, kap, sa, a, dx = S["z"], S["zp"], S["kappa"], S["sqrt_a"], S["a"], S["dx"]
    p = g.p
    avg = lambda q: np.sum(q * dx) / p
    lam, J1, J2 = avg(sa), avg(1 / sa), avg(z**2 * sa)
    corr = 1 + h**2 * kap**2 / 48
    C = -(12 * J2 / h**2 + J1)
    c1 = -(exx + nu * eyy) / C
    c4 = (kxx + nu * kyy) / lam
    alpha1 = 1 / avg(sa / corr)
    c2 = alpha1 * gxy
    zpp = kap * sa**3                          # z'' = κ a^{3/2}
    g11 = c1 * sa - nu * a * (eyy + z * kyy)
    g12x2 = (sa * c2 - h**2 * zpp * kxy / (12 * a)) / corr
    g22 = eyy + z * kyy
    r11 = a * (12 * c1 * z / h**2 + c4) - nu * sa * kyy
    r12x2 = (-2 * sa * kxy + zpp * c2 / (2 * a)) / corr
    r22 = kyy / sa
    eps_s, eps_y, gam = g11 / a, g22, g12x2 / sa
    kap_s, kap_y, kap_sy2 = r11 / a, r22, r12x2 / sa
    Q = E / (1 - nu**2)
    N_s, N_y = Q * h * (eps_s + nu * eps_y), Q * h * (eps_y + nu * eps_s)
    M_s, M_y = Q * h**3 / 12 * (kap_s + nu * kap_y), Q * h**3 / 12 * (kap_y + nu * kap_s)
    sig_s, sig_y, tau, vm = (np.zeros((2, len(z))) for _ in range(4))
    for i, zeta in enumerate((+h / 2, -h / 2)):
        es, ey, gs = eps_s - zeta * kap_s, eps_y - zeta * kap_y, gam - zeta * kap_sy2
        sig_s[i], sig_y[i], tau[i] = Q * (es + nu * ey), Q * (ey + nu * es), G * gs
        vm[i] = np.sqrt(sig_s[i] ** 2 - sig_s[i] * sig_y[i] + sig_y[i] ** 2 + 3 * tau[i] ** 2)
    return RecoveryResult(S["x"], z, eps_s, eps_y, gam, kap_s, kap_y, kap_sy2, sig_s, sig_y, tau, vm, N_s, N_y, M_s, M_y)


def cell_energy(g: CorrugationGeometry, m: Sheet, r: RecoveryResult) -> float:
    """복원 국부 변형률로 계산한 단위 투영면적당 변형에너지 (한 주기 평균)."""
    E, nu, G, h = m.E, m.nu, m.G, g.h
    S = profile_samples(g, len(r.x))
    Q = E / (1 - nu**2)
    u = 0.5 * (Q * h * (r.eps_s**2 + 2 * nu * r.eps_s * r.eps_y + r.eps_y**2) + G * h * r.gam_sy**2
               + Q * h**3 / 12 * (r.kap_s**2 + 2 * nu * r.kap_s * r.kap_y + r.kap_y**2) + G * h**3 / 12 * r.kap_sy2**2)
    return float(np.sum(u * S["ds"]) / g.p)


def macro_energy(K: Dict[str, float], macro: Dict[str, float]) -> float:
    e = np.array([macro.get("exx", 0), macro.get("eyy", 0), macro.get("gxy", 0)])
    k = np.array([macro.get("kxx", 0), macro.get("kyy", 0), 2 * macro.get("kxy", 0)])
    A = np.array([[K["A11"], K["A12"], 0], [K["A12"], K["A22"], 0], [0, 0, K["A66"]]])
    D = np.array([[K["D11"], K["D12"], 0], [K["D12"], K["D22"], 0], [0, 0, K["D66"]]])
    return 0.5 * (e @ A @ e + k @ D @ k)


def recovery_from_resultants(g: CorrugationGeometry, m: Sheet, N_loc: np.ndarray, M_loc: np.ndarray, n: int = 800) -> RecoveryResult:
    """국부축 합력 {N_x,N_y,N_xy}, {M_x,M_y,M_xy} → 거시 변형률(채택식 역행렬) → 복원."""
    K = adopted(g, m)
    A = np.array([[K["A11"], K["A12"], 0], [K["A12"], K["A22"], 0], [0, 0, K["A66"]]])
    D = np.array([[K["D11"], K["D12"], 0], [K["D12"], K["D22"], 0], [0, 0, K["D66"]]])
    e = np.linalg.solve(A, np.asarray(N_loc, float)); k = np.linalg.solve(D, np.asarray(M_loc, float))
    return ye_recovery(g, m, dict(exx=e[0], eyy=e[1], gxy=e[2], kxx=k[0], kyy=k[1], kxy=k[2] / 2), n)


# ============================================================================= 자체 검증
def _selfcheck() -> None:
    ok = True

    def chk(name, got, ref, tol):
        nonlocal ok
        rel = abs(got - ref) / max(abs(ref), 1e-30)
        ok &= rel < tol
        print(f"  [{'OK ' if rel < tol else 'BAD'}] {name}: {got:.6g} vs {ref:.6g} (rel {rel:.1e})")

    # [1] Briassoulis 폐형식 기본 검사
    assert abs(kt_x_tension(8.0, 1.0) - 49.0) < 1e-12
    o, i = sigma_y_curvature(1.0, 0.0, 1.0, 0.3, 0.0); assert abs(o + 6) < 1e-12 and abs(i - 6) < 1e-12

    g = arc_tangent(9.0, 3.2, 0.6, R_c=1.5); m = Sheet(193000.0, 0.3)
    Kf = ye_full(g, m); f, h = g.f, g.h
    print("[A] 평판 극한 (H→0): 복원 변형률 = 고전 판 변형률")
    gf = sinusoidal(10.0, 1e-7, 1.0)
    r = ye_recovery(gf, m, dict(exx=1e-3, eyy=-2e-4, gxy=3e-4, kxx=1e-4, kyy=-5e-5, kxy=2e-5))
    chk("eps_s", r.eps_s.mean(), 1e-3, 1e-6); chk("eps_y", r.eps_y.mean(), -2e-4, 1e-6); chk("gam_sy", r.gam_sy.mean(), 3e-4, 1e-6)
    chk("kap_s", r.kap_s.mean(), 1e-4, 1e-6); chk("kap_y", r.kap_y.mean(), -5e-5, 1e-6); chk("2kap_sy", r.kap_sy2.mean(), -2 * 2e-5, 1e-6)

    print("[B] Case 1 (ε_xx): 능선 N_s = A11 ε_xx, M_s = N_x f, K_t = 1+6f/h (Xia / Briassoulis)")
    r = ye_recovery(g, m, dict(exx=1e-3)); ic = int(np.argmax(r.z))
    Nx = Kf["A11"] * 1e-3
    chk("N_s(ridge)", r.N_s[ic], Nx, 2e-3); chk("M_s(ridge)", r.M_s[ic], Nx * f, 2e-3)
    chk("K_t inner", r.sig_s[1, ic] / (Nx / h), 1 + 6 * f / h, 2e-3)

    print("[C] Case 5 (κ_yy): 능선 N_y = E h f κ_yy (Briassoulis Eq. 12) — Xia 가정은 1/(1-ν²) 과대")
    r = ye_recovery(g, m, dict(kyy=1e-3)); ic = int(np.argmax(r.z))
    chk("N_y(ridge)", r.N_y[ic], m.E * h * f * 1e-3, 2e-3)

    print("[D] 변형에너지 일치: 복원 국부 변형률의 셀 에너지 == ½ mᵀ[A,D]m (ye_full)")
    for name, mac in (("exx", dict(exx=1e-3)), ("eyy", dict(eyy=1e-3)), ("gxy", dict(gxy=1e-3)), ("kxx", dict(kxx=1e-3)),
                      ("kyy", dict(kyy=1e-3)), ("kxy", dict(kxy=1e-3)), ("exx+eyy", dict(exx=1e-3, eyy=5e-4)),
                      ("kxx+kyy", dict(kxx=1e-3, kyy=-4e-4)), ("all", dict(exx=1e-3, eyy=-3e-4, gxy=4e-4, kxx=2e-4, kyy=-1e-4, kxy=3e-4))):
        r = ye_recovery(g, m, mac)
        tol = 5e-3 if ("gxy" in mac or "kxy" in mac) else 2e-3
        chk(f"energy {name}", cell_energy(g, m, r), macro_energy(Kf, mac), tol)

    print("[E] Briassoulis Fig. 6 경향: 순수 κ_yy 에서 능선 |σ_y|max / (6 M_y*/t²) 는 얇은 판(f/t 大)에서 1 미만, 두꺼운 판에서 ~1.13")
    for ft in (0.25, 0.5, 1.0, 2.0, 4.0):
        gg = arc_tangent(9.0, 2 * ft * 0.6, 0.6, R_c=1.2) if ft * 0.6 * 2 < 9 else None
        rr = ye_recovery(gg, m, dict(kyy=1e-3)); Kk = ye_full(gg, m)
        ratio = np.max(np.abs(rr.sig_y)) / (6 * Kk["D22"] * 1e-3 / gg.h**2)
        print(f"     f/t={ft:4.2f}: ratio={ratio:.3f}")

    print("[F] 합력 → 복원 왕복: N,M 으로 복원한 거시 변형률이 원래와 같음")
    mac = dict(exx=1e-3, eyy=-3e-4, gxy=4e-4, kxx=2e-4, kyy=-1e-4, kxy=3e-4)
    K = adopted(g, m)
    e = np.array([mac["exx"], mac["eyy"], mac["gxy"]]); k = np.array([mac["kxx"], mac["kyy"], 2 * mac["kxy"]])
    A = np.array([[K["A11"], K["A12"], 0], [K["A12"], K["A22"], 0], [0, 0, K["A66"]]]); D = np.array([[K["D11"], K["D12"], 0], [K["D12"], K["D22"], 0], [0, 0, K["D66"]]])
    r1, r2 = ye_recovery(g, m, mac), recovery_from_resultants(g, m, A @ e, D @ k)
    chk("roundtrip max vm", r2.max_vm()[0], r1.max_vm()[0], 1e-9)
    print("ALL OK" if ok else "SOME CHECKS FAILED")


if __name__ == "__main__":
    _selfcheck()
    g = arc_tangent(9.0, 3.2, 0.6, R_c=1.5); m = Sheet(193000.0, 0.3)
    print("\n=== 예시: 국부축 합력 N_x=10 N/mm, M_y=100 N·mm/mm (가정 기본값 판) ===")
    r = recovery_from_resultants(g, m, [10.0, 0, 0], [0, 100.0, 0])
    vmax, side, idx = r.max_vm()
    print(f"  max von Mises = {vmax:.1f} MPa at x={r.x[idx]:.2f} mm (z={r.z[idx]:.2f}), {'outer' if side == 0 else 'inner'} fibre")
    ic = int(np.argmax(r.z)); iflank = int(np.argmin(np.abs(r.z)))
    print(f"  ridge (z=+f): sig_s outer/inner = {r.sig_s[0,ic]:.1f}/{r.sig_s[1,ic]:.1f}, sig_y outer/inner = {r.sig_y[0,ic]:.1f}/{r.sig_y[1,ic]:.1f} MPa")
    print(f"  flank (z≈0):  sig_s outer/inner = {r.sig_s[0,iflank]:.1f}/{r.sig_s[1,iflank]:.1f}, sig_y = {r.sig_y[0,iflank]:.1f}/{r.sig_y[1,iflank]:.1f}, tau = {r.tau[0,iflank]:.1f} MPa")
