"""C-4: 압력차 → 판 합력·응력 폐형식 (FEA 미사용 경로).  docs/12_pressure_closed_form.md

적층 스택에서 판 n 은 이웃 판과 능선 교차점(접촉점 격자)에서만 지지된다. 압력차 ΔP 가 걸리면 판은
접촉점 격자 위의 연속 판이 되고, 하중은 두 경로로 지지점에 전달된다.
  (A) 가로(주름 단면) 경로 — 지지선(고압면 반대쪽 능선) 사이 한 피치 폭의 주름 단면이 '곡선 프레임'으로
      압력을 받는다. 스팬이 정확히 1 피치이므로 균질화 판(D11)이 아니라 단면 자체의 곡선보 정역학으로 푼다
      (단위 폭 스트립, y 방향 평면변형, Castigliano 로 잉여력 해결). 출력: 단면 위치별 N_s, M_s → 표면 응력.
  (B) 능선 방향 경로 — 지지선 위의 접촉점 간격 a = p / sin 2β 를 스팬으로 하는 연속보. 단위 폭당 M_y =
      −q a²/12 (지지점, 고압면 인장), +q a²/24 (중앙). D22 ≫ D11 이므로 하중은 능선 방향으로 흐른다고 본다.
      VAM 복원(κ_yy = M_y/D22)으로 산·골 응력.
  (C) 접촉점 — 접촉력 F = ΔP·p²/sin 2β (한 격자 셀의 압력). 교차 원통 Hertz 로 최대 접촉압 p0 (1차 추정,
      국부 항복 검토용).
(A)+(B) 를 산·골 × +z/−z 면 × (지지점, 스팬 중앙) 8 점에서 중첩하여 최대 von Mises 를 구한다.

부호·좌표: 국부축 x ⊥ 능선, z 판 법선(+z = 산 쪽). 고압면 '+z' 이면 압력이 −z 로 밀어 골 선(−z 접촉)이 지지선,
'−z' 이면 산 선이 지지선(단면을 z 반전하여 같은 절차). 면내 구속: 'free'(판 가로 방향 자유, 거시 N_x = 0) 또는
'fixed'(가로 변위 0, 추력 P 발생). 표면 응력: σ(+n 면) = N/t − 6M/t², σ(−n 면) = N/t + 6M/t² (n = +z 쪽 법선).
검증: 평판 극한에서 M(중앙) = q p²/24, M(지지) = −q p²/12, N = 0 ; 수직 반력 = q c ; 회전 적합 ∫M/EI ds = 0.
한계: 두께 방향 응력 분포는 얇은 곡선보 가정 (t/R_min ≤ 0.2 밖이면 09 §5 와 같은 주의), 접촉 마찰·판 간 미끄럼
무시, 테두리·포트 근처 비주기 영역 제외.
"""
from __future__ import annotations

import math
import os
import sys
from dataclasses import dataclass
from typing import Dict, Optional

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from geometry import CorrugationGeometry, arc_tangent, sinusoidal   # noqa: E402
from stiffness import Sheet, adopted                               # noqa: E402
from stress_recovery import profile_samples, ye_recovery           # noqa: E402


# ============================================================================= (C) 접촉 격자
def contact_lattice(p: float, beta_deg: float) -> Dict[str, float]:
    """접촉점 격자: 두 판 능선군(±β)의 교차. 능선을 따라 간격 a = p/sin2β, 셀 면적 p²/sin2β."""
    s2b = math.sin(math.radians(2 * beta_deg))
    if s2b < 1e-9:
        raise ValueError("β = 0 또는 90°: 능선이 교차하지 않아 접촉 격자가 없음")
    return dict(a=p / s2b, A_cell=p**2 / s2b, crossing_deg=2 * beta_deg)


def hertz_crossed_cylinders(F: float, R_outer: float, crossing_deg: float, E: float, nu: float) -> Dict[str, float]:
    """같은 반경 R 두 원통이 각 φ 로 교차하는 Hertz 접촉 (Johnson §4.2, F1≈1 근사, 오차 <5 %).
    A=(1−|cosφ|)/2R, B=(1+|cosφ|)/2R, R_e = 1/(2√(AB)), E* = E/(2(1−ν²)) (같은 재료), c=√(ab)=(3FR_e/4E*)^{1/3}, p0=3F/(2πc²)."""
    phi = math.radians(crossing_deg); cphi = abs(math.cos(phi))
    A, B = (1 - cphi) / (2 * R_outer), (1 + cphi) / (2 * R_outer)
    Re = 1 / (2 * math.sqrt(A * B)); Es = E / (2 * (1 - nu**2))
    c = (3 * F * Re / (4 * Es)) ** (1 / 3)
    p0 = 3 * F / (2 * math.pi * c**2)
    return dict(F=F, R_outer=R_outer, crossing_deg=crossing_deg, R_e=Re, c_mean_semi_axis=c, p0=p0, p_mean=2 * p0 / 3,
                axis_ratio_BA=B / A, yield_onset_p0_over_Y=1.6)


# ============================================================================= (A) 가로 곡선 프레임
@dataclass
class FrameResult:
    s: np.ndarray; x: np.ndarray; z: np.ndarray; theta: np.ndarray; ds: np.ndarray
    N: np.ndarray; M: np.ndarray                       # 단위 폭당 [N/mm], [N·mm/mm]; M>0 = −n(−z) 면 인장
    sig_pos: np.ndarray; sig_neg: np.ndarray           # +n 면 / −n 면 σ_s [MPa]
    P: float; M_c: float; R_support: float; inplane: str
    i_crest: int; i_valley: int          # 표준 좌표: i_crest = 하중 쪽 능선(x=0), i_valley = 지지선(x=c)
    M_loaded: float = 0.0; N_loaded: float = 0.0; M_support: float = 0.0; N_support: float = 0.0   # 정확한 끝점 값 (표본 중점 아님)

    def at(self, i: int) -> Dict[str, float]:
        return dict(N=float(self.N[i]), M=float(self.M[i]), sig_pos=float(self.sig_pos[i]), sig_neg=float(self.sig_neg[i]))


def transverse_frame(g: CorrugationGeometry, m: Sheet, q: float, face: str = "+z", inplane: str = "free", n: int = 800) -> FrameResult:
    """한 피치 폭 주름 단면(단위 폭 스트립)이 압력 q [MPa] 를 받아 지지선(고압면 반대쪽 능선)으로 전달하는 곡선 프레임 해.
    face: 고압면 '+z' (지지 = 골 선) | '−z' (지지 = 산 선).  inplane: 'free' (P=0) | 'fixed' (가로 변위 0)."""
    S = profile_samples(g, n)
    half = S["x"] <= g.c + 1e-12
    x, z, zp, ds = S["x"][half], S["z"][half], S["zp"][half], S["ds"][half]
    if face in ("-z", "−z", "valley"):
        # 고압면이 −z: 압력은 +z 로 밀고 지지선은 산. 좌표를 z → −z, x → c − x 로 바꾸면(순서 반전) '압력이 −z 로 밀고
        # x=0 이 하중 쪽 능선(z=+f), x=c 가 지지선(z=−f)' 인 표준 문제가 된다. dz'/dx' = (−dz)/(−dx) = z' 그대로.
        x, z, zp, ds = (g.c - x)[::-1], (-z)[::-1], zp[::-1], ds[::-1]
    # 표준 문제: 압력 q 가 −z 로 밀고, x=0 (대칭선, 하중 쪽 능선) 에서 전단 0, x=c (대칭선, 지지선) 에서 수직 반력 q·c 가 나감.
    h, E, nu = g.h, m.E, m.nu
    EI = E * h**3 / (12 * (1 - nu**2)); EA = E * h / (1 - nu**2)
    th = np.arctan(zp); dx = ds * np.cos(th); dz = ds * np.sin(th)
    z0 = z[0] + 0.0                              # 시작 단면(x=0, 대칭선) 높이
    # 분포 압력 (단위 폭, x 구간 dx): F = q (z', −1) dx = q (dz, −dx)  [−z 로 미는 압력]
    Fx, Fz = q * dz, -q * dx
    cFx, cFz = np.cumsum(Fx) - 0.5 * Fx, np.cumsum(Fz) - 0.5 * Fz          # 표본 중점까지의 누적 (반 셀 보정)
    cxFz, czFx = np.cumsum(x * Fz) - 0.5 * x * Fz, np.cumsum(z * Fx) - 0.5 * z * Fx
    # 절단면 i 에 대한 하중 모멘트 Σ r_j × F_j = Σ[(x_j − x_i) F_zj − (z_j − z_i) F_xj]
    S_load = (cxFz - x * cFz) - (czFx - z * cFx)
    # 내력 = 선형 조합:  Q = −(P,0) − ΣF ;  M = −M_c + (z0 − z_i)·P·(−1)?  → 유도: M_i = −M_c − [r_c × (P,0)] − S_load,
    # r_c = (0 − x_i, z0 − z_i), (P,0): cross = r_x·0 − r_z·P = −(z0 − z_i) P  →  M_i = −M_c + (z0 − z_i) P − S_load
    t_vec_x, t_vec_z = np.cos(th), np.sin(th)
    N0 = (-cFx) * t_vec_x + (-cFz) * t_vec_z;  NP = -t_vec_x;              M0 = -S_load; MP = (z0 - z); MM = -np.ones_like(x)
    # Castigliano: ∂U/∂M_c = 0 (항상), ∂U/∂P = 0 (fixed 일 때만); U = Σ (M²/2EI + N²/2EA) ds
    w_M, w_N = ds / EI, ds / EA
    a_MM = np.sum(MM * MM * w_M); a_MP = np.sum(MM * MP * w_M); b_M = -np.sum(MM * M0 * w_M)
    if inplane == "fixed":
        a_PP = np.sum(MP * MP * w_M + NP * NP * w_N); b_P = -np.sum(MP * M0 * w_M + NP * N0 * w_N)
        Mc, P = np.linalg.solve(np.array([[a_MM, a_MP], [a_MP, a_PP]]), np.array([b_M, b_P]))
    else:
        P = 0.0; Mc = b_M / a_MM
    N = N0 + P * NP; M = M0 + P * MP + Mc * MM
    sig_pos, sig_neg = N / h - 6 * M / h**2, N / h + 6 * M / h**2
    s = np.cumsum(ds) - 0.5 * ds
    # 끝점 정확값: x=0 → S_load=0, MP=0 : M=−Mc, N=−P ;  x=c (θ=0, z=−z0) → 전체 합
    SFx, SFz = Fx.sum(), Fz.sum(); Sc = (np.sum(x * Fz) - g.c * SFz) - (np.sum(z * Fx) - (-z0) * SFx)
    M_sup = -Sc + P * (z0 - (-z0)) - Mc; N_sup = -SFx - P
    fr = FrameResult(s, x, z, th, ds, N, M, sig_pos, sig_neg, float(P), float(Mc), float(q * g.c), inplane, 0, len(x) - 1,
                     float(-Mc), float(-P), float(M_sup), float(N_sup))
    return fr


# ============================================================================= (B) 능선 방향 연속보
def ridge_beam_moments(q: float, p: float, beta_deg: float) -> Dict[str, float]:
    """지지선 위 접촉점 간격 a 의 연속보(단위 폭, 등분포 q): 지지점 M_y = −q a²/12 (고압면 인장), 중앙 +q a²/24."""
    a = contact_lattice(p, beta_deg)["a"]
    return dict(a=a, My_support=-q * a**2 / 12, My_mid=q * a**2 / 24)


# ============================================================================= 종합
@dataclass
class PressureCase:
    q: float; face: str; inplane: str; beta_deg: float
    lattice: Dict[str, float]; frame: FrameResult; ridge: Dict[str, float]; contact: Dict[str, float]
    table: list                                       # 8 점 중첩 결과 dict 목록
    max_vm: float; max_where: str
    frame_M_loaded: float; frame_M_support: float; frame_N_loaded: float; frame_N_support: float   # 하중 쪽 능선 / 지지선 (단위 폭)
    homogenized_Mx_check: Dict[str, float]


def pressure_case(g: CorrugationGeometry, m: Sheet, dP: float, beta_deg: float, face: str = "+z", inplane: str = "free",
                  n: int = 800) -> PressureCase:
    """압력차 dP [MPa] 에 대한 (A)+(B)+(C) 종합. face = 고압면 ('+z' 산 쪽 / '−z' 골 쪽)."""
    fr = transverse_frame(g, m, dP, face, inplane, n)
    lat = contact_lattice(g.p, beta_deg)
    rb = ridge_beam_moments(dP, g.p, beta_deg)
    K = adopted(g, m)
    # 반전 좌표계에서 fr.x=0 은 '하중 쪽 능선(자유)', fr.x=c 는 '지지선'. 물리 명칭으로 되돌림
    loaded_name, support_name = ("crest", "valley") if face in ("+z", "crest") else ("valley", "crest")
    zsign = 1.0 if face in ("+z", "crest") else -1.0   # 반전 좌표 +n 면 = 물리 +z 면 (face=+z) / 물리 −z 면 (face=−z)
    rows = []
    h = g.h
    ends = {loaded_name: dict(N=fr.N_loaded, M=fr.M_loaded), support_name: dict(N=fr.N_support, M=fr.M_support)}
    for pos_name in (loaded_name, support_name):
        e_ = ends[pos_name]; fa = dict(N=e_["N"], M=e_["M"], sig_pos=e_["N"] / h - 6 * e_["M"] / h**2, sig_neg=e_["N"] / h + 6 * e_["M"] / h**2)
        for ypos, My in (("contact (support)", rb["My_support"]), ("mid-span", rb["My_mid"])):
            My = zsign * My                                   # 물리 부호: 지지점에서 고압면이 인장 (face=+z → M_y<0)
            r = ye_recovery(g, m, dict(kyy=My / K["D22"]), n=400)   # VAM: κ_yy = M_y/D22 (M_y>0 = +z 면 압축)
            ic = int(np.argmax(r.z)) if pos_name == "crest" else int(np.argmin(r.z))
            for fname, fidx, fsig in (("+z", 0, None), ("-z", 1, None)):
                # 프레임 응력: 반전 좌표 +n 면 ↔ 물리 (+z if zsign>0 else −z)
                phys_pos_is_plus = (zsign > 0)
                sig_frame = fa["sig_pos"] if (fname == "+z") == phys_pos_is_plus else fa["sig_neg"]
                ss = sig_frame + r.sig_s[fidx, ic]
                sy = m.nu * sig_frame + r.sig_y[fidx, ic]        # 프레임: y 평면변형 → σ_y = ν σ_s
                tau = r.tau[fidx, ic]
                vm = math.sqrt(ss**2 - ss * sy + sy**2 + 3 * tau**2)
                rows.append(dict(x_pos=pos_name, y_pos=ypos, face=fname, sig_s_frame=sig_frame, sig_s_ridge=float(r.sig_s[fidx, ic]),
                                 sig_y_ridge=float(r.sig_y[fidx, ic]), sig_s=ss, sig_y=sy, tau=float(tau), vm=vm, My=My))
    best = max(rows, key=lambda d: d["vm"])
    # 접촉: 산 외면 반경 (접촉은 항상 산-산)
    Rc_out = g.extra.get("R_c", g.extra.get("R", float("nan"))) + g.h / 2
    contact = hertz_crossed_cylinders(dP * lat["A_cell"], Rc_out, lat["crossing_deg"], m.E, m.nu) if Rc_out == Rc_out else {}
    # 균질화 D11 보 근사와의 비교 (참고): 연속보 스팬 p → M(중앙)=q p²/24, M(지지)=−q p²/12 을 단위 폭 M_x 로 보고 VAM 복원 max vM
    hom = {}
    for nm, Mx in (("mid", dP * g.p**2 / 24), ("support", -dP * g.p**2 / 12)):
        rr = ye_recovery(g, m, dict(kxx=zsign * Mx / K["D11"]), n=400)
        hom[f"vm_max_{nm}"] = rr.max_vm()[0]
    hom["frame_sig_max"] = max(abs(fr.sig_pos).max(), abs(fr.sig_neg).max(), abs(fr.N_support / h) + 6 * abs(fr.M_support) / h**2)
    return PressureCase(dP, face, inplane, beta_deg, lat, fr, rb, contact, rows, best["vm"],
                        f"{best['x_pos']}, {best['y_pos']}, {best['face']} face",
                        fr.M_loaded, fr.M_support, fr.N_loaded, fr.N_support, hom)


# ============================================================================= 자체 검증
def _selfcheck() -> None:
    ok = True

    def chk(name, got, ref, tol):
        nonlocal ok
        rel = abs(got - ref) / max(abs(ref), 1e-30); ok &= rel < tol
        print(f"  [{'OK ' if rel < tol else 'BAD'}] {name}: {got:.6g} vs {ref:.6g} (rel {rel:.1e})")

    m = Sheet(193000.0, 0.3)
    print("[A] 평판 극한: 연속보 M(중앙)=q p²/24, M(지지)=−q p²/12, N≈0")
    gf = sinusoidal(10.0, 1e-6, 1.0); q = 0.5
    for ip in ("free", "fixed"):
        fr = transverse_frame(gf, m, q, "+z", ip)
        chk(f"M crest ({ip})", fr.M[0], q * 10**2 / 24 - q * fr.x[0]**2 / 2, 1e-4)
        chk(f"M near valley ({ip})", fr.M[-1], q * 10**2 / 24 - q * fr.x[-1]**2 / 2, 1e-4)      # 표본은 셀 중점 → 해석해를 같은 x 에서 비교
        assert abs(fr.N).max() < 1e-3 * q * 10 and abs(fr.P) < 1e-3 * q * 10, (abs(fr.N).max(), fr.P)   # H→0: 축력·추력은 H 에 비례해 사라짐
        chk(f"support reaction ({ip})", fr.R_support, q * 5.0, 1e-12)
    print("[B] 기본값 판: 평형·적합 검사")
    g = arc_tangent(9.0, 3.2, 0.6, R_c=1.5); q = 1.0
    for ip in ("free", "fixed"):
        fr = transverse_frame(g, m, q, "+z", ip)
        EI = m.E * g.h**3 / (12 * (1 - m.nu**2))
        rot = float(np.sum(fr.M / EI * fr.ds))
        assert abs(rot) < 1e-9 * abs(fr.M).max() / EI * fr.s[-1], rot
        # 골 지지선 수직력: 절단면 x=c 에서 Q_z = q c
        Qz_end = -(np.sum(-q * fr.ds * np.cos(fr.theta)))
        chk(f"vertical shear at support ({ip})", Qz_end, q * g.c, 1e-6)
        print(f"     {ip}: P={fr.P:.4f} N/mm  M_c={fr.M[0]:.4f}  M_v={fr.M[-1]:.4f} N·mm/mm  |σ|max={max(abs(fr.sig_pos).max(), abs(fr.sig_neg).max()):.1f} MPa")
    print("[C] 고압면 반전 대칭: R_c=R_v 이면 +z 와 −z 결과가 산↔골 거울")
    pa, pb = pressure_case(g, m, 1.0, 60.0, "+z"), pressure_case(g, m, 1.0, 60.0, "-z")
    chk("max vM +z vs −z", pb.max_vm, pa.max_vm, 1e-9)
    chk("frame M loaded ridge", pb.frame_M_loaded, pa.frame_M_loaded, 1e-9)
    print("[D] 격자·접촉")
    lat = contact_lattice(9.0, 60.0); chk("a = p/sin120°", lat["a"], 9.0 / math.sin(math.radians(120)), 1e-12)
    hz = hertz_crossed_cylinders(100.0, 1.8, 90.0, 193000.0, 0.3)
    # φ=90°: A=B=1/2R → R_e = R, 구-평면 Hertz: a=(3FR/4E*)^{1/3}
    Es = 193000 / (2 * 0.91); a_ref = (3 * 100 * 1.8 / (4 * Es)) ** (1 / 3)
    chk("Hertz 90° = sphere-plane", hz["c_mean_semi_axis"], a_ref, 1e-12)
    print("ALL OK" if ok else "SOME CHECKS FAILED")
    print("\n=== 기본값 판, ΔP = 1 MPa, β = 60°, 고압면 +z, 면내 자유 ===")
    pc = pa
    print(f"  격자: 접촉 간격 a = {pc.lattice['a']:.2f} mm, 셀 면적 {pc.lattice['A_cell']:.1f} mm², 접촉력 F = {pc.contact['F']:.1f} N")
    print(f"  (A) 프레임: M(하중측 능선) = {pc.frame_M_loaded:.3f}, M(지지선) = {pc.frame_M_support:.3f} N·mm/mm; N = {pc.frame_N_loaded:.3f}/{pc.frame_N_support:.3f} N/mm; P = {pc.frame.P:.3f}")
    print(f"      균질화 D11 보 근사 대비: frame |σ|max {pc.homogenized_Mx_check['frame_sig_max']:.1f} vs D11-beam VAM {pc.homogenized_Mx_check['vm_max_mid']:.1f}/{pc.homogenized_Mx_check['vm_max_support']:.1f} MPa")
    print(f"  (B) 능선보 (스팬 a): |M_y| 지지 = {abs(pc.ridge['My_support']):.3f}, 중앙 = {abs(pc.ridge['My_mid']):.3f} N·mm/mm")
    print(f"  (C) Hertz: p0 = {pc.contact['p0']:.0f} MPa (반타원 평균 반축 {pc.contact['c_mean_semi_axis']:.3f} mm)")
    print(f"  최대 von Mises = {pc.max_vm:.1f} MPa at {pc.max_where}")
    for r in pc.table:
        print(f"     {r['x_pos']:6s} {r['y_pos']:18s} {r['face']:2s}: σ_s={r['sig_s']:8.1f} (frame {r['sig_s_frame']:7.1f} + ridge {r['sig_s_ridge']:6.1f})  σ_y={r['sig_y']:7.1f}  vM={r['vm']:7.1f}")


if __name__ == "__main__":
    _selfcheck()
