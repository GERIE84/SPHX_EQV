"""국부축(x ⊥ 능선, y ∥ 능선) 등가 강성 → 판축(x_p, y_p) 회전 변환 및 쉐브론 영역 조합 — B-3 구현.

docs/07_stiffness_transformation.md 의 수식을 구현한다.

각도 규약 (docs/05 §1, 그림 (c)):
    β  = 능선(y)과 유동축(y_p) 사이 각.  apex 우측 영역 +β, 좌측 −β.
    θ  = 판축 x_p 에서 국부축 x 까지의 회전각 (반시계 +).  우측 영역 θ = −β, 좌측 영역 θ = +β.
    (국부 x = (cos θ, sin θ) 를 판축 성분으로 쓴 것.  우측: x = (cos β, −sin β) → θ = −β.)

행렬 규약: 3×3, 성분 순서 (xx, yy, xy),  N = A {ε_x, ε_y, γ_xy},  M = D {κ_x, κ_y, 2κ_xy}  (공학 전단·비틀림).
A 와 D 는 같은 4차 텐서 변환을 따른다 (고전 적층판 이론의 Q̄ 변환과 동일).
"""
from __future__ import annotations

import math
from typing import Dict, Tuple

import numpy as np

from geometry import CorrugationGeometry, arc_tangent
from stiffness import Sheet, adopted

IDX = ("11", "22", "12", "66", "16", "26")


# ----------------------------------------------------------------------------- dict <-> 3x3
def to_matrix(K: Dict[str, float], prefix: str) -> np.ndarray:
    """{'A11','A12','A22','A66'[,'A16','A26']} → 3×3 (xx, yy, xy)."""
    g = lambda k: K.get(prefix + k, 0.0)
    return np.array([[g("11"), g("12"), g("16")],
                     [g("12"), g("22"), g("26")],
                     [g("16"), g("26"), g("66")]])


def to_dict(M: np.ndarray, prefix: str) -> Dict[str, float]:
    return {prefix + "11": M[0, 0], prefix + "22": M[1, 1], prefix + "12": M[0, 1],
            prefix + "66": M[2, 2], prefix + "16": M[0, 2], prefix + "26": M[1, 2]}


# ----------------------------------------------------------------------------- 회전 변환
def strain_transform(theta: float) -> np.ndarray:
    """T_ε: 판축 공학변형률 {ε_xp, ε_yp, γ} → 국부축 {ε_x, ε_y, γ_xy}.  θ = x_p → x 회전각."""
    c, s = math.cos(theta), math.sin(theta)
    return np.array([[c * c, s * s, s * c],
                     [s * s, c * c, -s * c],
                     [-2 * s * c, 2 * s * c, c * c - s * s]])


def stress_transform(theta: float) -> np.ndarray:
    """T_σ: 판축 합력 {N_xp, N_yp, N_xpyp} → 국부축 {N_x, N_y, N_xy}."""
    c, s = math.cos(theta), math.sin(theta)
    return np.array([[c * c, s * s, 2 * s * c],
                     [s * s, c * c, -2 * s * c],
                     [-s * c, s * c, c * c - s * s]])


def rotate_stiffness(K_local: np.ndarray, theta: float) -> np.ndarray:
    """국부축 강성 K (3×3) → 판축 강성 K̄ = T_σ⁻¹ K T_ε.  (N_p = T_σ⁻¹ N_local, ε_local = T_ε ε_p)"""
    return np.linalg.inv(stress_transform(theta)) @ K_local @ strain_transform(theta)


def rotate_stiffness_closed(K: Dict[str, float], prefix: str, theta: float) -> Dict[str, float]:
    """고전 적층판 이론 폐형식 (Jones 1975 Eq. 2.80 형태). 검증용 — rotate_stiffness 와 일치해야 함."""
    q11, q12, q22, q66 = (K[prefix + "11"], K[prefix + "12"], K[prefix + "22"], K[prefix + "66"])
    m, n = math.cos(theta), math.sin(theta)
    m2, n2, m4, n4 = m * m, n * n, m**4, n**4
    return {
        prefix + "11": q11 * m4 + 2 * (q12 + 2 * q66) * n2 * m2 + q22 * n4,
        prefix + "22": q11 * n4 + 2 * (q12 + 2 * q66) * n2 * m2 + q22 * m4,
        prefix + "12": (q11 + q22 - 4 * q66) * n2 * m2 + q12 * (n4 + m4),
        prefix + "66": (q11 + q22 - 2 * q12 - 2 * q66) * n2 * m2 + q66 * (n4 + m4),
        prefix + "16": (q11 - q12 - 2 * q66) * n * m**3 + (q12 - q22 + 2 * q66) * n**3 * m,
        prefix + "26": (q11 - q12 - 2 * q66) * n**3 * m + (q12 - q22 + 2 * q66) * n * m**3,
    }


def zone_theta(beta_deg: float, zone: str) -> float:
    """apex 우측(+β 영역) θ = −β, 좌측(−β 영역) θ = +β  [rad]."""
    b = math.radians(beta_deg)
    return -b if zone == "right" else b


def zone_stiffness(K_local: Dict[str, float], beta_deg: float, zone: str) -> Dict[str, float]:
    """국부축 채택식 dict (A.., D..) → 해당 영역의 판축 강성 dict (A16, A26, D16, D26 포함)."""
    th = zone_theta(beta_deg, zone)
    out = {}
    for pre in ("A", "D"):
        out.update(to_dict(rotate_stiffness(to_matrix(K_local, pre), th), pre))
    return out


def resultants_to_local(N_plate: np.ndarray, beta_deg: float, zone: str) -> np.ndarray:
    """판축 합력 {N_xp, N_yp, N_xpyp} (또는 M) → 국부축 {N_x, N_y, N_xy}.  국부 응력 복원(stress_recovery) 입력용."""
    return stress_transform(zone_theta(beta_deg, zone)) @ np.asarray(N_plate, float)


# ----------------------------------------------------------------------------- 영역 조합 (±β 스트립이 x_p 방향으로 교대)
def voigt(K1: np.ndarray, K2: np.ndarray, v1: float = 0.5) -> np.ndarray:
    """균일 변형률 가정 (상한). ±β 대칭 조합이면 16, 26 항이 정확히 상쇄."""
    return v1 * K1 + (1 - v1) * K2


def reuss(K1: np.ndarray, K2: np.ndarray, v1: float = 0.5) -> np.ndarray:
    """균일 합력 가정 (하한)."""
    return np.linalg.inv(v1 * np.linalg.inv(K1) + (1 - v1) * np.linalg.inv(K2))


def strip_laminate(K1: np.ndarray, K2: np.ndarray, v1: float, cont_strain: Tuple[int, ...]) -> np.ndarray:
    """경계가 y_p 에 평행한 두 스트립의 정확한 균질화 (rank-1 laminate).
    cont_strain: 경계를 가로질러 연속인 변형률 성분 인덱스 (막: (1,)=ε_yp ; 굽힘: (1,2)=κ_yp, 2κ_xpyp).
    나머지 성분은 합력(견인력)이 연속: 막 N_xp, N_xpyp ; 굽힘 M_xp."""
    b = list(cont_strain)
    a = [i for i in range(3) if i not in b]
    Ks, vs = (K1, K2), (v1, 1 - v1)
    Sbar = sum(v * np.linalg.inv(K[np.ix_(a, a)]) for K, v in zip(Ks, vs))                     # Σ v K_aa⁻¹
    P = sum(v * np.linalg.inv(K[np.ix_(a, a)]) @ K[np.ix_(a, b)] for K, v in zip(Ks, vs))     # Σ v K_aa⁻¹ K_ab
    Q = sum(v * K[np.ix_(b, a)] @ np.linalg.inv(K[np.ix_(a, a)]) for K, v in zip(Ks, vs))     # Σ v K_ba K_aa⁻¹
    R = sum(v * (K[np.ix_(b, b)] - K[np.ix_(b, a)] @ np.linalg.inv(K[np.ix_(a, a)]) @ K[np.ix_(a, b)]) for K, v in zip(Ks, vs))
    Sinv = np.linalg.inv(Sbar)
    Keff = np.zeros((3, 3))
    Keff[np.ix_(a, a)] = Sinv                      # N_a = S̄⁻¹ ε̄_a + S̄⁻¹ P ε_b
    Keff[np.ix_(a, b)] = Sinv @ P
    Keff[np.ix_(b, a)] = Q @ Sinv                  # N̄_b = Q N_a + R ε_b
    Keff[np.ix_(b, b)] = Q @ Sinv @ P + R
    return 0.5 * (Keff + Keff.T)                   # 대칭화 (수치 오차 제거; 이론적으로 대칭)


def chevron_homogenized(K_local: Dict[str, float], beta_deg: float, v_right: float = 0.5) -> Dict[str, Dict[str, float]]:
    """±β 영역 조합의 세 추정치 (voigt / reuss / strip) 를 A, D 각각 반환."""
    out = {}
    for pre, cont in (("A", (1,)), ("D", (1, 2))):
        KR = rotate_stiffness(to_matrix(K_local, pre), zone_theta(beta_deg, "right"))
        KL = rotate_stiffness(to_matrix(K_local, pre), zone_theta(beta_deg, "left"))
        out.setdefault("voigt", {}).update(to_dict(voigt(KR, KL, v_right), pre))
        out.setdefault("reuss", {}).update(to_dict(reuss(KR, KL, v_right), pre))
        out.setdefault("strip", {}).update(to_dict(strip_laminate(KR, KL, v_right, cont), pre))
    return out


# ----------------------------------------------------------------------------- 자체 검증
def _selfcheck() -> None:
    ok = True

    def chk(name, cond):
        nonlocal ok
        ok &= bool(cond)
        print(f"  [{'OK ' if cond else 'BAD'}] {name}")

    g = arc_tangent(9.0, 3.2, 0.6, R_c=1.5); m = Sheet(193000.0, 0.3)
    K = adopted(g, m)
    A, D = to_matrix(K, "A"), to_matrix(K, "D")
    # 1) 폐형식 == 행렬식
    for th in (0.3, -1.1, 2.0):
        Ac = to_matrix(rotate_stiffness_closed(K, "A", th), "A")
        chk(f"closed-form == matrix rotation (θ={th})", np.allclose(Ac, rotate_stiffness(A, th), rtol=1e-12, atol=1e-9))
    # 2) θ=0 identity, θ=90° swap
    chk("θ=0 → identity", np.allclose(rotate_stiffness(A, 0.0), A))
    R90 = rotate_stiffness(A, math.pi / 2)
    chk("θ=90° → 11↔22 swap, 66 unchanged, 16=26=0", abs(R90[0, 0] - A[1, 1]) < 1e-6 and abs(R90[1, 1] - A[0, 0]) < 1e-6
        and abs(R90[2, 2] - A[2, 2]) < 1e-6 and abs(R90[0, 2]) < 1e-6 and abs(R90[1, 2]) < 1e-6)
    # 3) 불변량: U1 = A11+A22+2A12 (trace-like), U2 = A66 - A12  는 회전 불변
    for th in (0.4, 1.3):
        R = rotate_stiffness(A, th)
        chk(f"invariants A11+A22+2A12, A66−A12 (θ={th})",
            abs((R[0, 0] + R[1, 1] + 2 * R[0, 1]) - (A[0, 0] + A[1, 1] + 2 * A[0, 1])) < 1e-6 and abs((R[2, 2] - R[0, 1]) - (A[2, 2] - A[0, 1])) < 1e-6)
    # 4) 대칭성 유지, ±β 에서 16/26 부호 반전
    KR, KL = rotate_stiffness(A, zone_theta(60, "right")), rotate_stiffness(A, zone_theta(60, "left"))
    chk("rotated matrix symmetric", np.allclose(KR, KR.T))
    chk("±β zones: A16, A26 opposite sign, others equal", abs(KR[0, 2] + KL[0, 2]) < 1e-6 and abs(KR[1, 2] + KL[1, 2]) < 1e-6
        and np.allclose(KR[:2, :2], KL[:2, :2]) and abs(KR[2, 2] - KL[2, 2]) < 1e-6)
    # 5) 에너지 보존: ε_pᵀ K̄ ε_p == ε_localᵀ K ε_local
    th = zone_theta(60, "right"); ep = np.array([1e-3, -4e-4, 7e-4])
    el = strain_transform(th) @ ep
    chk("strain energy invariant under rotation", abs(ep @ rotate_stiffness(A, th) @ ep - el @ A @ el) < 1e-9 * abs(el @ A @ el))
    # 6) 합력 변환 일관성: N_local = T_σ N_p  with N_p = K̄ ε_p, N_local = K ε_local
    chk("resultant transform consistent", np.allclose(stress_transform(th) @ (rotate_stiffness(A, th) @ ep), A @ el))
    # 7) 조합 규칙: β=0 → voigt = reuss = strip = K ; 등방 평판 → 회전 불변
    H = chevron_homogenized(K, 0.0)
    chk("β=0: voigt=reuss=strip=K", all(abs(H[k]["A11"] - K["A11"]) < 1e-6 and abs(H[k]["D22"] - K["D22"]) < 1e-6 for k in H))
    iso = dict(A11=1.0, A12=0.3, A22=1.0, A66=0.35, D11=1.0, D12=0.3, D22=1.0, D66=0.35)
    chk("isotropic K rotation-invariant", np.allclose(rotate_stiffness(to_matrix(iso, "A"), 0.7), to_matrix(iso, "A")))
    # 8) 스트립 정확해는 Voigt/Reuss 사이 (에너지 부등식)
    HA = chevron_homogenized(K, 60.0)
    for comp in ("A11", "A22", "A66", "D11", "D22", "D66"):
        chk(f"reuss ≤ strip ≤ voigt for {comp}", HA["reuss"][comp] - 1e-6 <= HA["strip"][comp] <= HA["voigt"][comp] + 1e-6)
    print("ALL OK" if ok else "SOME CHECKS FAILED")


def _table(K: Dict[str, float], beta_deg: float) -> str:
    Z = zone_stiffness(K, beta_deg, "right")
    Hm = chevron_homogenized(K, beta_deg)
    rows = []
    for pre in ("A", "D"):
        for k in IDX:
            key = pre + k
            rows.append(f"  {key:4} {K.get(key, 0.0):12.4g} {Z[key]:12.4g} {Hm['voigt'][key]:12.4g} {Hm['strip'][key]:12.4g} {Hm['reuss'][key]:12.4g}")
    head = f"  {'':4} {'local':>12} {'right(+β)':>12} {'voigt':>12} {'strip':>12} {'reuss':>12}"
    return head + "\n" + "\n".join(rows)


if __name__ == "__main__":
    _selfcheck()
    g = arc_tangent(9.0, 3.2, 0.6, R_c=1.5); m = Sheet(193000.0, 0.3)
    K = adopted(g, m)
    for beta in (30.0, 45.0, 60.0):
        print(f"\n=== β = {beta:.0f}°  (가정 기본값 판, N·mm 단위; 좌측 영역은 16/26 부호 반전) ===")
        print(_table(K, beta))
