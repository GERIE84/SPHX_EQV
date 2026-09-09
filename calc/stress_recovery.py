"""등가 판 응력합력 → 주름 원판 국부 응력 1차 복원 (B-4 예비, Briassoulis 1986 Appendix B 기반).

등가 직교이방성 판 해석으로 얻은 단위폭당 합력 N_x*, M_y* 로부터, 능선(z=±f)에서의 원판 섬유 응력을 추정한다.
국부축: x ⊥ 능선(파형 진행), y ∥ 능선.  기호는 docs/05, docs/06 을 따른다 (f = H/2, t = 판 두께).

[A] 파형 가로 인장 N_x*  (Briassoulis Eq. 6, B1–B3; Xia 2012 Case 1 과 동일 물리)
    국부 모멘트 M_x(z) = N_x* z  →  σ_x(z) = N_x*/t · (1 ∓ 6z/t)   (∓: 외측/내측 섬유)
    능선 내측 최대: σ_x,in = N_x*/t (1 + 6f/t)  → K_t = 1 + 6f/t
    능선 외측:      σ_x,out = N_x*/t (1 − 6f/t)  (f/t > 1/3 이면 부호 반전)
[B] 능선 방향 굽힘 M_y*  (Briassoulis A5, A6, B4–B7)
    균일 모멘트 M = M_y* / (1 + A f²/(2I)) ,  A = t (단위폭 단면적), I = t³/12
    축력     N_y(z) = M_y* · z / (I/A − f²/2)      ( = M_y* f sin(πx/c) / (I/A − f²/2) 을 z = f sin(πx/c) 로 쓴 것)
    σ_y(z) = ∓ N_y/t ∓ 6M/t²  →  Eq. (B5):
        σ_y = (6 M_y*/t²) · [ ∓ (z/t)·2(1−μ²) ∓ 1 ] / [ (f/t)² 6(1−μ²) + 1 ]
    (원문은 사인형 가정; 임의 단면에도 z = 국부 높이로 두면 1차 근사로 사용 가능)

주의: Kirchhoff 박판 가정. t/R_min 이 0.2 를 넘으면 원호부 두께 방향 응력 분포가 비선형 → 상세 FE 로 보정 (docs/05 §4).
"""
from __future__ import annotations

from dataclasses import dataclass


@dataclass
class RidgeStress:
    inner: float   # 오목(내측) 섬유 응력
    outer: float   # 볼록(외측) 섬유 응력
    nominal: float # 등가 판 공칭응력


def kt_x_tension(f: float, t: float) -> float:
    """능선 내측 섬유 응력집중계수 (Briassoulis Eq. B2): K_t = 1 + 6 f/t."""
    return 1.0 + 6.0 * f / t


def sigma_x_ridge(N_x: float, f: float, t: float) -> RidgeStress:
    """파형 가로 인장 N_x* [N/mm] 에 의한 능선(z=f) 섬유 응력 [MPa] (Briassoulis B2, B3)."""
    s0 = N_x / t
    return RidgeStress(inner=s0 * (1 + 6 * f / t), outer=s0 * (1 - 6 * f / t), nominal=s0)


def sigma_y_curvature(M_y: float, f: float, t: float, mu: float, z: float) -> tuple[float, float]:
    """능선 방향 굽힘 M_y* [N·mm/mm] 상태에서 높이 z 인 지점의 (외측, 내측) 섬유 응력 (Briassoulis B5).
    z = +f (상부 능선), z = -f (하부 능선)."""
    k = 6 * M_y / t**2 / ((f / t) ** 2 * 6 * (1 - mu**2) + 1)
    a = (z / t) * 2 * (1 - mu**2)
    return (k * (-a - 1), k * (-a + 1))          # (outer, inner) — 부호는 M_y* 양의 방향 기준


def uniform_moment_and_axial(M_y: float, f: float, t: float) -> tuple[float, float]:
    """Briassoulis A5/A6: 균일 모멘트 M 와 능선 축력 N_y(z=f)."""
    A, I = t, t**3 / 12
    M = M_y / (1 + A * f**2 / (2 * I))
    N_ridge = M_y * f / (I / A - f**2 / 2) if abs(I / A - f**2 / 2) > 1e-30 else float("inf")
    return M, N_ridge


def _selfcheck() -> None:
    # Briassoulis Fig. 6: σ_x,max/σ_x,orth 는 f/t 에 선형, f/t=8 에서 약 49~50
    assert abs(kt_x_tension(8.0, 1.0) - 49.0) < 1e-12
    # Briassoulis §4.1 판: f=0.21875 in, t=0.25 in → K_t = 6.25 ; 외측은 1-6f/t = -4.25 (부호 반전, f/t>1/3)
    r = sigma_x_ridge(1.0, 0.21875, 0.25)
    assert abs(r.inner / r.nominal - 6.25) < 1e-12 and abs(r.outer / r.nominal + 4.25) < 1e-12
    # 평판 극한 f→0: K_t → 1, σ_y → ±6M/t² (순수 판 굽힘)
    assert abs(kt_x_tension(0.0, 1.0) - 1.0) < 1e-12
    o, i = sigma_y_curvature(1.0, 0.0, 1.0, 0.3, 0.0)
    assert abs(o + 6.0) < 1e-12 and abs(i - 6.0) < 1e-12
    print("stress_recovery selfcheck OK")


if __name__ == "__main__":
    _selfcheck()
    print("\nPHE 가정 기본값 (f=1.6, t=0.6 mm): K_t(x-tension, ridge inner) =", round(kt_x_tension(1.6, 0.6), 2))
    r = sigma_x_ridge(N_x=10.0, f=1.6, t=0.6)
    print(f"  N_x* = 10 N/mm → 공칭 {r.nominal:.1f} MPa, 능선 내측 {r.inner:.1f} MPa, 외측 {r.outer:.1f} MPa")
    M, Nr = uniform_moment_and_axial(M_y=100.0, f=1.6, t=0.6)
    o, i = sigma_y_curvature(100.0, 1.6, 0.6, 0.3, 1.6)
    print(f"  M_y* = 100 N·mm/mm → 균일 M = {M:.3f} N·mm/mm, 능선 축력 N_y = {Nr:.1f} N/mm, 상부 능선 섬유응력 (외측 {o:.1f}, 내측 {i:.1f}) MPa")
    print(f"  비교: 등가 판 공칭 굽힘응력 6M_y*/t_e² (t_e=t) = {6*100/0.36:.0f} MPa → 실제 원판 응력은 훨씬 작음 (하중을 단면 2차 모멘트가 부담)")
