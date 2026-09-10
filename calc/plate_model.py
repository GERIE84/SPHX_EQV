"""C-1: YAML 입력 → 등가 판 계산 모듈 (geometry → stiffness → transform → equivalent_plate → stress_recovery) → 보고서.

사용:
    python calc/plate_model.py calc/plate_input_template.yaml -o out/example
    python calc/plate_model.py --selfcheck
출력 (out_dir):
    report.md      사람이 읽는 계산 보고서 (입력 전처리·형상·강성·등가 두께·영역 강성·단위 응답·하중 케이스·경고)
    results.json   같은 내용의 기계 판독 값 (엑셀 시트 C-2, 파라메트릭 C-3 의 기준값)
    sections.inp   ANSYS preintegrated general shell section APDL 블록 (docs/09 §2)

입력 전처리 규칙 (docs/05 §2 표기 기준 → 계산 기준):
    dims_reference=outer_surface  → H = H_input − t  (외면 산-골 높이 H_o = H + t)
    pitch_reference=along_plate_axis → p = p_input · cos β  (x_p 방향 측정 피치를 능선 수직 피치로 환산)
    beta_reference=horizontal      → β = 90° − β_input
    계산 두께 t_calc = (t_min 또는 t) − corrosion_allowance ; 강성·응력 모두 t_calc 사용 (05 §2.5)
    R_c, R_v 는 중앙면 값으로 입력된 것으로 본다 (outer 기준이면 사용자가 R_o − t/2 로 환산해 기입)
하중 케이스 (선택, `load_cases:` 목록): axes=local|plate, zone=right|left, N=[Nx,Ny,Nxy] N/mm, M=[Mx,My,Mxy] N·mm/mm,
    allowable(MPa, 선택; 없으면 material.S_allow). plate 축 합력은 transform.resultants_to_local 로 영역 국부축으로 변환.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
from dataclasses import dataclass, field, asdict
from typing import Any, Dict, List, Optional

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import geometry as geo                                   # noqa: E402
from geometry import CorrugationGeometry                 # noqa: E402
from stiffness import Sheet, MODELS, adopted             # noqa: E402
from transform import zone_stiffness, chevron_homogenized, resultants_to_local, IDX  # noqa: E402
from equivalent_plate import (equivalent_thicknesses, homogeneous_orthotropic, transverse_shear_estimate,
                              mass_per_area, chevron_gens_deck)  # noqa: E402
from stress_recovery import ye_recovery, recovery_from_resultants, kt_x_tension, RecoveryResult  # noqa: E402

try:
    import yaml
except ImportError:  # pragma: no cover
    yaml = None

KEYS_A = [f"A{k}" for k in ("11", "12", "22", "66", "16", "26")]
KEYS_D = [f"D{k}" for k in ("11", "12", "22", "66", "16", "26")]


# ============================================================================= 입력 전처리
@dataclass
class PreparedInput:
    raw: Dict[str, Any]
    plate_id: str
    profile: str
    p: float
    H: float
    t_nominal: float
    t_calc: float
    section_kw: Dict[str, Any]
    beta_deg: float
    E: float
    nu: float
    rho: float
    S_allow: Optional[float]
    notes: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    load_cases: List[Dict[str, Any]] = field(default_factory=list)


def _g(d: Optional[dict], k: str, default=None):
    v = (d or {}).get(k)
    return default if v is None else v


def prepare(raw: Dict[str, Any]) -> PreparedInput:
    sec, chev, mat = raw.get("section", {}), raw.get("chevron", {}), raw.get("material", {})
    notes, warns = [], []
    beta = float(chev["beta_deg"])
    if _g(chev, "beta_reference", "flow_axis") == "horizontal":
        beta = 90.0 - beta; notes.append(f"β 는 수평 기준으로 입력됨 → 유동축 기준 β = {beta:g}°")
    t_nom = float(sec["t"])
    t_min = _g(mat, "t_min", t_nom); ca = float(_g(mat, "corrosion_allowance", 0.0))
    t_calc = float(t_min) - ca
    if t_min != t_nom or ca:
        notes.append(f"계산 두께 t_calc = t_min({t_min:g}) − CA({ca:g}) = {t_calc:g} mm (공칭 {t_nom:g})")
    if _g(mat, "t_min") is None:
        warns.append("t_min 미입력 → 공칭 두께로 계산. 성형 후 두께 감소(능선부 5~15% 통상)를 반영하지 않음")
    H = float(sec["H"])
    if _g(sec, "dims_reference", "mid_surface") == "outer_surface":
        H = H - t_nom; notes.append(f"H 는 외면 기준 입력 → 중앙면 H = {H:g} mm")
    p = float(sec["p"])
    if _g(sec, "pitch_reference", "normal_to_ridge") == "along_plate_axis":
        p = p * math.cos(math.radians(beta)); notes.append(f"피치는 x_p 방향 측정 → 능선 수직 피치 p = p_in·cosβ = {p:.4g} mm")
    kw = dict(profile=sec["profile"], p=p, H=H, t=t_calc)
    for k in ("R_c", "R_v", "alpha_deg", "R", "L"):
        if _g(sec, k) is not None:
            kw[k] = float(sec[k])
    rho = float(_g(mat, "rho", 7.9e-9))
    if _g(mat, "rho") is None:
        notes.append("밀도 미입력 → 7.9e-9 tonne/mm³ (스테인리스) 가정")
    if _g(mat, "T_design") is not None and _g(mat, "E_source") is None:
        warns.append(f"설계온도 {mat['T_design']} °C 입력됨 — E={mat['E']:g} MPa 가 해당 온도 값인지 확인")
    return PreparedInput(raw=raw, plate_id=str(_g(raw.get("meta"), "plate_id", "plate")), profile=sec["profile"], p=p, H=H,
                         t_nominal=t_nom, t_calc=t_calc, section_kw=kw, beta_deg=beta, E=float(mat["E"]), nu=float(mat["nu"]),
                         rho=rho, S_allow=_g(mat, "S_allow"), notes=notes, warnings=warns,
                         load_cases=list(raw.get("load_cases") or []))


# ============================================================================= 계산
@dataclass
class LoadCaseResult:
    name: str
    axes: str
    zone: str
    N_plate: List[float]; M_plate: List[float]
    N_local: List[float]; M_local: List[float]
    macro: Dict[str, float]
    max_vm: float; at_x: float; at_z: float; fibre: str
    ridge_sig_s: List[float]; ridge_sig_y: List[float]; ridge_tau: List[float]
    ridge_N_s: float; ridge_N_y: float; ridge_M_s: float; ridge_M_y: float
    membrane_vm_max: float
    allowable: Optional[float]; ratio: Optional[float]; verdict: Optional[str]


def _z(v: float, tol: float = 1e-9) -> float:
    """수치 잡음(±1e-17 등) 을 0 으로 정리."""
    return 0.0 if abs(v) < tol else float(v)


def _recovery_summary(r: RecoveryResult) -> Dict[str, Any]:
    vm, side, i = r.max_vm(); ic = int(np.argmax(r.z))
    Q_vm_m = np.sqrt(r.N_s**2 - r.N_s * r.N_y + r.N_y**2)  # 막 합력 등가 [N/mm]
    return dict(max_vm=vm, at_x=float(r.x[i]), at_z=float(r.z[i]), fibre="+z" if side == 0 else "-z", ic=ic,
                ridge_sig_s=[_z(r.sig_s[0, ic]), _z(r.sig_s[1, ic])], ridge_sig_y=[_z(r.sig_y[0, ic]), _z(r.sig_y[1, ic])],
                ridge_tau=[_z(r.tau[0, ic]), _z(r.tau[1, ic])], ridge_N_s=_z(r.N_s[ic]), ridge_N_y=_z(r.N_y[ic]),
                ridge_M_s=_z(r.M_s[ic]), ridge_M_y=_z(r.M_y[ic]), membrane_eq_max=float(np.max(Q_vm_m)))


class PlateModel:
    def __init__(self, inp: PreparedInput):
        self.inp = inp
        self.g: CorrugationGeometry = geo.from_dict(inp.section_kw)
        self.m = Sheet(inp.E, inp.nu)
        self.K = adopted(self.g, self.m)
        self.models = {name: fn(self.g, self.m) for name, fn in MODELS.items()}
        self.th = equivalent_thicknesses(self.K, self.g)
        self.homog = {basis: homogeneous_orthotropic(self.K, self.th["t_b"], basis) for basis in ("A", "D")}
        self.zones = {z: zone_stiffness(self.K, inp.beta_deg, z) for z in ("right", "left")}
        self.homogenized = chevron_homogenized(self.K, inp.beta_deg)
        self.E_ts = transverse_shear_estimate(self.g, self.m)
        self.dens_area = mass_per_area(inp.rho, self.g)
        self.kt = kt_x_tension(self.g.f, self.g.h)
        self.unit_response = self._unit_response()
        self.load_results = [self._load_case(lc) for lc in inp.load_cases]
        self.warnings = list(inp.warnings) + self._range_checks()

    # --- 검사
    def _range_checks(self) -> List[str]:
        w = []; d = self.g.dimensionless()
        if d.get("t_over_Rmin", 0) > 0.2:
            w.append(f"t/R_min = {d['t_over_Rmin']:.2f} > 0.2: 쉘 이론 유효범위 밖. 강성은 유지, 능선 응력은 D단계 솔리드 FE 보정계수 필요 (09 §5)")
        if d["H_over_p"] > 1.0:
            w.append(f"H/p = {d['H_over_p']:.2f} > 1: 깊은 주름. 채택식(VAM)은 유효하나 문헌 검증 범위(≤0.8) 밖")
        if self.th["t_b"] > 0.3 * self.g.p:
            w.append(f"t_b/p = {self.th['t_b'] / self.g.p:.2f}: 등가 판 두께가 피치와 같은 차수. 등가 판 FE 요소 크기 ≥ 2p, 스팬/t_b ≫ 10 을 확인하고 횡전단(SSPE) 감도를 본다")
        if self.inp.raw.get("chevron", {}).get("n_apex", 1) and self.inp.raw.get("plate", {}).get("D_e") is None and self.inp.raw.get("plate", {}).get("W_e") is None:
            w.append("판 유효 주름 영역 크기(D_e 또는 W_e) 미입력 → apex/테두리 보정, 전체 판 합력 산정은 이 보고서 범위 밖")
        return w

    # --- 단위 합력 응답 (국부축)
    def _unit_response(self) -> List[Dict[str, Any]]:
        out = []
        for name, N, M in (("N_x", [1, 0, 0], [0, 0, 0]), ("N_y", [0, 1, 0], [0, 0, 0]), ("N_xy", [0, 0, 1], [0, 0, 0]),
                           ("M_x", [0, 0, 0], [1, 0, 0]), ("M_y", [0, 0, 0], [0, 1, 0]), ("M_xy", [0, 0, 0], [0, 0, 1])):
            r = recovery_from_resultants(self.g, self.m, N, M)
            s = _recovery_summary(r)
            out.append(dict(resultant=name, unit="N/mm" if name.startswith("N") else "N·mm/mm", max_vm_per_unit=s["max_vm"],
                            at_x=s["at_x"], at_z=s["at_z"], fibre=s["fibre"], ridge_sig_s=s["ridge_sig_s"], ridge_sig_y=s["ridge_sig_y"]))
        return out

    # --- 하중 케이스
    def _load_case(self, lc: Dict[str, Any]) -> LoadCaseResult:
        axes = lc.get("axes", "local"); zone = lc.get("zone", "right")
        N = [float(v) for v in lc.get("N", [0, 0, 0])]; M = [float(v) for v in lc.get("M", [0, 0, 0])]
        if axes == "plate":
            Nl = resultants_to_local(np.array(N), self.inp.beta_deg, zone).tolist()
            Ml = resultants_to_local(np.array(M), self.inp.beta_deg, zone).tolist()
        else:
            Nl, Ml = N, M
        A = np.array([[self.K["A11"], self.K["A12"], 0], [self.K["A12"], self.K["A22"], 0], [0, 0, self.K["A66"]]])
        D = np.array([[self.K["D11"], self.K["D12"], 0], [self.K["D12"], self.K["D22"], 0], [0, 0, self.K["D66"]]])
        e, k = np.linalg.solve(A, Nl), np.linalg.solve(D, Ml)
        macro = dict(exx=e[0], eyy=e[1], gxy=e[2], kxx=k[0], kyy=k[1], kxy=k[2] / 2)
        r = ye_recovery(self.g, self.m, macro); s = _recovery_summary(r)
        # 막 응력만의 von Mises (1차 일반 막응력 판정용): 국부 합력/h
        vm_m = float(np.max(np.sqrt(r.N_s**2 - r.N_s * r.N_y + r.N_y**2 + 3 * (self.m.G * self.g.h * r.gam_sy) ** 2)) / self.g.h)
        allow = lc.get("allowable", self.inp.S_allow)
        ratio = s["max_vm"] / allow if allow else None
        verdict = None if ratio is None else ("OK" if ratio <= 1.0 else "NG")
        return LoadCaseResult(str(lc.get("name", "case")), axes, zone, N, M, [_z(v) for v in Nl], [_z(v) for v in Ml],
                              {k_: float(v) for k_, v in macro.items()}, s["max_vm"], s["at_x"], s["at_z"], s["fibre"],
                              s["ridge_sig_s"], s["ridge_sig_y"], s["ridge_tau"], s["ridge_N_s"], s["ridge_N_y"], s["ridge_M_s"], s["ridge_M_y"],
                              vm_m, allow, ratio, verdict)

    # --- 출력
    def results(self) -> Dict[str, Any]:
        i = self.inp
        return dict(
            plate_id=i.plate_id, units="mm-N-MPa-tonne-deg",
            input_prepared=dict(profile=i.profile, p=i.p, H=i.H, t_nominal=i.t_nominal, t_calc=i.t_calc, beta_deg=i.beta_deg,
                                E=i.E, nu=i.nu, G=self.m.G, rho=i.rho, S_allow=i.S_allow, section_kw=i.section_kw, notes=i.notes),
            geometry=dict(c=self.g.c, f=self.g.f, h=self.g.h, l=self.g.l, I1=self.g.I1, I2=self.g.I2, lam=self.g.l / self.g.c,
                          J1=self.g.avg_inv_sqrt_a, J2=self.g.eps2_avg_phi2_sqrt_a, A_u=self.g.A_u, I_u=self.g.I_u, W_u=self.g.W_u,
                          extra=dict(self.g.extra), dimensionless=self.g.dimensionless()),
            stiffness_local=dict(adopted={k: float(self.K[k]) for k in KEYS_A[:4] + KEYS_D[:4]},
                                 models={n: {k: float(v[k]) for k in KEYS_A[:4] + KEYS_D[:4]} for n, v in self.models.items()}),
            equivalent_thickness={k: float(v) for k, v in self.th.items()},
            homogeneous_single_layer={b: {k: (v if not isinstance(v, dict) else {kk: float(vv) for kk, vv in v.items()}) for k, v in h.items()}
                                      for b, h in self.homog.items()},
            zones={z: {k: float(v[k]) for k in KEYS_A + KEYS_D} for z, v in self.zones.items()},
            homogenized={n: {k: float(v.get(k, 0.0)) for k in KEYS_A + KEYS_D} for n, v in self.homogenized.items()},
            fe_input=dict(method="ANSYS preintegrated general shell section (SECTYPE,,GENS)", transverse_shear=self.E_ts,
                          dens_unit_thickness=self.dens_area, twist_convention="to verify (E-1)"),
            stress=dict(K_t_ridge=self.kt, unit_response=self.unit_response, load_cases=[asdict(r) for r in self.load_results]),
            warnings=self.warnings,
        )

    def apdl(self) -> str:
        return chevron_gens_deck(self.g, self.m, self.inp.rho, self.inp.beta_deg)

    def report_md(self) -> str:
        i, g, K, th = self.inp, self.g, self.K, self.th
        L: List[str] = []
        a = L.append
        a(f"# 등가 판 계산 보고서 — {i.plate_id}")
        a(f"\n생성: `calc/plate_model.py`, 입력: {i.raw.get('meta', {}).get('source', '-')}, 날짜 {i.raw.get('meta', {}).get('date', '-')}. 단위 mm–N–MPa. 정의·근거: docs/05 (형상), 06 (채택식), 07 (β 회전), 09 (등가 두께·GENS·응력 복원).\n")
        if self.warnings:
            a("## ⚠ 경고·확인 필요\n")
            for w in self.warnings: a(f"- {w}")
            a("")
        a("## 1. 입력 (전처리 후)\n")
        a("| 항목 | 값 |\n|---|---|")
        a(f"| 단면 | {i.profile}, p = {i.p:g}, H = {i.H:g} mm |")
        ex = ", ".join(f"{k} = {v:g}" for k, v in i.section_kw.items() if k not in ("profile", "p", "H", "t"))
        if ex: a(f"| 단면 추가 | {ex} |")
        a(f"| 두께 | 공칭 t = {i.t_nominal:g}, 계산 t = {i.t_calc:g} mm |")
        a(f"| 쉐브론 각 | β = {i.beta_deg:g}° (유동축 기준) |")
        a(f"| 재료 | {i.raw.get('material', {}).get('name', '-')}: E = {i.E:g} MPa, ν = {i.nu:g}, G = {self.m.G:.0f} MPa, ρ = {i.rho:g} tonne/mm³ |")
        if i.S_allow: a(f"| 허용응력 | S_allow = {i.S_allow:g} MPa |")
        for n in i.notes: a(f"\n- 전처리: {n}")
        a("\n## 2. 형상 파생량 (docs/05 §3)\n")
        d = g.dimensionless()
        a("| 량 | 값 | 량 | 값 |\n|---|---|---|---|")
        a(f"| c = p/2 | {g.c:.4g} mm | f = H/2 | {g.f:.4g} mm |")
        a(f"| l (반주기 호장) | {g.l:.5g} mm | λ = l/c | {g.l / g.c:.5g} |")
        a(f"| J1 = ⟨1/√a⟩ | {g.avg_inv_sqrt_a:.5g} | J2 = ε²⟨φ²√a⟩ | {g.eps2_avg_phi2_sqrt_a:.5g} mm² |")
        a(f"| A_u (단위폭 단면적) | {g.A_u:.5g} mm²/mm | I_u | {g.I_u:.5g} mm⁴/mm |")
        a(f"| H/p | {d['H_over_p']:.3f} | t/p | {d['t_over_p']:.4f} |")
        if "t_over_Rmin" in d: a(f"| t/R_min | {d['t_over_Rmin']:.3f} | R_min | {g.extra['R_min']:.4g} mm |")
        for k, v in g.extra.items():
            if k not in ("R_min",): a(f"| {k} | {v:.5g} | | |")
        a("\n## 3. 국부축 등가 강성 (x ⊥ 능선, y ∥ 능선; N = A ε, M = D κ, κ_xy 공학 비틀림)\n")
        a("### 3.1 채택식 (Ye 2014 Eq. 19 기반, docs/06 §2)\n")
        a("| A11 [N/mm] | A12 | A22 | A66 | D11 [N·mm] | D12 | D22 | D66 |\n|---|---|---|---|---|---|---|---|")
        a("| " + " | ".join(f"{K[k]:.5g}" for k in KEYS_A[:4] + KEYS_D[:4]) + " |")
        a(f"\n직교이방성 비: A22/A11 = {K['A22'] / K['A11']:.1f}, D22/D11 = {K['D22'] / K['D11']:.1f}, A66/A11 = {K['A66'] / K['A11']:.2f}\n")
        a("### 3.2 모델 비교 (채택식 대비 %)\n")
        a("| 모델 | " + " | ".join(KEYS_A[:4] + KEYS_D[:4]) + " |\n|---|" + "---|" * 8)
        for n, v in self.models.items():
            a(f"| {n} | " + " | ".join(f"{(v[k] / K[k] - 1) * 100:+.1f}" for k in KEYS_A[:4] + KEYS_D[:4]) + " |")
        a("\n## 4. 등가 두께와 균질 단일층 공학상수 (docs/09 §1)\n")
        a("| 정의 | 값 [mm] | 용도 |\n|---|---|---|")
        a(f"| t_b = √(12D11/A11) | {th['t_b']:.4g} | 굽힘–막 정합 (11/12/22 공통, 편차 {th['t_b_spread']:.1e}) |")
        a(f"| t_s = √(12D66/A66) = hλ | {th['t_s']:.4g} | 전단–비틀림 정합 |")
        a(f"| t_m = hλ | {th['t_m']:.4g} | 질량·전개면적 |")
        hA, hD = self.homog["A"], self.homog["D"]
        a(f"\n균질 직교이방성 단일층(t_e = t_b, 간이 검토용): E_x = {hA['E_x']:.4g} MPa, E_y = {hA['E_y']:.5g} MPa, ν_xy = {hA['nu_xy']:.4f}, "
          f"G_xy = {hA['G_from_A']:.4g} (A66 기준, D66 오차 {hA['resid']['D66'] * 100:+.0f} %) 또는 {hD['G_from_D']:.4g} MPa (D66 기준, A66 오차 {hD['resid']['A66'] * 100:+.0f} %). "
          "→ 조립체 FE는 preintegrated section 사용.\n")
        a(f"## 5. 판축 영역 강성 (β = {i.beta_deg:g}°, docs/07)\n")
        a("| 성분 | 국부축 | 우측 영역 (θ = −β) | 좌측 영역 (θ = +β) | Strip 조합 | Voigt | Reuss |\n|---|---|---|---|---|---|---|")
        for k in KEYS_A + KEYS_D:
            a(f"| {k} | {K.get(k, 0.0):.5g} | {self.zones['right'][k]:.5g} | {self.zones['left'][k]:.5g} | {self.homogenized['strip'].get(k, 0.0):.5g} | "
              f"{self.homogenized['voigt'].get(k, 0.0):.5g} | {self.homogenized['reuss'].get(k, 0.0):.5g} |")
        a("\nFE 입력(`sections.inp`): 영역별 GENS 절 2개(판축, 16/26 포함). 횡전단 추정 "
          f"E11 = E22 = {self.E_ts['E11']:.4g} N/mm, 밀도(단위두께) {self.dens_area:.4g} tonne/mm². 비틀림 D33 규약은 E-1 단일 요소 시험으로 확인.\n")
        a("## 6. 응력 복원 — 단위 합력 응답 (국부축, VAM 복원식 docs/09 §3)\n")
        a(f"능선 가로 인장 응력집중계수 K_t = 1 + 6f/t = {self.kt:.2f}\n")
        a("| 합력 (단위) | max σ_vM / 단위 [MPa] | 위치 x, z [mm] | 면 (+z/−z) | 능선 σ_s +z/−z | 능선 σ_y +z/−z |\n|---|---|---|---|---|---|")
        for u in self.unit_response:
            a(f"| {u['resultant']} = 1 {u['unit']} | {u['max_vm_per_unit']:.4g} | {u['at_x']:.2f}, {u['at_z']:+.2f} | {u['fibre']} | "
              f"{u['ridge_sig_s'][0]:.3g} / {u['ridge_sig_s'][1]:.3g} | {u['ridge_sig_y'][0]:.3g} / {u['ridge_sig_y'][1]:.3g} |")
        if self.load_results:
            a("\n## 7. 하중 케이스\n")
            a("| 케이스 | 축/영역 | N (국부) [N/mm] | M (국부) [N·mm/mm] | max σ_vM [MPa] | 위치 | 막 σ_vM max | 허용 | 비 | 판정 |\n|---|---|---|---|---|---|---|---|---|---|")
            for r in self.load_results:
                a(f"| {r.name} | {r.axes}/{r.zone} | {', '.join(f'{v:.3g}' for v in r.N_local)} | {', '.join(f'{v:.3g}' for v in r.M_local)} | "
                  f"{r.max_vm:.4g} | x={r.at_x:.2f}, z={r.at_z:+.2f}, {r.fibre} | {r.membrane_vm_max:.4g} | "
                  f"{'-' if r.allowable is None else f'{r.allowable:g}'} | {'-' if r.ratio is None else f'{r.ratio:.2f}'} | {r.verdict or '-'} |")
            a("\n능선(산) 상세 (+z면/−z면; +z = 판의 산 쪽 표면): σ_s, σ_y, τ [MPa]; 국부 합력 N_s, N_y [N/mm], M_s, M_y [N·mm/mm]\n")
            a("| 케이스 | σ_s | σ_y | τ | N_s | N_y | M_s | M_y |\n|---|---|---|---|---|---|---|---|")
            for r in self.load_results:
                a(f"| {r.name} | {r.ridge_sig_s[0]:.3g} / {r.ridge_sig_s[1]:.3g} | {r.ridge_sig_y[0]:.3g} / {r.ridge_sig_y[1]:.3g} | "
                  f"{r.ridge_tau[0]:.3g} / {r.ridge_tau[1]:.3g} | {r.ridge_N_s:.4g} | {r.ridge_N_y:.4g} | {r.ridge_M_s:.4g} | {r.ridge_M_y:.4g} |")
            a("\n판정은 선형 복원 응력(막+굽힘, 1차 국부)의 von Mises 대 허용응력 단순 비교이며, apex·테두리·포트 보정계수(D단계)와 t/R 보정은 포함하지 않는다.")
        a("\n## 8. 한계\n")
        a("- 판 본체(주기 영역) 전용. apex·테두리·포트·용접부는 상세 FE 보정계수 필요 (docs/08 §2).")
        a("- 선형 탄성, Kirchhoff 등가 판. 국부 좌굴·접촉·소성은 범위 밖.")
        a("- 판축 합력은 사용자가 별도 산정(폐형식 또는 등가 판 FE)하여 `load_cases` 로 입력한다.")
        return "\n".join(L) + "\n"


# ============================================================================= 실행
def load_yaml(path: str) -> Dict[str, Any]:
    if yaml is None:
        raise RuntimeError("pyyaml 필요: pip install pyyaml")
    with open(path, encoding="utf-8") as f:
        return yaml.safe_load(f)


def run(path: str, out_dir: str) -> PlateModel:
    pm = PlateModel(prepare(load_yaml(path)))
    os.makedirs(out_dir, exist_ok=True)
    with open(os.path.join(out_dir, "report.md"), "w", encoding="utf-8") as f: f.write(pm.report_md())
    with open(os.path.join(out_dir, "results.json"), "w", encoding="utf-8") as f: json.dump(pm.results(), f, ensure_ascii=False, indent=1, default=float)
    with open(os.path.join(out_dir, "sections.inp"), "w", encoding="utf-8") as f: f.write(pm.apdl() + "\n")
    return pm


def _selfcheck() -> None:
    here = os.path.dirname(os.path.abspath(__file__))
    raw = load_yaml(os.path.join(here, "plate_input_template.yaml"))
    pm = PlateModel(prepare(raw))
    K = pm.K
    assert abs(K["A11"] / 2603.0 - 1) < 2e-3 and abs(K["D22"] / 170110 - 1) < 2e-3, K   # docs/06 §5 기본값
    assert abs(pm.th["t_b"] - 3.734) < 2e-3
    # 전처리 규칙
    raw2 = json.loads(json.dumps(raw)); raw2["section"]["dims_reference"] = "outer_surface"; raw2["section"]["H"] = 3.8
    assert abs(prepare(raw2).H - 3.2) < 1e-12
    raw2 = json.loads(json.dumps(raw)); raw2["chevron"]["beta_reference"] = "horizontal"; raw2["chevron"]["beta_deg"] = 30.0
    assert abs(prepare(raw2).beta_deg - 60.0) < 1e-12
    raw2 = json.loads(json.dumps(raw)); raw2["section"]["pitch_reference"] = "along_plate_axis"; raw2["section"]["p"] = 18.0
    assert abs(prepare(raw2).p - 9.0) < 1e-9
    raw2 = json.loads(json.dumps(raw)); raw2["material"]["t_min"] = 0.55; raw2["material"]["corrosion_allowance"] = 0.05
    p2 = prepare(raw2); assert abs(p2.t_calc - 0.5) < 1e-12 and PlateModel(p2).g.h == 0.5
    # 하중 케이스: 판축 β=0 → 국부축과 동일 / 단위 응답 일관성
    raw2 = json.loads(json.dumps(raw)); raw2["chevron"]["beta_deg"] = 0.0
    raw2["load_cases"] = [dict(name="a", axes="plate", zone="right", N=[10, 0, 0], M=[0, 0, 0]),
                          dict(name="b", axes="local", N=[10, 0, 0], M=[0, 0, 0], allowable=200.0)]
    pm2 = PlateModel(prepare(raw2)); ra, rb = pm2.load_results
    assert abs(ra.max_vm - rb.max_vm) < 1e-9 and rb.verdict == "NG" and abs(rb.ratio - rb.max_vm / 200) < 1e-12
    u = {x["resultant"]: x["max_vm_per_unit"] for x in pm2.unit_response}
    assert abs(u["N_x"] * 10 - rb.max_vm) < 1e-6
    assert abs(rb.ridge_sig_s[1] / (10 / pm2.g.h) - pm2.kt) < 1e-3
    # 판축 케이스가 국부 변환을 거치는지 (β=60, N_xp 만 → 국부 N_x, N_y, N_xy 모두 비영)
    raw3 = json.loads(json.dumps(raw)); raw3["load_cases"] = [dict(name="c", axes="plate", zone="left", N=[10, 0, 0])]
    rc = PlateModel(prepare(raw3)).load_results[0]; assert all(abs(v) > 1e-6 for v in rc.N_local)
    out = os.path.join(here, "..", "out", "_selfcheck"); run(os.path.join(here, "plate_input_template.yaml"), out)
    for fn in ("report.md", "results.json", "sections.inp"): assert os.path.getsize(os.path.join(out, fn)) > 500
    import shutil; shutil.rmtree(out)
    print("plate_model selfcheck OK")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", nargs="?", help="YAML 입력 파일 (calc/plate_input_template.yaml 형식)")
    ap.add_argument("-o", "--out", default=None, help="출력 디렉터리 (기본: out/<plate_id>)")
    ap.add_argument("--selfcheck", action="store_true")
    args = ap.parse_args()
    if args.selfcheck or not args.input:
        _selfcheck()
    else:
        pm_ = PlateModel(prepare(load_yaml(args.input)))
        out_ = args.out or os.path.join("out", pm_.inp.plate_id)
        run(args.input, out_)
        print(pm_.report_md())
        print(f"[written] {out_}/report.md, results.json, sections.inp")
