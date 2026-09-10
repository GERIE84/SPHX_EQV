"""C-2: 엑셀 계산 시트 생성기 → calc/SPHX_EQV_calc.xlsx  (설계 현장용, 수식 셀 = 라이브 계산).

plate_model.py 와 같은 로직을 엑셀 수식으로 구현한다 (원호+접선 단면, 채택식 강성, 등가 두께, ±β 영역 강성,
ANSYS GENS 문자열, 하중 케이스 → 능선(산·골) VAM 응력 복원·판정).  범위 차이:
  - 단면: arc_tangent 전용 (모드 A: R_c[,R_v] 입력 → α 를 Newton 반복 20회로 해결; 모드 B: α 입력, R_c=R_v).
    sinusoidal / trapezoidal / round 은 파이썬 전용.
  - 응력 복원: 산(crest)·골(valley) 두 점 × 외/내 면 4곳만 (플랭크 전체 분포는 파이썬 전용). 단위 응답 검증에서
    최대 von Mises 는 6 하중 성분 모두 산 또는 골에서 발생함을 확인했다 (docs/example_report).
사용:  python calc/make_excel_sheet.py            → 시트 생성 후 LibreOffice 재계산, 파이썬 결과와 대조
"""
from __future__ import annotations

import json
import math
import os
import re
import subprocess
import sys

from openpyxl import Workbook, load_workbook
from openpyxl.comments import Comment
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
from openpyxl.utils import get_column_letter
from openpyxl.workbook.defined_name import DefinedName
from openpyxl.worksheet.datavalidation import DataValidation

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "SPHX_EQV_calc.xlsx")
REF = os.path.join(HERE, "..", "docs", "example_report", "results.json")

FONT = "Arial"
F_BASE = Font(name=FONT, size=10)
F_BOLD = Font(name=FONT, size=10, bold=True)
F_H1 = Font(name=FONT, size=14, bold=True)
F_H2 = Font(name=FONT, size=11, bold=True, color="FFFFFF")
F_IN = Font(name=FONT, size=10, color="0000FF")
F_LINK = Font(name=FONT, size=10, color="008000")
F_NOTE = Font(name=FONT, size=9, italic=True, color="666666")
FILL_IN = PatternFill("solid", fgColor="FFFF99")
FILL_H2 = PatternFill("solid", fgColor="1F4E78")
FILL_OUT = PatternFill("solid", fgColor="E2EFDA")
FILL_WARN = PatternFill("solid", fgColor="FCE4D6")
THIN = Side(style="thin", color="BFBFBF")
BORDER = Border(left=THIN, right=THIN, top=THIN, bottom=THIN)

# 정의 이름 규칙: 셀 주소와 같은 이름(A11, D66, I1, J2, C1, Z1 …) 은 엑셀에서 불가, 이름은 대소문자 구분 없음(H_ ≠ h_ 불가).
# 소스에서는 읽기 쉬운 짧은 이름을 쓰고, 셀에 쓸 때 아래 맵으로 치환한다.
RENAME = {"H_": "H_mid", "h_": "h_calc", "q_": "q_int", "I1": "I_1", "I2": "I_2", "J1": "J_1", "J2": "J_2",
          "I2c": "I_2c", "I2v": "I_2v", "I2f": "I_2f", "c1": "c_1", "c2": "c_2", "c4": "c_4", "z1": "z_1", "z2": "z_2"}
for _pre in "AD":
    for _ij in ("11", "12", "22", "66", "16", "26"):
        RENAME[f"{_pre}{_ij}"] = f"{_pre}_{_ij}"
        for _z in "RL":
            RENAME[f"{_pre}{_ij}_{_z}"] = f"{_pre}_{_ij}_{_z}"
_RE_NAME = re.compile(r"\b(" + "|".join(sorted(map(re.escape, RENAME), key=len, reverse=True)) + r")\b")


def fix_name(name: str) -> str:
    return RENAME.get(name, name)


def fix_formula(v):
    """주의: 이름 치환은 토큰 단위이므로 수식 안의 상대 셀참조(D22 등)와 충돌한다 → 생성기에서는 셀참조를 항상 $D$22 절대형으로 쓴다."""
    """수식 문자열 안의 이름 토큰을 치환 (문자열 리터럴 내부는 건드리지 않음)."""
    if not (isinstance(v, str) and v.startswith("=")):
        return v
    parts = v.split('"')
    for i in range(0, len(parts), 2):          # 짝수 인덱스 = 문자열 리터럴 밖
        parts[i] = _RE_NAME.sub(lambda m: RENAME[m.group(1)], parts[i])
    return '"'.join(parts)


FMT_MM = "0.0000"
FMT_K = "#,##0.0"
FMT_SCI = "0.0000E+00"
FMT_G = "General"


class SheetBuilder:
    """라벨 | 값/수식 | 단위 | 비고 4열 표를 한 행씩 추가하고, 값 셀에 정의 이름을 등록한다."""

    def __init__(self, wb: Workbook, title: str, heading: str, widths=(34, 18, 12, 90)):
        self.wb, self.ws, self.row = wb, wb.create_sheet(title), 1
        for i, w in enumerate(widths, 1):
            self.ws.column_dimensions[get_column_letter(i)].width = w
        self.cell(1, 1, heading, F_H1); self.row = 3

    def cell(self, r, c, v, font=F_BASE, fmt=None, fill=None, border=False, align=None):
        x = self.ws.cell(row=r, column=c, value=fix_formula(v)); x.font = font
        if fmt: x.number_format = fmt
        if fill: x.fill = fill
        if border: x.border = BORDER
        if align: x.alignment = align
        return x

    def section(self, text):
        self.row += 1
        for c in range(1, 5):
            self.cell(self.row, c, text if c == 1 else None, F_H2, fill=FILL_H2)
        self.row += 1
        for c, t in enumerate(("항목", "값", "단위", "비고 / 수식 근거"), 1):
            self.cell(self.row, c, t, F_BOLD, border=True)
        self.row += 1

    def add(self, name, label, value, unit="", note="", fmt=FMT_G, kind="formula", comment=None):
        """kind: 'input' (파랑+노랑), 'formula' (검정), 'link' (녹색: 다른 시트 참조), 'text'."""
        r = self.row
        self.cell(r, 1, label, F_BASE, border=True)
        font = {"input": F_IN, "link": F_LINK}.get(kind, F_BASE)
        v = self.cell(r, 2, value, font, fmt=fmt, fill=FILL_IN if kind == "input" else None, border=True)
        self.cell(r, 3, unit, F_BASE, border=True); self.cell(r, 4, note, F_NOTE, border=True)
        if comment: v.comment = Comment(comment, "SPHX_EQV")
        if name:
            self.wb.defined_names[fix_name(name)] = DefinedName(fix_name(name), attr_text=f"'{self.ws.title}'!$B${r}")
        self.row += 1
        return f"'{self.ws.title}'!$B${r}"

    def note(self, text):
        self.cell(self.row, 1, text, F_NOTE); self.row += 1

    def blank(self, n=1):
        self.row += n


def build() -> Workbook:
    wb = Workbook(); wb.remove(wb.active)

    # ------------------------------------------------------------------ README
    rd = wb.create_sheet("README")
    rd.column_dimensions["A"].width = 120
    lines = [
        ("SPHX_EQV — 쉐브론 전열판 등가 판 계산 시트 (C-2)", F_H1),
        ("", F_BASE),
        ("목적: 원호+접선 주름 단면의 치수와 재료를 입력하면 등가 직교이방성 판 강성(채택식), 등가 두께, ±β 영역 강성, ANSYS 절 입력 문자열,", F_BASE),
        ("      하중 케이스에 대한 능선(산·골) 응력과 허용응력 판정을 수식으로 계산한다. 파이썬 모듈 calc/plate_model.py 와 같은 로직이다.", F_BASE),
        ("", F_BASE),
        ("사용법", F_BOLD),
        ("  1. 'Input' 시트의 노란 셀(파란 글씨)만 채운다. 나머지 시트는 전부 수식이다 (검정 = 수식, 녹색 = 다른 시트 참조).", F_BASE),
        ("  2. 단면 입력 모드: A = R_c(,R_v) 입력 → 플랭크 각 α 를 Newton 반복으로 계산 / B = α 입력, R_c=R_v 계산.", F_BASE),
        ("  3. 'Geometry' 상단 경고 행과 'Check' 시트(파이썬 대조)를 확인한다. Check 의 오차는 기본값 입력에서만 의미가 있다.", F_BASE),
        ("  4. 'Rotation' 시트 하단의 SSPA/SSPD/SSPE/SSPM 문자열을 ANSYS 입력파일에 복사한다 (비틀림 규약은 E-1 확인 전까지 'verify').", F_BASE),
        ("  5. 'Stress' 시트에 판축 또는 국부축 합력을 넣으면 산·골 4개 표면점의 응력과 최대 von Mises, 허용응력 비가 나온다.", F_BASE),
        ("", F_BASE),
        ("범위·한계", F_BOLD),
        ("  - 단면은 원호+접선(arc_tangent)만. 사인형/사다리꼴/반원 단면은 파이썬 모듈 사용.", F_BASE),
        ("  - 응력은 산(crest)·골(valley) 위치만 복원한다. 6개 단위 하중 모두에서 최대 응력이 산 또는 골에 생기지만, 플랭크 분포는 파이썬 전용.", F_BASE),
        ("  - 판 본체(주기 영역) 전용. apex·테두리·포트 보정계수, t/R_min>0.2 보정, 좌굴·접촉은 범위 밖 (docs/08, 09).", F_BASE),
        ("  - 판축 합력(N, M)은 사용자가 별도 산정(폐형식 또는 등가 판 FE)하여 입력한다.", F_BASE),
        ("", F_BASE),
        ("근거 문서: docs/05 (형상 정의), 06 (채택식 = Ye 2014 Eq. 19), 07 (β 회전), 09 (등가 두께·GENS·VAM 복원), 10 (계산 모듈).", F_BASE),
        ("단위: mm, N, MPa, N/mm (A, N), N·mm/mm (D, M), tonne/mm³ (ρ). 각도 입력은 도(°), 내부 계산은 rad.", F_BASE),
        ("기호: c=p/2, f=H/2, h=t_calc, λ=l/c, J1=I1/(2c), J2=I2/(2c). 국부축 x ⊥ 능선, y ∥ 능선. N=A·{εx,εy,γxy}, M=D·{κx,κy,2κxy}.", F_BASE),
        ("생성: calc/make_excel_sheet.py (재생성 시 이 파일을 덮어쓴다 — 수정은 생성기에서).", F_NOTE),
    ]
    for i, (t, f) in enumerate(lines, 1):
        rd.cell(row=i, column=1, value=t).font = f

    # ------------------------------------------------------------------ Input
    s = SheetBuilder(wb, "Input", "입력 (노란 셀만 편집)")
    s.section("판 식별")
    s.add("in_id", "판 식별자", "SPHX-example", "", "도면 번호 등", kind="input")
    s.add(None, "값 출처", "assumed defaults (docs/05 §5)", "", "예시값은 05 §5 가정 기본값 — 실제 도면값으로 교체", kind="input")
    s.section("단면 (docs/05 §2.1, 중앙면 기준)")
    s.add("in_p", "피치 p", 9.0, "mm", "능선 수직 산-산 거리 (pitch_reference 에 따라 환산)", FMT_MM, "input")
    s.add("in_H", "깊이 H", 3.2, "mm", "산-골 높이차 (dims_reference 에 따라 환산)", FMT_MM, "input")
    s.add("in_t", "공칭 두께 t", 0.6, "mm", "", FMT_MM, "input")
    s.add("in_mode", "단면 입력 모드", "A", "", "A: R_c(,R_v) 입력 → α 계산 | B: α 입력 → R_c=R_v 계산", kind="input")
    s.add("in_Rc", "산 원호 반경 R_c", 1.5, "mm", "모드 A. 외면 반경으로 측정했으면 R_o − t/2 로 환산해 기입", FMT_MM, "input")
    s.add("in_Rv", "골 원호 반경 R_v", None, "mm", "모드 A. 비우면 R_c 와 같게 처리", FMT_MM, "input")
    s.add("in_alpha", "플랭크 각 α (입력)", None, "°", "모드 B 전용", FMT_MM, "input")
    s.add("in_dims", "치수 기준 dims_reference", "mid", "", "mid: 중앙면 | outer: 외면 (H = H_in − t 로 환산)", kind="input")
    s.add("in_pitchref", "피치 기준 pitch_reference", "normal", "", "normal: 능선 수직 | axis: x_p 방향 측정 (p = p_in·cosβ)", kind="input")
    s.section("쉐브론 (docs/05 §2.2)")
    s.add("in_beta", "쉐브론 각 β (입력)", 60.0, "°", "", FMT_MM, "input")
    s.add("in_betaref", "β 기준 beta_reference", "flow", "", "flow: 유동축 기준 | horizontal: 수평 기준 (β = 90 − 입력)", kind="input")
    s.section("재료 (docs/05 §2.5)")
    s.add(None, "재질", "SA-240 316L", "", "", kind="input")
    s.add("in_E", "탄성계수 E", 193000.0, "MPa", "설계온도 값 사용", FMT_K, "input")
    s.add("in_nu", "포아송비 ν", 0.3, "", "", FMT_MM, "input")
    s.add("in_rho", "밀도 ρ", 7.9e-9, "tonne/mm³", "SSPM 입력용", FMT_SCI, "input")
    s.add("in_tmin", "성형 후 최소 두께 t_min", None, "mm", "비우면 공칭 t 사용 (경고 표시)", FMT_MM, "input")
    s.add("in_CA", "부식 여유 CA", 0.0, "mm", "", FMT_MM, "input")
    s.add("in_Sallow", "허용응력 S_allow", None, "MPa", "설계코드·온도 기준. 비우면 판정 생략", FMT_K, "input")
    s.section("하중 케이스 (docs/10 §2.2) — 합력은 별도 산정값")
    s.add("lc_name", "케이스 이름", "plate-axis My (right zone)", "", "", kind="input")
    s.add("lc_axes", "합력 좌표축 axes", "plate", "", "plate: 판축 (x_p, y_p) | local: 국부축 (x ⊥ 능선)", kind="input")
    s.add("lc_zone", "영역 zone", "right", "", "right: θ = −β | left: θ = +β (plate 축일 때만 사용)", kind="input")
    s.add("lc_Nx", "N_x (또는 N_xp)", 0.0, "N/mm", "", FMT_MM, "input")
    s.add("lc_Ny", "N_y (또는 N_yp)", 0.0, "N/mm", "", FMT_MM, "input")
    s.add("lc_Nxy", "N_xy", 0.0, "N/mm", "", FMT_MM, "input")
    s.add("lc_Mx", "M_x (또는 M_xp)", 0.0, "N·mm/mm", "", FMT_MM, "input")
    s.add("lc_My", "M_y (또는 M_yp)", 100.0, "N·mm/mm", "", FMT_MM, "input")
    s.add("lc_Mxy", "M_xy", 0.0, "N·mm/mm", "공학 비틀림 κ_xy 에 대응 (M_xy = D66·2κ_xy)", FMT_MM, "input")
    s.add("lc_allow", "케이스 허용응력", None, "MPa", "비우면 S_allow 사용", FMT_K, "input")
    ws = s.ws
    for name, rng, opts in (("in_mode", "B9", '"A,B"'), ("in_dims", "B13", '"mid,outer"'), ("in_pitchref", "B14", '"normal,axis"'),
                            ("in_betaref", "B18", '"flow,horizontal"'), ("lc_axes", "B33", '"plate,local"'), ("lc_zone", "B34", '"right,left"')):
        pass  # 드롭다운은 아래에서 정의 이름 → 셀 주소로 추가
    for name, opts in (("in_mode", '"A,B"'), ("in_dims", '"mid,outer"'), ("in_pitchref", '"normal,axis"'),
                       ("in_betaref", '"flow,horizontal"'), ("lc_axes", '"plate,local"'), ("lc_zone", '"right,left"')):
        addr = wb.defined_names[fix_name(name)].attr_text.split("!")[1].replace("$", "")
        dv = DataValidation(type="list", formula1=opts, allow_blank=False); ws.add_data_validation(dv); dv.add(addr)

    # ------------------------------------------------------------------ Geometry
    g = SheetBuilder(wb, "Geometry", "형상 전처리·파생량 (docs/05 §3, 원호+접선 폐형식)")
    g.section("전처리 (표기 기준 → 계산 기준, docs/10 §2.1)")
    g.add("beta_deg", "β (유동축 기준)", '=IF(in_betaref="horizontal",90-in_beta,in_beta)', "°", "beta_reference=horizontal → 90 − 입력")
    g.add("beta", "β [rad]", "=RADIANS(beta_deg)", "rad")
    g.add("p_", "피치 p (능선 수직)", '=IF(in_pitchref="axis",in_p*COS(beta),in_p)', "mm", "pitch_reference=axis → p_in·cosβ")
    g.add("H_", "깊이 H (중앙면)", '=IF(in_dims="outer",in_H-in_t,in_H)', "mm", "dims_reference=outer → H_in − t")
    g.add("h_", "계산 두께 h = t_calc", "=IF(ISBLANK(in_tmin),in_t,in_tmin)-in_CA", "mm", "t_min(없으면 t) − CA. 강성·응력 모두 이 값 사용")
    g.add("c_", "반피치 c = p/2", "=p_/2", "mm")
    g.add("f_", "반깊이 f = H/2", "=H_/2", "mm")
    g.section("원호+접선 단면 해 (모드 A: Newton 반복 / 모드 B: 선형)")
    g.add("Rc_in", "R_c (모드 A)", "=in_Rc", "mm", kind="link")
    g.add("Rv_in", "R_v (모드 A)", "=IF(ISBLANK(in_Rv),in_Rc,in_Rv)", "mm", "비우면 R_c")
    g.add("Rs", "R_s = R_c + R_v", "=Rc_in+Rv_in", "mm")
    g.add("a_hi", "α 상한 (T_L ≥ 0)", "=MIN(ASIN(MIN(1,c_/Rs)),ACOS(MAX(-1,1-H_/Rs)),PI()/2)-0.000001", "rad", "수평식·수직식에서 T_L ≥ 0 이 되는 α 범위의 상한")
    g.note("Newton 반복: F(α) = R_s(1−cosα) + (c − R_s sinα)·tanα − H = 0,  F'(α) = (c − R_s sinα)/cos²α.  각 반복값은 (0.001, a_hi) 로 제한.")
    r0 = g.row
    g.cell(r0, 1, "반복 n", F_BOLD, border=True); g.cell(r0, 2, "α_n [rad]", F_BOLD, border=True)
    g.cell(r0, 3, "F(α_n)", F_BOLD, border=True); g.cell(r0, 4, "F'(α_n)", F_BOLD, border=True)
    g.row += 1
    n_iter = 20
    first = g.row
    for n in range(n_iter + 1):
        r = g.row
        g.cell(r, 1, n, F_BASE, border=True)
        if n == 0:
            g.cell(r, 2, "=MIN(a_hi,MAX(0.001,ATAN(H_/MAX(c_-Rs*0.5,0.001))))", F_BASE, FMT_MM, border=True)
        else:
            g.cell(r, 2, f"=MIN(a_hi,MAX(0.001,$B${r-1}-$C${r-1}/$D${r-1}))", F_BASE, FMT_MM, border=True)
        g.cell(r, 3, f"=Rs*(1-COS($B${r}))+(c_-Rs*SIN($B${r}))*TAN($B${r})-H_", F_BASE, FMT_SCI, border=True)
        g.cell(r, 4, f"=(c_-Rs*SIN($B${r}))/COS($B${r})^2", F_BASE, FMT_MM, border=True)
        g.row += 1
    last = g.row - 1
    g.blank()
    g.add("alpha_A", "α (모드 A 해)", f"=$B${last}", "rad", "Newton 20회 최종값")
    g.add("resid_A", "잔차 |F(α)| (모드 A)", f"=ABS($C${last})", "mm", "1e-9 이하이면 수렴. 크면 R_c 가 너무 커서 접선 구간이 없음")
    g.add("alpha_B", "α (모드 B 입력)", "=IF(ISBLANK(in_alpha),0,RADIANS(in_alpha))", "rad")
    g.add("R_B", "R (모드 B: R_c = R_v)", "=IF(alpha_B>0,(c_*SIN(alpha_B)-H_*COS(alpha_B))/(2*(1-COS(alpha_B))),0)", "mm",
          "2R sinα + T cosα = c,  2R(1−cosα) + T sinα = H  →  det = 2(1−cosα)")
    g.add("alpha", "α (채택)", '=IF(in_mode="B",alpha_B,alpha_A)', "rad")
    g.add("alpha_deg", "α [°]", "=DEGREES(alpha)", "°")
    g.add("Rc", "R_c (채택)", '=IF(in_mode="B",R_B,Rc_in)', "mm")
    g.add("Rv", "R_v (채택)", '=IF(in_mode="B",R_B,Rv_in)', "mm")
    g.add("TL", "접선(플랭크) 길이 T_L", "=MAX(0,(c_-(Rc+Rv)*SIN(alpha))/COS(alpha))", "mm", "수평식. 반주기당 직선 길이")
    g.add("Rmin", "R_min", "=MIN(Rc,Rv)", "mm")
    g.add("chk_geom", "형상 정합 검사 (수직식 잔차)", "=(Rc+Rv)*(1-COS(alpha))+TL*SIN(alpha)-H_", "mm", "0 이어야 함 (모드 A 수렴, 모드 B 물리해 확인)", FMT_SCI)
    g.section("파생량 (docs/05 §3, Lang & Su Eq. 27 일반화)")
    g.add("l_", "반주기 호장 l", "=(Rc+Rv)*alpha+TL", "mm")
    g.add("lam", "λ = l/c", "=l_/c_", "")
    g.add("q_", "q = ∫₀^α cos²θ dθ", "=alpha/2+SIN(2*alpha)/4", "rad")
    g.add("I1", "I1 = ∫(dx/ds)² ds (1주기)", "=2*((Rc+Rv)*q_+TL*COS(alpha)^2)", "mm")
    g.add("a_c", "산 원호 중심 |z| = f − R_c", "=f_-Rc", "mm")
    g.add("a_v", "골 원호 중심 |z| = f − R_v", "=f_-Rv", "mm")
    g.add("z1", "플랭크 시작 z1", "=f_-Rc*(1-COS(alpha))", "mm")
    g.add("z2", "플랭크 끝 z2", "=z1-TL*SIN(alpha)", "mm")
    g.add("I2c", "I2 산 원호", "=Rc*(a_c^2*alpha+2*a_c*Rc*SIN(alpha)+Rc^2*q_)", "mm³")
    g.add("I2v", "I2 골 원호", "=Rv*(a_v^2*alpha+2*a_v*Rv*SIN(alpha)+Rv^2*q_)", "mm³")
    g.add("I2f", "I2 플랭크", "=IF(TL>0,(z1^3-z2^3)/(3*SIN(alpha)),0)", "mm³")
    g.add("I2", "I2 = ∫ z² ds (1주기)", "=2*(I2c+I2f+I2v)", "mm³")
    g.add("J1", "J1 = I1/(2c) = ⟨1/√a⟩", "=I1/(2*c_)", "")
    g.add("J2", "J2 = I2/(2c) = ε²⟨φ²√a⟩", "=I2/(2*c_)", "mm²")
    g.add("A_u", "단위폭 단면적 A_u = hλ", "=h_*lam", "mm²/mm")
    g.add("I_u", "단위폭 2차모멘트 I_u = h·J2", "=h_*J2", "mm⁴/mm")
    g.add("W_u", "단위폭 단면계수 W_u", "=I_u/(f_+h_/2)", "mm³/mm")
    g.section("무차원·경고")
    g.add("Hp", "H/p", "=H_/p_", "")
    g.add("tp", "t/p", "=h_/p_", "")
    g.add("tR", "t/R_min", "=h_/Rmin", "")
    g.add("warn1", "경고: 두께", '=IF(ISBLANK(in_tmin),"t_min 미입력 → 공칭 두께로 계산 (성형 두께 감소 미반영)","-")', "", kind="text")
    g.add("warn2", "경고: t/R_min", '=IF(tR>0.2,"t/R_min > 0.2: 쉘 이론 범위 밖. 능선 응력은 D단계 솔리드 FE 보정 (09 §5)","-")', "", kind="text")
    g.add("warn3", "경고: H/p", '=IF(Hp>1,"H/p > 1: 문헌 검증 범위 밖","-")', "", kind="text")
    g.add("warn4", "경고: 수렴", '=IF(AND(in_mode="A",resid_A>0.000000001),"모드 A Newton 미수렴: R_c/R_v 가 커서 접선 구간이 없음","-")', "", kind="text")
    for nm in ("warn1", "warn2", "warn3", "warn4"):
        addr = wb.defined_names[fix_name(nm)].attr_text.split("!")[1].replace("$", "")
        g.ws[addr].fill = FILL_WARN

    # ------------------------------------------------------------------ Stiffness
    k = SheetBuilder(wb, "Stiffness", "국부축 등가 강성 — 채택식 (docs/06 §2, Ye 2014 Eq. 19) 및 등가 두께 (docs/09 §1)")
    k.section("재료")
    k.add("E_", "E", "=in_E", "MPa", kind="link"); k.add("nu_", "ν", "=in_nu", "", kind="link")
    k.add("G_", "G = E/(2(1+ν))", "=E_/(2*(1+nu_))", "MPa")
    k.add("Q_", "Q = E/(1−ν²)", "=E_/(1-nu_^2)", "MPa", "평면응력 강성")
    k.section("채택식 A (막) [N/mm], D (굽힘) [N·mm]  — x ⊥ 능선, y ∥ 능선")
    k.add("A11", "A11", "=E_*h_^3/((1-nu_^2)*(12*J2+h_^2*J1))", "N/mm", "E h³ / ((1−ν²)(12J2 + h²J1))", FMT_K)
    k.add("A12", "A12", "=nu_*A11", "N/mm", "ν A11", FMT_K)
    k.add("A22", "A22", "=E_*h_*lam+nu_^2*A11", "N/mm", "E h λ + ν² A11", FMT_K)
    k.add("A66", "A66", "=G_*h_/lam", "N/mm", "G h / λ", FMT_K)
    k.add("D11", "D11", "=E_*h_^3/(12*(1-nu_^2)*lam)", "N·mm", "E h³ / (12(1−ν²) λ)", FMT_K)
    k.add("D12", "D12", "=nu_*D11", "N·mm", "ν D11", FMT_K)
    k.add("D22", "D22", "=E_*h_*J2+E_*h_^3*J1/12+nu_^2*D11", "N·mm", "E h J2 + E h³ J1/12 + ν² D11", FMT_K)
    k.add("D66", "D66", "=G_*h_^3*lam/12", "N·mm", "G h³ λ / 12 (κ_xy 공학 비틀림 관례: M_xy = D66·2κ_xy)", FMT_K)
    k.add(None, "A22/A11", "=A22/A11", "", "직교이방성 비", FMT_K)
    k.add(None, "D22/D11", "=D22/D11", "", "", FMT_K)
    k.add(None, "A66/A11", "=A66/A11", "", "", FMT_K)
    k.section("모델 비교 (채택식 대비 %, docs/06 §3)")
    k.add(None, "Xia 2012 D22 (Dirichlet 상한)", "=((E_*h_*J2+E_*h_^3*J1/12)/(1-nu_^2)/D22-1)*100", "%", "1/(1−ν²) 과대", FMT_MM)
    k.add(None, "Lang & Su 2022 A22", "=((nu_^2*A11+(1-nu_^2)*E_*h_*lam)/A22-1)*100", "%", "평판 극한 불환원", FMT_MM)
    k.add(None, "Lang & Su 2022 A11", "=((E_*h_^3/(12*J2*(1-nu_^2)+h_^2*J1))/A11-1)*100", "%", "", FMT_MM)
    k.add(None, "Ye 얕은 주름식 A11", "=((E_*h_^3/(12*(1-nu_^2)*J2))/A11-1)*100", "%", "J1 항 생략", FMT_MM)
    k.add(None, "Ye 얕은 주름식 D22", "=((E_*h_*J2)/D22-1)*100", "%", "", FMT_MM)
    k.section("등가 두께 (docs/09 §1)")
    k.add("t_b", "t_b = √(12 D11/A11)", "=SQRT(12*D11/A11)", "mm", "굽힘–막 정합 (11/12/22 공통)", FMT_MM)
    k.add("t_b22", "√(12 D22/A22) (검사)", "=SQRT(12*D22/A22)", "mm", "t_b 와 같아야 함", FMT_MM)
    k.add("t_s", "t_s = √(12 D66/A66) = hλ", "=SQRT(12*D66/A66)", "mm", "전단–비틀림 정합", FMT_MM)
    k.add("t_m", "t_m = hλ", "=h_*lam", "mm", "질량·전개면적", FMT_MM)
    k.section("균질 직교이방성 단일층 (t_e = t_b, 간이 검토용 — 조립체 FE 에는 GENS 절 사용)")
    k.add("Ex_h", "E_x", "=(A11*A22-A12^2)/(A22*t_b)", "MPa", "", FMT_K)
    k.add("Ey_h", "E_y", "=(A11*A22-A12^2)/(A11*t_b)", "MPa", "", FMT_K)
    k.add("nuxy_h", "ν_xy", "=A12/A22", "", "", FMT_MM)
    k.add("GA_h", "G_xy (A66 기준)", "=A66/t_b", "MPa", "D66 오차 → 아래", FMT_K)
    k.add(None, "  → D66 재현 오차", "=(GA_h*t_b^3/12/D66-1)*100", "%", "", FMT_MM)
    k.add("GD_h", "G_xy (D66 기준)", "=12*D66/t_b^3", "MPa", "A66 오차 → 아래", FMT_K)
    k.add(None, "  → A66 재현 오차", "=(GD_h*t_b/A66-1)*100", "%", "", FMT_MM)
    k.section("응력집중 (Briassoulis 1986 B2)")
    k.add("Kt", "K_t = 1 + 6f/h (능선 가로 인장, 내측)", "=1+6*f_/h_", "", "σ_s,inner / (N_x/h)", FMT_MM)

    # ------------------------------------------------------------------ Rotation
    t = SheetBuilder(wb, "Rotation", "판축 영역 강성 (docs/07) 및 ANSYS preintegrated general shell section 입력 (docs/09 §2)", (34, 18, 18, 90))
    t.section("회전각")
    t.add("th_R", "우측 영역 θ_R = −β", "=-beta", "rad", "능선이 +β 방향인 영역")
    t.add("th_L", "좌측 영역 θ_L = +β", "=beta", "rad")
    t.note("Q̄ 폐형식 (Jones 1975): m=cosθ, n=sinθ.  11: q11m⁴+2(q12+2q66)n²m²+q22n⁴ ; 22: q11n⁴+…+q22m⁴ ; 12: (q11+q22−4q66)n²m²+q12(n⁴+m⁴)")
    t.note("66: (q11+q22−2q12−2q66)n²m²+q66(n⁴+m⁴) ; 16: (q11−q12−2q66)nm³+(q12−q22+2q66)n³m ; 26: (q11−q12−2q66)n³m+(q12−q22+2q66)nm³")
    t.blank()
    hdr = t.row
    for c, txt in enumerate(("성분", "우측 영역 (θ=−β)", "좌측 영역 (θ=+β)", "비고"), 1):
        t.cell(hdr, c, txt, F_BOLD, border=True)
    t.row += 1
    rot_names = {}
    for pre in ("A", "D"):
        q11, q12, q22, q66 = (f"{pre}11", f"{pre}12", f"{pre}22", f"{pre}66")
        forms = {
            "11": "{q11}*m^4+2*({q12}+2*{q66})*n^2*m^2+{q22}*n^4",
            "22": "{q11}*n^4+2*({q12}+2*{q66})*n^2*m^2+{q22}*m^4",
            "12": "({q11}+{q22}-4*{q66})*n^2*m^2+{q12}*(n^4+m^4)",
            "66": "({q11}+{q22}-2*{q12}-2*{q66})*n^2*m^2+{q66}*(n^4+m^4)",
            "16": "({q11}-{q12}-2*{q66})*n*m^3+({q12}-{q22}+2*{q66})*n^3*m",
            "26": "({q11}-{q12}-2*{q66})*n^3*m+({q12}-{q22}+2*{q66})*n*m^3",
        }
        for ij, fm in forms.items():
            r = t.row
            t.cell(r, 1, f"{pre}{ij}", F_BASE, border=True)
            for col, th in ((2, "th_R"), (3, "th_L")):
                expr = fm.format(q11=q11, q12=q12, q22=q22, q66=q66).replace("m", f"COS({th})").replace("n", f"SIN({th})")
                t.cell(r, col, "=" + expr, F_BASE, FMT_K, border=True)
            t.cell(r, 4, "N/mm" if pre == "A" else "N·mm", F_NOTE, border=True)
            wb.defined_names[fix_name(f"{pre}{ij}_R")] = DefinedName(fix_name(f"{pre}{ij}_R"), attr_text=f"'Rotation'!$B${r}")
            wb.defined_names[fix_name(f"{pre}{ij}_L")] = DefinedName(fix_name(f"{pre}{ij}_L"), attr_text=f"'Rotation'!$C${r}")
            t.row += 1
    t.section("횡전단·질량 (Mindlin 쉘 보조량, docs/09 §2 항 5–6)")
    t.add("E_ts", "횡전단 E11 = E22 = (5/6) G h λ", "=5/6*G_*h_*lam", "N/mm", "1차 추정 — E단계 ×0.1/×10 감도 확인", FMT_K)
    t.add("dens_a", "단위두께 가정 밀도 = ρ h λ", "=in_rho*h_*lam", "tonne/mm²", "SSPM", FMT_SCI)
    t.section("ANSYS APDL 문자열 (복사용; 인자 순서 SSPA/SSPD: 11,21,31,22,32,33 — 3 = xy; 판축 요소좌표계)")
    fmt = '"0.000000E+00"'
    def T(name): return f"TEXT({name},{fmt})"
    for zone, sfx, secid, thtxt in (("우측", "R", 11, "-beta"), ("좌측", "L", 12, "+beta")):
        t.add(None, f"{zone} SECTYPE", f'="SECTYPE,{secid},GENS,,SPHX_{sfx}   ! zone {sfx}, theta = {thtxt}, element CS = plate axes"', "", "", kind="text")
        t.add(None, f"{zone} SSPA", f'="SSPA,"&{T(f"A11_{sfx}")}&","&{T(f"A12_{sfx}")}&","&{T(f"A16_{sfx}")}&","&{T(f"A22_{sfx}")}&","&{T(f"A26_{sfx}")}&","&{T(f"A66_{sfx}")}', "", "A11,A21,A31,A22,A32,A33", kind="text")
        t.add(None, f"{zone} SSPB", '="SSPB,0,0,0,0,0,0   ! B = 0 (symmetric section)"', "", "", kind="text")
        t.add(None, f"{zone} SSPD", f'="SSPD,"&{T(f"D11_{sfx}")}&","&{T(f"D12_{sfx}")}&","&{T(f"D16_{sfx}")}&","&{T(f"D22_{sfx}")}&","&{T(f"D26_{sfx}")}&","&{T(f"D66_{sfx}")}&"   ! D33/16/26 twist convention: verify (E-1)"', "", "D11,D21,D31,D22,D32,D33", kind="text")
        t.add(None, f"{zone} SSPE", f'="SSPE,"&{T("E_ts")}&",0,"&{T("E_ts")}&"   ! transverse shear estimate"', "", "", kind="text")
        t.add(None, f"{zone} SSPM", f'="SSPM,"&{T("dens_a")}&"   ! density assuming unit thickness = rho*h*lambda"', "", "", kind="text")
        t.blank()

    # ------------------------------------------------------------------ Stress
    st = SheetBuilder(wb, "Stress", "하중 케이스 → 국부축 합력 → VAM 복원 (Ye 2014 Eq. 14, 20; docs/09 §3) — 산·골 응력과 판정", (36, 18, 18, 90))
    st.section("하중 케이스 합력 → 국부축 (transform.resultants_to_local, T_σ)")
    st.add("th_lc", "θ (plate 축일 때)", '=IF(lc_axes="plate",IF(lc_zone="right",-beta,beta),0)', "rad", "local 축이면 0 (변환 없음)")
    st.add("mm_", "m = cosθ", "=COS(th_lc)", ""); st.add("nn_", "n = sinθ", "=SIN(th_lc)", "")
    st.add("Nx", "N_x (국부)", "=lc_Nx*mm_^2+lc_Ny*nn_^2+2*lc_Nxy*mm_*nn_", "N/mm", "T_σ 1행: m², n², 2mn", FMT_MM)
    st.add("Ny", "N_y (국부)", "=lc_Nx*nn_^2+lc_Ny*mm_^2-2*lc_Nxy*mm_*nn_", "N/mm", "n², m², −2mn", FMT_MM)
    st.add("Nxy", "N_xy (국부)", "=-lc_Nx*mm_*nn_+lc_Ny*mm_*nn_+lc_Nxy*(mm_^2-nn_^2)", "N/mm", "−mn, mn, m²−n²", FMT_MM)
    st.add("Mx", "M_x (국부)", "=lc_Mx*mm_^2+lc_My*nn_^2+2*lc_Mxy*mm_*nn_", "N·mm/mm", "", FMT_MM)
    st.add("My", "M_y (국부)", "=lc_Mx*nn_^2+lc_My*mm_^2-2*lc_Mxy*mm_*nn_", "N·mm/mm", "", FMT_MM)
    st.add("Mxy", "M_xy (국부)", "=-lc_Mx*mm_*nn_+lc_My*mm_*nn_+lc_Mxy*(mm_^2-nn_^2)", "N·mm/mm", "", FMT_MM)
    st.section("거시 변형률 (채택식 역행렬)")
    st.add("detA", "det A (11,22 블록)", "=A11*A22-A12^2", "", "", FMT_K)
    st.add("exx", "ε_xx", "=(A22*Nx-A12*Ny)/detA", "", "", FMT_SCI)
    st.add("eyy", "ε_yy", "=(-A12*Nx+A11*Ny)/detA", "", "", FMT_SCI)
    st.add("gxy", "γ_xy = 2ε_xy", "=Nxy/A66", "", "", FMT_SCI)
    st.add("detD", "det D (11,22 블록)", "=D11*D22-D12^2", "", "", FMT_K)
    st.add("kxx", "κ_xx", "=(D22*Mx-D12*My)/detD", "1/mm", "", FMT_SCI)
    st.add("kyy", "κ_yy", "=(-D12*Mx+D11*My)/detD", "1/mm", "", FMT_SCI)
    st.add("kxy", "κ_xy (텐서)", "=Mxy/(2*D66)", "1/mm", "M_xy = D66·2κ_xy", FMT_SCI)
    st.section("VAM 복원 상수 (대칭 단면)")
    st.add("CC", "𝒞 = −(12J2/h² + J1)", "=-(12*J2/h_^2+J1)", "", "", FMT_MM)
    st.add("c1", "c1 = −(ε_xx + ν ε_yy)/𝒞", "=-(exx+nu_*eyy)/CC", "", "", FMT_SCI)
    st.add("c4", "c4 = (κ_xx + ν κ_yy)/λ", "=(kxx+nu_*kyy)/lam", "1/mm", "", FMT_SCI)
    st.add("ch_c", "c_h(산) = 1 + h²/(48R_c²)", "=1+h_^2/(48*Rc^2)", "", "두께 보정")
    st.add("ch_v", "c_h(골) = 1 + h²/(48R_v²)", "=1+h_^2/(48*Rv^2)", "")
    st.add("inva1", "⟨√a/c_h⟩·2c = 2R_cα/c_h,c + 2T_L + 2R_vα/c_h,v", "=2*Rc*alpha/ch_c+2*TL+2*Rv*alpha/ch_v", "mm", "1주기 적분")
    st.add("alpha1", "α1 = 2c / 위", "=2*c_/inva1", "", "")
    st.add("c2", "c2 = α1 γ_xy", "=alpha1*gxy", "", "", FMT_SCI)
    st.section("산(crest, z=+f, κ=−1/R_c) 과 골(valley, z=−f, κ=+1/R_v) 의 물리 변형률 (√a = 1)")
    hdr = st.row
    for c, txt in enumerate(("량", "산 (crest)", "골 (valley)", "식"), 1):
        st.cell(hdr, c, txt, F_BOLD, border=True)
    st.row += 1
    rows = [
        ("eps_s", "ε_s", "=c1-nu_*(eyy+f_*kyy)", "=c1-nu_*(eyy-f_*kyy)", "γ11/a = c1√a − νa(ε_yy + zκ_yy)"),
        ("eps_y", "ε_y", "=eyy+f_*kyy", "=eyy-f_*kyy", "γ22 = ε_yy + zκ_yy"),
        ("gam_sy", "γ_sy", "=(c2+h_^2*kxy/(12*Rc))/ch_c", "=(c2-h_^2*kxy/(12*Rv))/ch_v", "2γ12/√a = (√a c2 − h² z'' κ_xy/(12a))/c_h,  z''=∓1/R"),
        ("kap_s", "κ_s", "=12*c1*f_/h_^2+c4-nu_*kyy", "=-12*c1*f_/h_^2+c4-nu_*kyy", "ρ11/a = 12c1 z/h² + c4 − νκ_yy/√a"),
        ("kap_y", "κ_y", "=kyy", "=kyy", "ρ22 = κ_yy/√a"),
        ("kap_sy2", "2κ_sy", "=(-2*kxy-c2/(2*Rc))/ch_c", "=(-2*kxy+c2/(2*Rv))/ch_v", "2ρ12/√a = (−2√a κ_xy + z'' c2/(2a))/c_h"),
    ]
    for nm, lab, fc, fv, note in rows:
        r = st.row
        st.cell(r, 1, lab, F_BASE, border=True)
        st.cell(r, 2, fc, F_BASE, FMT_SCI, border=True); st.cell(r, 3, fv, F_BASE, FMT_SCI, border=True)
        st.cell(r, 4, note, F_NOTE, border=True)
        wb.defined_names[nm + "_c"] = DefinedName(nm + "_c", attr_text=f"'Stress'!$B${r}")
        wb.defined_names[nm + "_v"] = DefinedName(nm + "_v", attr_text=f"'Stress'!$C${r}")
        st.row += 1
    st.section("표면 응력 (ζ = +h/2 외측 / −h/2 내측;  e = ε − ζκ;  평면응력)  [MPa]")
    hdr = st.row
    cols = ("산 외측", "산 내측", "골 외측", "골 내측")
    st.ws.column_dimensions["E"].width = 16; st.ws.column_dimensions["F"].width = 16
    for c, txt in enumerate(("량",) + cols + ("식",), 1):
        st.cell(hdr, c, txt, F_BOLD, border=True)
    st.row += 1
    pts = (("c", "+"), ("c", "-"), ("v", "+"), ("v", "-"))  # (위치, ζ 부호)
    def zeta(sg): return f"({sg}h_/2)"
    srows = {
        "es": ("e_s", lambda p, sg: f"=eps_s_{p}-{zeta(sg)}*kap_s_{p}", "ε_s − ζκ_s"),
        "ey": ("e_y", lambda p, sg: f"=eps_y_{p}-{zeta(sg)}*kap_y_{p}", "ε_y − ζκ_y"),
        "gs": ("e_sy (공학)", lambda p, sg: f"=gam_sy_{p}-{zeta(sg)}*kap_sy2_{p}", "γ_sy − ζ·2κ_sy"),
        "ss": ("σ_s", lambda p, sg: None, "Q (e_s + ν e_y)"),
        "sy": ("σ_y", lambda p, sg: None, "Q (e_y + ν e_s)"),
        "tt": ("τ_sy", lambda p, sg: None, "G e_sy"),
        "vm": ("σ_vM", lambda p, sg: None, "√(σ_s² − σ_sσ_y + σ_y² + 3τ²)"),
    }
    addr = {}
    for key, (lab, fn, note) in srows.items():
        r = st.row
        st.cell(r, 1, lab, F_BOLD if key == "vm" else F_BASE, border=True)
        for i, (p, sg) in enumerate(pts):
            col = 2 + i; L = "$" + get_column_letter(col) + "$"
            if key in ("es", "ey", "gs"):
                f = fn(p, sg)
            elif key == "ss":
                f = f"=Q_*({L}{addr['es']}+nu_*{L}{addr['ey']})"
            elif key == "sy":
                f = f"=Q_*({L}{addr['ey']}+nu_*{L}{addr['es']})"
            elif key == "tt":
                f = f"=G_*{L}{addr['gs']}"
            else:
                f = f"=SQRT({L}{addr['ss']}^2-{L}{addr['ss']}*{L}{addr['sy']}+{L}{addr['sy']}^2+3*{L}{addr['tt']}^2)"
            st.cell(r, col, f, F_BASE, FMT_SCI if key in ("es", "ey", "gs") else FMT_K, border=True,
                    fill=FILL_OUT if key == "vm" else None)
        st.cell(r, 6, note, F_NOTE, border=True)
        addr[key] = r; st.row += 1
    vm_r = addr["vm"]
    st.section("결과 요약·판정")
    st.add("vm_max", "최대 von Mises (산·골 4점)", f"=MAX($B${vm_r}:$E${vm_r})", "MPa", "플랭크 분포는 파이썬 전용", FMT_K)
    st.add("vm_where", "최대 위치", f'=CHOOSE(MATCH(vm_max,$B${vm_r}:$E${vm_r},0),"산 외측","산 내측","골 외측","골 내측")', "", "", kind="text")
    st.add("Ns_c", "산 막력 N_s = Q h (ε_s + ν ε_y)", "=Q_*h_*(eps_s_c+nu_*eps_y_c)", "N/mm", "능선 국부 합력", FMT_MM)
    st.add("Ny_c", "산 막력 N_y", "=Q_*h_*(eps_y_c+nu_*eps_s_c)", "N/mm", "", FMT_MM)
    st.add("Ms_c", "산 모멘트 M_s = Q h³/12 (κ_s + ν κ_y)", "=Q_*h_^3/12*(kap_s_c+nu_*kap_y_c)", "N·mm/mm", "", FMT_MM)
    st.add("My_c", "산 모멘트 M_y", "=Q_*h_^3/12*(kap_y_c+nu_*kap_s_c)", "N·mm/mm", "", FMT_MM)
    st.add("vm_mem", "산 막 응력 von Mises", "=SQRT(Ns_c^2-Ns_c*Ny_c+Ny_c^2+3*(G_*h_*gam_sy_c)^2)/h_", "MPa", "1차 일반 막응력 판정용", FMT_K)
    st.add("allow_", "허용응력 (케이스 → S_allow)", "=IF(ISBLANK(lc_allow),IF(ISBLANK(in_Sallow),0,in_Sallow),lc_allow)", "MPa", "0 이면 판정 생략", FMT_K)
    st.add("ratio_", "σ_vM,max / 허용", '=IF(allow_>0,vm_max/allow_,"-")', "", "", FMT_MM)
    st.add("verdict", "판정", '=IF(allow_>0,IF(vm_max<=allow_,"OK","NG"),"허용응력 미입력")', "", "선형 복원 응력 단순 비교. apex·테두리·t/R 보정 미포함", kind="text")
    st.note("특수해 검사: N_x 만 → 산 내측 σ_s/(N_x/h) = K_t = 1+6f/h ;  κ_yy 만 → 산 N_y = E h f κ_yy (Briassoulis Eq. 12).")

    # ------------------------------------------------------------------ Check
    ref = json.load(open(REF, encoding="utf-8"))
    ck = SheetBuilder(wb, "Check", "파이썬 대조 (docs/example_report/results.json, 가정 기본값 입력에서만 유효)", (30, 18, 18, 14))
    ck.ws.column_dimensions["D"].width = 14; ck.ws.column_dimensions["E"].width = 60
    ck.section("기본값: p=9, H=3.2, t=0.6, R_c=1.5, β=60°, E=193000, ν=0.3, 하중 케이스 plate My=100 (right)")
    hdr = ck.row - 1
    for c, txt in enumerate(("항목", "엑셀", "파이썬 (하드코딩 참조값)", "상대오차", "출처"), 1):
        ck.cell(hdr, c, txt, F_BOLD, border=True)
    gref, kref, zr, th_ = ref["geometry"], ref["stiffness_local"]["adopted"], ref["zones"]["right"], ref["equivalent_thickness"]
    lc = ref["stress"]["load_cases"][1]
    items = [
        ("α [°]", "=alpha_deg", gref["extra"]["alpha_deg"], "geometry.arc_tangent"),
        ("T_L", "=TL", gref["extra"]["T_L"], ""), ("l", "=l_", gref["l"], ""), ("I1", "=I1", gref["I1"], ""), ("I2", "=I2", gref["I2"], ""),
        ("A11", "=A11", kref["A11"], "stiffness.adopted"), ("A22", "=A22", kref["A22"], ""), ("A66", "=A66", kref["A66"], ""),
        ("D11", "=D11", kref["D11"], ""), ("D22", "=D22", kref["D22"], ""), ("D66", "=D66", kref["D66"], ""),
        ("t_b", "=t_b", th_["t_b"], "equivalent_plate"), ("t_s", "=t_s", th_["t_s"], ""),
        ("A11 우측", "=A11_R", zr["A11"], "transform.zone_stiffness"), ("A16 우측", "=A16_R", zr["A16"], ""), ("A26 우측", "=A26_R", zr["A26"], ""),
        ("D12 우측", "=D12_R", zr["D12"], ""), ("D16 우측", "=D16_R", zr["D16"], ""), ("D66 우측", "=D66_R", zr["D66"], ""),
        ("M_x 국부", "=Mx", lc["M_local"][0], "resultants_to_local"), ("M_xy 국부", "=Mxy", lc["M_local"][2], ""),
        ("산 외측 σ_s", f"=Stress!$B${addr['ss']}", lc["ridge_sig_s"][0], "ye_recovery (능선)"),
        ("산 내측 σ_s", f"=Stress!$C${addr['ss']}", lc["ridge_sig_s"][1], ""),
        ("산 외측 σ_y", f"=Stress!$B${addr['sy']}", lc["ridge_sig_y"][0], ""), ("산 내측 σ_y", f"=Stress!$C${addr['sy']}", lc["ridge_sig_y"][1], ""),
        ("산 외측 τ", f"=Stress!$B${addr['tt']}", lc["ridge_tau"][0], ""), ("산 내측 τ", f"=Stress!$C${addr['tt']}", lc["ridge_tau"][1], ""),
        ("최대 von Mises", "=vm_max", lc["max_vm"], "파이썬 전 단면 최대 (산·골에서 발생)"),
    ]
    first = ck.row
    for lab, f, v, src in items:
        r = ck.row
        ck.cell(r, 1, lab, F_BASE, border=True); ck.cell(r, 2, f, F_LINK, FMT_G, border=True)
        ck.cell(r, 3, v, F_IN, FMT_G, border=True); ck.cell(r, 4, f"=IF(ABS($C${r})>1E-9,ABS($B${r}/$C${r}-1),ABS($B${r}-$C${r}))", F_BASE, FMT_SCI, border=True)
        ck.cell(r, 5, src or "〃", F_NOTE, border=True); ck.row += 1
    last = ck.row - 1
    ck.blank()
    ck.add("chk_max", "최대 상대오차", f"=MAX($D${first}:$D${last})", "", "기본값 입력에서 1e-6 이하이면 시트 로직이 파이썬과 일치", FMT_SCI)
    ck.add(None, "판정", '=IF(chk_max<0.000001,"일치","불일치 — 입력이 기본값과 다르거나 시트 수정됨")', "", "", kind="text")
    ck.note("참조값은 calc/plate_model.py 실행 결과(docs/example_report/results.json)를 생성 시점에 복사한 하드코딩 값이다 (파랑).")

    wb.calculation.fullCalcOnLoad = True          # openpyxl 은 캐시값을 쓰지 않으므로 열 때 전체 재계산
    wb.move_sheet("README", offset=-len(wb.sheetnames))
    for ws_ in wb.worksheets:
        ws_.sheet_view.showGridLines = False
        ws_.freeze_panes = None
    return wb


def recalc_and_verify(path: str) -> None:
    """이 환경에는 LibreOffice Calc 가 없어 `formulas` 라이브러리(pip)로 전체 워크북을 평가하고 Check 시트를 대조한다."""
    import formulas
    xl = formulas.ExcelModel().loads(path).finish()
    xl.calculate()
    import shutil
    outdir = os.path.join(HERE, "..", "out", "_xlsx_eval"); shutil.rmtree(outdir, ignore_errors=True); os.makedirs(outdir)
    xl.write(dirpath=outdir)
    evaluated = [f for f in os.listdir(outdir) if f.upper().startswith("SPHX_EQV_CALC")][0]
    wb = load_workbook(os.path.join(outdir, evaluated), data_only=True)
    ck = wb["CHECK"]
    rows = [(ck.cell(r, 1).value, ck.cell(r, 2).value, ck.cell(r, 3).value, ck.cell(r, 4).value)
            for r in range(1, ck.max_row + 1) if isinstance(ck.cell(r, 3).value, (int, float)) and ck.cell(r, 1).value]
    if not rows:
        print("평가 값 없음 — 검증 실패"); sys.exit(1)
    worst, bad = 0.0, 0
    for lab, x, ref, d in rows:
        ok = isinstance(x, (int, float)) and isinstance(d, (int, float))
        if not ok: bad += 1
        else: worst = max(worst, d)
        print(f"  {str(lab):14s} excel={x!s:>22} python={ref:14.6g} rel={d if ok else 'ERR'}")
    # 오류 셀 스캔
    errs = []
    for ws_ in wb.worksheets:
        for row in ws_.iter_rows():
            for c in row:
                if isinstance(c.value, str) and c.value.startswith("#"):
                    errs.append(f"{ws_.title}!{c.coordinate}={c.value}")
    print(f"formula error cells: {len(errs)}", errs[:20])
    print(f"max rel diff = {worst:.2e}", "OK" if worst < 1e-6 and bad == 0 and not errs else "MISMATCH")
    if worst >= 1e-6 or bad or errs:
        sys.exit(1)
    # 시나리오 2: 모드 B(α 입력) + 국부축 N_x=10 → 산 내측 σ_s/(N_x/h) = K_t = 17, max vM = 252.75 (docs/09 §3.3)
    wbf = load_workbook(path); ws = wbf["Input"]
    def setn(name, val):
        a = wbf.defined_names[name].attr_text.split("!")[1].replace("$", ""); ws[a] = val
    setn("in_mode", "B"); setn("in_alpha", 44.30459982548396); setn("in_Rc", None)
    setn("lc_axes", "local"); setn("lc_Nx", 10.0); setn("lc_My", 0.0); setn("lc_allow", 200.0)
    tmp = os.path.join(outdir, "scenario2.xlsx"); wbf.save(tmp)
    xl2 = formulas.ExcelModel().loads(tmp).finish(); xl2.calculate(); xl2.write(dirpath=os.path.join(outdir, "s2"))
    wv = load_workbook(os.path.join(outdir, "s2", "SCENARIO2.XLSX"), data_only=True)
    def getn(name):
        sh, a = wbf.defined_names[name].attr_text.split("!"); return wv[sh.strip("'").upper()][a.replace("$", "")].value
    kt, ss_in, vm, verdict, alpha_deg, tb = getn("Kt"), None, getn("vm_max"), getn("verdict"), getn("alpha_deg"), getn("t_b")
    st = wv["STRESS"]
    for r in range(1, st.max_row + 1):
        if st.cell(r, 1).value == "σ_s": ss_in = st.cell(r, 3).value
    print(f"scenario2: alpha={alpha_deg:.6f} t_b={tb:.6f} Kt={kt:.4f} crest inner sig_s/(Nx/h)={ss_in / (10 / 0.6):.4f} vm_max={vm:.4f} verdict={verdict}")
    assert abs(alpha_deg - 44.30459982548396) < 1e-9 and abs(tb - 3.7341051493604263) < 1e-9
    assert abs(ss_in / (10 / 0.6) - 17.0) < 1e-6 and abs(vm / 252.75037359510597 - 1) < 1e-4 and verdict == "NG"  # 파이썬은 800점 표본, 엑셀은 정확한 산 위치
    print("scenario2 OK")


if __name__ == "__main__":
    wb = build(); wb.save(OUT); print("saved", OUT)
    recalc_and_verify(OUT)
