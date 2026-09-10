# 10. 계산 모듈 (C-1) — YAML 입력 → 등가 판 보고서

**목적**: `calc/`의 개별 모듈(geometry → stiffness → transform → equivalent_plate → stress_recovery)을 하나의 실행 경로로 묶어, 04 양식의 치수를 YAML에 기입하면 등가 강성·등가 두께·영역 강성·ANSYS 입력·응력 복원 결과가 보고서로 나오게 한다. 엑셀 시트(C-2)와 파라메트릭 스터디(C-3)는 이 모듈의 결과를 기준값으로 삼는다.

## 1. 실행

```bash
python calc/plate_model.py calc/plate_input_template.yaml -o out/example   # 보고서 생성
python calc/plate_model.py --selfcheck                                      # 자체 검증
```

출력 디렉터리에 세 파일이 생긴다.

| 파일 | 내용 |
|---|---|
| `report.md` | 계산 보고서: 경고, 전처리된 입력, 형상 파생량, 국부축 채택식과 모델 비교, 등가 두께·균질 단일층 상수, ±β 영역 강성과 조합 추정, 단위 합력 응답, 하중 케이스 판정, 한계 |
| `results.json` | 같은 값의 기계 판독 형식 (C-2 엑셀 대조, C-3 파라메트릭 입력) |
| `sections.inp` | ANSYS preintegrated general shell section APDL 블록 (09 §2, 우·좌 영역 2절) |

가정 기본값에 대한 출력 예: `docs/example_report/`.

## 2. 입력 파일 (`calc/plate_input_template.yaml`)

04 양식과 05 정의서 항목을 그대로 따른다. 계산에 쓰는 키는 `section`, `chevron.beta_deg`, `material`, `load_cases`이고, 나머지(`plate`, `stack`, `loads`)는 기록용으로 보고서 경고에만 참조한다.

### 2.1 전처리 규칙 (표기 기준 → 계산 기준)

| 입력 플래그 | 규칙 | 근거 |
|---|---|---|
| `dims_reference: outer_surface` | $H = H_{in} - t$ | 05 §2.1 ($H_o = H + t$) |
| `pitch_reference: along_plate_axis` | $p = p_{in}\cos\beta$ | 능선 수직 방향은 $x_p$와 $\beta$를 이룬다 |
| `beta_reference: horizontal` | $\beta = 90° - \beta_{in}$ | 05 §2.2 |
| `material.t_min`, `corrosion_allowance` | $t_{calc} = t_{min} - CA$; 강성·응력 모두 $t_{calc}$ | 05 §2.5 |
| `material.rho` 미입력 | $7.9\times10^{-9}$ tonne/mm³ 가정 | SSPM 입력용 |

$R_c, R_v$는 중앙면 값으로 본다. 외면 반경으로 측정했으면 $R_o - t/2$로 환산해 기입한다.

### 2.2 하중 케이스

```yaml
load_cases:
  - name: "plate-axis My (right zone)"
    axes: plate          # local | plate
    zone: right          # right(θ=−β) | left(θ=+β) — plate 축일 때 국부축 변환에 사용
    N: [0, 0, 0]         # N_x, N_y, N_xy [N/mm]
    M: [0, 100, 0]       # M_x, M_y, M_xy [N·mm/mm], M_xy ↔ 공학 비틀림 κ_xy
    allowable: null      # MPa; null → material.S_allow
```

처리: 판축 합력 → `transform.resultants_to_local` → 채택식 역행렬로 거시 변형률 → `stress_recovery.ye_recovery` → 단면 위치·양면 응력 → 최대 von Mises, 능선 상세, 막 성분 최대값. `allowable`이 있으면 비와 OK/NG를 표시한다. **판축 합력 자체는 이 모듈이 산정하지 않는다** — 폐형식(등가 판 이론) 또는 등가 판 FE에서 얻어 입력한다.

## 3. 경고 규칙

| 조건 | 메시지 취지 |
|---|---|
| `t_min` 미입력 | 성형 두께 감소 미반영 |
| $t/R_{min} > 0.2$ | 쉘 이론 범위 밖, 능선 응력은 D단계 솔리드 FE 보정 (09 §5) |
| $H/p > 1$ | 문헌 검증 범위 밖 |
| $t_b/p > 0.3$ | 등가 판 요소 크기·횡전단 감도 확인 |
| `D_e`/`W_e` 미입력 | apex·테두리 보정, 전체 판 합력 산정은 범위 밖 |
| `T_design` 입력 & `E_source` 없음 | E가 설계온도 값인지 확인 |

## 4. 자체 검증 (`--selfcheck`)

- 기본값 YAML에서 $A_{11}=2603$, $D_{22}=170\,110$ N·mm, $t_b=3.734$ mm (06 §5, 09 §1 재현).
- 전처리 4규칙: outer $H$=3.8→3.2, horizontal 30°→60°, along-axis $p$=18→9, $t_{min}$ 0.55−0.05→0.5 (형상 객체 두께까지 전달).
- 하중 케이스: $\beta=0$에서 plate 축 = local 축 결과 동일; 단위 응답 × 10 = 케이스 결과; 능선 내측 $\sigma_s/(N_x/h)=K_t$; $\beta=60°$ plate 축 $N_{x_p}$만 부과 시 국부 $N_x, N_y, N_{xy}$ 모두 비영; `allowable` 판정.
- 세 출력 파일 생성 확인.

## 5. 가정 기본값 결과 요지 (`docs/example_report/report.md`)

- 채택식 국부축: $A_{11}=2603$, $A_{22}=146\,400$, $A_{66}=35\,287$ N/mm; $D_{11}=3025$, $D_{22}=170\,110$, $D_{66}=1687$ N·mm.
- $t_b=3.73$, $t_s=t_m=0.757$ mm. 균질 단일층은 $G$ 선택에 따라 $D_{66}$ +2331 % 또는 $A_{66}$ −96 %.
- 단위 응답(국부축, 최대 von Mises/단위): $N_x$ 25.3, $N_y$ 1.32, $N_{xy}$ 3.17 MPa per N/mm; $M_x$ 15.0, $M_y$ 2.16, $M_{xy}$ 23.6 MPa per N·mm/mm.
- 판축 $M_{y_p}=100$ N·mm/mm을 우측 영역 국부축으로 옮기면 $M_x=75$, $M_y=25$, $M_{xy}=-43.3$이 되어 최대 1509 MPa — 회전으로 생기는 $M_x$, $M_{xy}$ 성분이 지배한다. 이 값은 "영역이 판축 모멘트를 그대로 받는다"는 가정이며, 실제 영역별 모멘트 분배는 zone-wise 등가 판 FE(E단계)에서 얻는다.

## 7. 엑셀 계산 시트 (C-2, `calc/SPHX_EQV_calc.xlsx`)

설계 현장용으로 같은 로직을 엑셀 수식으로 구현했다. 생성기는 `calc/make_excel_sheet.py`이며, 시트를 고칠 때는 생성기를 수정하고 재생성한다.

| 시트 | 내용 |
|---|---|
| README | 사용법, 색 규칙(파랑+노랑 = 입력, 검정 = 수식, 녹색 = 시트 간 참조), 범위·한계 |
| Input | 단면(p, H, t, R_c, R_v 또는 α), 기준 플래그(dims/pitch/beta reference), β, 재료(E, ν, ρ, t_min, CA, S_allow), 하중 케이스 1건(axes, zone, N, M, allowable). 드롭다운 검증 포함 |
| Geometry | 전처리 4규칙, 원호+접선 해(모드 A: Newton 20회 반복표 — $F(\alpha)=R_s(1-\cos\alpha)+(c-R_s\sin\alpha)\tan\alpha-H$, $F'=(c-R_s\sin\alpha)/\cos^2\alpha$, 반복값을 $(0.001, \alpha_{hi})$로 제한; 모드 B: α 입력 → R 선형해), $l, I_1, I_2, \lambda, J_1, J_2, A_u, I_u$, 무차원, 경고 4행 |
| Stiffness | 채택식 $A, D$ 8성분, 직교이방성 비, 모델 비교(Xia $D_{22}$, Lang & Su $A_{11}, A_{22}$, Ye 얕은식), 등가 두께 $t_b, t_s, t_m$, 균질 단일층 상수와 66 잔차, $K_t$ |
| Rotation | 우·좌 영역 $\bar A,\bar D$ 12성분(Q̄ 폐형식), 횡전단 추정, 면적 밀도, **SSPA/SSPB/SSPD/SSPE/SSPM 문자열**(TEXT 함수로 조립, 복사용) |
| Stress | 판축 합력 → 국부축($T_\sigma$) → 거시 변형률(채택식 역행렬) → VAM 상수($\mathcal C, c_1, c_2, c_4, \alpha_1$) → 산·골 물리 변형률 → 4 표면점 응력·von Mises → 최대값·위치, 능선 합력, 막 von Mises, 허용응력 비·판정 |
| Check | 파이썬 `results.json` 참조값(하드코딩) 대비 상대오차 28항목과 최대값 — 기본값 입력에서만 의미 |

**파이썬 모듈과의 범위 차이**: 단면은 원호+접선만(사인형·사다리꼴·반원은 파이썬), 응력은 산·골 두 위치만(플랭크 분포는 파이썬). 6개 단위 하중 모두 최대 응력이 산 또는 골에 생기므로(§5) 두께 사이징 용도에는 충분하다.

**검증**: 이 환경에는 LibreOffice Calc가 없어 `formulas` 라이브러리(pip)로 워크북 전체를 평가했다. 기본값에서 Check 28항목 최대 상대오차 $1.5\times10^{-7}$(파이썬의 800점 표본 위치 차이), 시나리오 2(모드 B α 입력 + 국부 $N_x=10$, 허용 200 MPa)에서 $K_t=17.000$, 산 내측 $\sigma_s/(N_x/h)=17.000$, 최대 von Mises 252.75 MPa, 판정 NG를 재현. 수식 오류 셀 0. 엑셀에서는 열 때 전체 재계산된다(`fullCalcOnLoad`).

**정의 이름 규칙(생성기)**: 셀 주소와 같은 이름(A11, D66, I1, C1, Z1)과 대소문자만 다른 이름(H_, h_)은 엑셀에서 쓸 수 없어 `A_11`, `H_mid`, `h_calc` 등으로 치환하고, 생성 수식의 셀 참조는 모두 절대형(`$D$22`)으로 써서 이름 치환과 충돌하지 않게 했다.

## 8. 다음 단계

- **C-2** 엑셀 시트: 완료 (§7). 실제 판 치수 수령 후 Input 갱신, Check 참조값은 파이썬 재실행 후 생성기로 갱신.
- **C-3** 파라메트릭: `PlateModel`을 $\beta$, $H/p$, $t$, $R/t$ 격자에 대해 호출, $K_t$·$A_{22}/A_{11}$·단위 응답 경향.
- **D/E**: `sections.inp` 를 단일 요소 시험(비틀림 규약)과 전체 판 등가 FE에 투입.
