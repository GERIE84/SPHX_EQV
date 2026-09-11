# SPHX_EQV — 쉐브론 Shell & Plate 열교환기 등가모델 개발

쉐브론(chevron) 주름형상 박판으로 구성된 Shell & Plate 열교환기(SPHX)의 전열판에 대해, 판 두께의 강도를 **수식으로 평가할 수 있는 등가모델**을 개발하고, 이어서 도출된 **등가 물성치를 FEA에 적용**하는 것을 목표로 하는 프로젝트이다.

## 목표

| 단계 | 내용 | 산출물 |
|---|---|---|
| 1 | 쉐브론 형상 치수(피치, 진폭, 각도, 두께 등)를 반영한 등가 수식 모델을 개발하여, 판 두께의 강도를 수식/엑셀로 사전 평가한다. | 등가 강성·등가 응력 수식, 계산 시트 |
| 2 | 1단계에서 도출한 등가 물성치(직교이방성 E, G, ν 등)를 FEA에 적용하여, 대형 조립체 해석에서 주름판을 등가 평판으로 치환한다. | 등가 물성치 산정 절차, FE 모델 검증 결과 |

## 현재까지의 결정사항

- **용도 우선순위**: 1순위 사전 강도계산(수식/엑셀) → 2순위 대형 조립체 FE에서의 판 치환.
- **대상 범위**: 단일 판을 먼저 다루고, 교차 적층 스택(판 간 접촉 포함)은 확장 단계로 미룬다.
- **진행 방식**: 문헌조사가 어느 정도 정리된 뒤 이후 작업(수식 구현, 해석 자동화)은 Claude Code로 이관한다.
- **FEA 미사용 (2026-09-11)**: 당분간 유한요소해석을 쓰지 않는다. D·E 단계는 보류하고, 압력차로부터 합력·응력을 구하는 폐형식(C-4, `docs/12`)으로 1순위 목표를 FEA 없이 완결한다. 채택식 근거는 문헌 FE 벤치마크로 갈음한다.
- **기준 문헌**: Lang & Su, *Equivalent orthotropic model for corrugated plates based on simplified constitutive relation*, Heliyon (2022) — 등가 직교이방성 모델의 출발점으로 삼는다.

## 저장소 구성

```
README.md                 프로젝트 개요 및 결정사항
docs/
  ref/                    참고 문헌 원문 (오픈액세스)
  01_plan.md              개발 계획 및 로드맵
  02_literature.md        문헌조사 (등가 직교이방성 모델 중심)
  03_next_steps.md        Claude Code 이관 이후 작업 항목
  04_plate_dimension_form.md  전열판 치수 정의 그림 및 입력 양식
  05_geometry_parameters.md   형상 파라미터 정의서 (기호·측정기준·파생량·유효범위)
  06_equivalent_stiffness.md  등가 강성 수식 정리 및 채택식 (국부축, 모델 비교, 수치 예)
  07_stiffness_transformation.md  β 회전 변환 및 ±β 영역 조합 규칙
  08_method_review_rve_vs_vam.md  방법론 검토: RVE vs VAM, FEA 필요 범위, 응력 복원식 채택
  09_equivalent_thickness.md  등가 두께 정의, FE 입력 방식(GENS) 결정, VAM 응력 복원 구현·검증
  10_calc_module.md       계산 모듈 사용법 (YAML 입력 → 보고서), 전처리·경고 규칙, 하중 케이스, 엑셀 시트(§7)
  11_parametric_study.md  파라메트릭 스터디: H/p, t/p, R_c/t, β 가 강성·등가 두께·K_t·단위 응답에 미치는 영향
  12_pressure_closed_form.md  압력차 폐형식 (C-4): 접촉점 격자, 주름 단면 곡선 프레임, 능선 연속보, Hertz 접촉, 8점 중첩
  param_study/            파라메트릭 케이스별 결과 CSV
  example_report/         가정 기본값 입력에 대한 출력 예 (report.md, results.json, sections.inp)
  fig/                    그림 및 생성 스크립트
calc/                     계산 모듈: plate_model.py(YAML 입력 → 보고서·JSON·APDL, C-1 진입점), SPHX_EQV_calc.xlsx(현장용 엑셀 시트, 생성기 make_excel_sheet.py), parametric.py(C-3 파라메트릭), pressure_resultants.py(C-4 압력차 폐형식), geometry.py(단면 라이브러리), stiffness.py(등가 강성 5모델·채택식), stress_recovery.py(VAM 응력 복원 + Briassoulis 능선식), transform.py(β 회전·영역 조합), equivalent_plate.py(등가 두께·ANSYS GENS APDL 생성), plate_input_template.yaml(입력 템플릿),
                          verify_{langsu,ye,xia,briassoulis}_tables.py(문헌 수치 재현 검증)
fea/                      (예정) 등가 물성치 FE 검증 모델
```

## 상태

- 2026-09: 프로젝트 정의, 용도·범위 결정, 문헌조사 착수.
- 2026-09-08: 기준 문헌 Lang & Su (2022) 원문 확보, 등가 강성식·형상 파라미터 추출 완료 (`docs/02_literature.md`), 수식 전사 검증 스크립트 추가 (`calc/`).
- 2026-09-08: Ye et al. (2014) VAM 등가판 모델 전사·검증, 고전식(Seydel·Briassoulis) 통합 비교표 및 채택식 권고 작성. 국부 변형률 복원식 확보로 응력 평가 접근 방향 갱신.
- 2026-09-09: B-1 형상 파라미터 정의서(`docs/05`) 및 단면 라이브러리(`calc/geometry.py`) 작성. 실제 판 치수 입력 대기.
- 2026-09-09: 방법론 검토(`docs/08`) — 강성은 이미 VAM 채택식, 응력 복원도 VAM 복원식을 1차 수단으로 결정. FEA는 검증·비주기 영역용으로 범위 불변.
- 2026-09-11: FEA 미사용 결정 반영(D·E 보류). C-4 압력차 폐형식(`docs/12`, `calc/pressure_resultants.py`, 엑셀 Pressure 시트): 1피치 스팬은 균질화 판이 아닌 단면 곡선 프레임으로 풀어야 함을 확인(균질화 시 2~3배 과소). 기본값 판 156.6 MPa per MPa ΔP.
- 2026-09-10: C-3 파라메트릭 스터디(`docs/11`). 등가 굽힘 두께는 판 두께와 무관한 형상량, 가로 인장·비틀림 응력 ∝ t⁻², 반경은 2차 인자, β 전 범위에서 굽힘–비틀림 연성 0.9 수준 확인. C단계 종료.
- 2026-09-10: C-2 엑셀 계산 시트(`calc/SPHX_EQV_calc.xlsx`). 파이썬과 동일 로직(원호+접선, 채택식, 등가 두께, 영역 강성, APDL 문자열, 산·골 응력·판정), 수식 평가로 파이썬 대비 일치 확인.
- 2026-09-10: C-1 계산 모듈 래핑(`calc/plate_model.py`, `docs/10`). YAML 한 파일로 형상→강성→등가 두께→영역 강성→ANSYS 절→응력 복원·판정까지 보고서 생성. 실제 판 치수 입력 대기.
- 2026-09-09: B-4 등가 두께·FE 입력 방식·VAM 응력 복원(`docs/09`, `calc/equivalent_plate.py`, `calc/stress_recovery.ye_recovery`). 균질 단일층 불가 확인 → ANSYS preintegrated section 채택; 복원식은 평판 극한·Case 1/5·에너지 일치로 검증. B단계 종료.
- 2026-09-09: B-3 β 좌표 변환·영역 조합 규칙 정리(`docs/07`, `calc/transform.py`). 영역별 연성항의 크기와 조합 방식 민감도 확인, FE는 zone-wise 채택.
- 2026-09-09: Briassoulis (1986) 원문 확인 — $D_{66}$ 관례 문제 해결, 독립 FE 벤치마크로 채택식 재검증, 능선 응력집중식 $K_t=1+6f/t$ 확보.
- 2026-09-09: B-2 등가 강성 수식 정리(`docs/06`, `calc/stiffness.py`). 채택식 확정(Ye 2014 Eq. 19 기반), 모델 간 차이 정량화, 면내/굽힘 등가두께 불일치 확인.
- 2026-09-08: Xia et al. (2012) 원문 전사·검증. 3대 기준 문헌 확보 완료, 문헌조사 A-1·A-2 항목 종결. 원호 주름 폐형식·직교이방성 원판 일반형·케이스별 국부 내력 가정 확보.
