# 03. 다음 작업 항목

문헌조사 이후 Claude Code로 이관하여 진행할 작업을 정리한다.

## A. 문헌조사 마무리
- [x] Lang & Su (2022) 원문에서 등가 굽힘강성 수식과 단면 파라미터 정의를 추출하여 `docs/02_literature.md`에 기입 (2026-09-08. 원문 `docs/ref/`, 수식 검증 `calc/verify_langsu_tables.py`)
- [x] Briassoulis, Samanta & Mukhopadhyay, Xia et al. 등가 강성식 비교표 작성 (2026-09-08 완료 — `02_literature.md` §3, §4)
  - [x] Xia (2012) / Lang & Su (2022) / Ye (2014) 수식 비교표 및 수치표(Table 5·6) — `02_literature.md` §1.5, §1.7
  - [ ] Briassoulis (1986), Samanta & Mukhopadhyay (1999) 원문 수식 추가 — **원문 PDF 필요**
  - [x] Ye (2014) 원문 확보·전사·수치검증 (2026-09-08, `02_literature.md` §2, `calc/verify_ye_tables.py`)
  - [x] Seydel / 1960년대 통용식 / Briassoulis 수식 — Ye (2014) 2차 전사로 확보·검증 (`02_literature.md` §3)
  - [x] Xia (2012) 원문 확보·전사·수치검증 (2026-09-08, `02_literature.md` §3, `calc/verify_xia_tables.py`). $A_{22}$ 원식 확정, Ye 전사 오류 확인
  - [x] Samanta & Mukhopadhyay, Yokozeki 수식 — Xia (2012) Table 3·4 전사로 확보 (§3.4)
  - [x] Briassoulis (1986) 원문 확보·확인 (2026-09-09, `02_literature.md` §4.1, `calc/verify_briassoulis_tables.py`) — $D_{66}$ 관례 문제 해결, 독립 FE 벤치마크, 응력집중식 확보
  - [ ] (선택) Samanta & Mukhopadhyay 원문 — 계수 불일치 확인. 우선순위 낮음
- [ ] 쉐브론 판 강성/좌굴 관련 문헌 확보

## B. 형상 파라미터 및 수식 정리
- [x] 쉐브론 형상 파라미터 정의서 작성 (p, h, β, t, 단면 형상, 유효 판 크기) — 2026-09-09 `docs/05_geometry_parameters.md`
  - [x] 치수 정의 그림 및 입력 양식 작성 (2026-09-08, `docs/04_plate_dimension_form.md`, `docs/fig/plate_dimensions.png`)
  - [x] 단면 라이브러리 `calc/geometry.py` (사인형/원호+접선/사다리꼴/원호, 파생량 $c,l,I_1,I_2$, 문헌값 자체검증) 및 입력 템플릿 `calc/plate_input_template.yaml`
  - [ ] 실제 판 치수·재질·설계조건 값 수령 → 템플릿 기입 (그 전까지 05 §5 가정값 사용)
- [x] 주름 방향/직각 방향 등가 면내 강성(E_x, E_y, G_xy, ν_xy) 및 굽힘 강성(D_x, D_y, D_xy) 수식 정리 — 2026-09-09 `docs/06_equivalent_stiffness.md`, `calc/stiffness.py` (채택식 확정: Ye Eq. 19 기반, 문헌 표·평판 극한 검증)
  - [ ] 비대칭 단면($R_c\ne R_v$) 확인 시 연성 강성 $B_{ij}$ 구현 (06 §6)
- [x] 쉐브론 각도 β에 대한 좌표 변환 및 조합식 정리 — 2026-09-09 `docs/07_stiffness_transformation.md`, `calc/transform.py` (회전식 검증, zone-wise/Strip/Voigt/Reuss 조합 규칙, 채택: FE는 zone-wise, 수식은 겉보기 강성 + Strip/zone 범위)
  - [ ] 비대칭 단면 시 6×6 (A,B,D) 회전으로 확장
- [x] 등가 두께 정의(면내 기준/굽힘 기준) 및 응력 환산 계수 개념 정리 — 2026-09-09 `docs/09_equivalent_thickness.md`, `calc/equivalent_plate.py` ($t_b=\sqrt{12D_{11}/A_{11}}$=3.73 mm, $t_s=h\lambda$=0.76 mm; 균질 단일층 불가 → FE 입력은 ANSYS preintegrated general shell section(GENS) 채택, APDL 생성기 포함)
  - [x] Briassoulis (1986) 능선 응력집중식(B1–B7) 구현 — `calc/stress_recovery.py` (2026-09-09)
  - [x] 방법론 검토: RVE vs VAM, FEA 필요 범위 — `docs/08_method_review_rve_vs_vam.md` (2026-09-09). 결정: 응력 복원 1차 수단 = VAM 복원식(Ye Eq. 14–17, 20); Briassoulis·Xia Case 1은 교차검증용; Xia Case 5 가정 폐기
  - [x] `stress_recovery.ye_recovery()` 구현 — 6개 거시 변형률 → 단면 위치별 국부 막·굽힘 변형률·응력·합력, `recovery_from_resultants()` (2026-09-09)
  - [x] 구현 검증: 평판 극한, Case 1·5 특수해(08 §4), 단위셀 변형에너지 일치(9 조합), Briassoulis Fig. 6 경향, 왕복 — 전부 통과 (09 §3.2)
  - [x] $t/R_{min}>0.2$ 시 보정 방침 — 강성 유지, 능선 응력은 D단계 솔리드 FE 보정계수 (09 §5)
  - [ ] GENS 절 비틀림($D_{33}$, 16/26) 규약 단일 요소 시험으로 확인 → E-1
  - [ ] 비대칭 단면($B_{ij}\ne0$) 복원식 $\mathcal B,\alpha_2$ 항 추가 (06 §6과 함께)

## C. 계산 시트 구현 (`calc/`)
- [x] 파이썬 모듈로 등가 물성치 계산 함수 구현 (입력: 형상 파라미터, 재료; 출력: 등가 물성치, 등가 두께) — 2026-09-10 `calc/plate_model.py` (YAML → report.md / results.json / sections.inp; 전처리 규칙, 경고, 하중 케이스 응력 복원·판정), `docs/10_calc_module.md`, 예시 출력 `docs/example_report/`
  - [ ] 실제 판 치수 수령 후 템플릿 기입·실행 (경고 항목 해소: t_min, D_e/W_e, S_allow)
- [ ] 엑셀 계산 시트 동일 로직 구현 (설계 현장 사용용) — `results.json` 대조
- [ ] 파라메트릭 스터디 (β, h/p, t 변화에 따른 등가 강성 경향)
  - [ ] $K_t = 1+6f/t$ 및 $A_{22}/A_{11}$ 비의 $f/t$, $H/p$ 의존성 포함 (두께 결정 지배 인자 파악)

## D. 단일 판 상세 FE 검증 (`fea/`)
- [ ] 쉐브론 단위 셀/단일 판 상세 쉘 모델 구축 (ANSYS APDL, 파라메트릭)
  - [ ] 일방향 주름 단위셀(Dirichlet BC, Xia Table 1 방식)로 등가 강성 산출 루틴 검증 → 쉐브론 단위셀로 확장
  - [ ] 전체 판 굽힘 해석으로 $D_{22}$ 상한 문제(Xia vs Lang & Su/Ye) 재현 확인
  - [ ] 비틀림 케이스로 $D_{66}$의 $\lambda$ 계수 확정 (Ye/VAPAS vs Briassoulis FE, 06 §7-7)
- [ ] 면내 인장/전단, 굽힘, 압력차 하중 케이스 해석
  - [ ] 단위셀 6개 일정 변형률·곡률 케이스(Xia Table 1 / Briassoulis §4 방식) — 순수 비틀림 $\kappa_{xy}$ 케이스 필수
  - [ ] Briassoulis Table 1·2(원호+접선 얕은 주름)를 첫 검증 케이스로 재현해 FE 절차 자체를 검증
- [ ] 등가 강성 및 최대응력 추출 → 수식 결과와 비교, 환산 계수 산정
  - [ ] Briassoulis $K_t$·Ye 복원식으로 계산한 국부 응력과 상세 FE 응력 비교 (일방향 주름 단위셀 → 쉐브론 꺾임부 보정)
  - [ ] $t/R_{min}>0.2$ 형상은 솔리드 요소로 두께 방향 응력 분포 확인 (일방향 주름 단위 셀 → 쉐브론 꺾임부 보정)

## E. 등가 물성치 FE 적용
- [ ] 직교이방성 쉘 등가 평판 모델 구축 및 상세 모델과 변위·반력 비교
  - [ ] zone-wise 이방성(16/26 포함) 입력 방식 확정 (preintegrated section vs 회전 층상 쉘)
  - [ ] apex 경계의 전단·비틀림 구속 정도 확인 (Strip vs zone-wise 사이 실제 위치, 07 §4.2)
- [ ] 조립체 모델 적용 절차 정리

## F. 확장 단계
- [ ] 교차 적층 스택(접촉 포함) 등가 강성 검토
