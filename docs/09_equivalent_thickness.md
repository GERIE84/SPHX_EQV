# 09. 등가 두께 정의, FE 입력 방식 결정, VAM 응력 복원 구현 (B-4)

**목적**: (1) 국부축 채택식 $A_{ij}, D_{ij}$(`06`)로부터 "등가 두께"를 어떻게 정의하고 어디에 쓰는지 확정한다. (2) 2단계 조립체 FE에 등가 판을 넣는 입력 방식을 결정한다. (3) `08` 결정에 따라 VAM 응력 복원식 `calc/stress_recovery.ye_recovery()`를 구현·검증하고, 합력 → 최대응력 절차를 확정한다.

구현: `calc/equivalent_plate.py`(등가 두께, 균질 단일층 공학상수, ANSYS preintegrated section APDL 생성), `calc/stress_recovery.py`(Briassoulis 폐형식 + VAM 복원). 수치 예는 모두 05 §5 가정 기본값(원호+접선 $p=9$, $H=3.2$, $t=0.6$, $R_c=R_v=1.5$ mm, 316L $E=193$ GPa, $\nu=0.3$)이다.

---

## 1. 등가 두께 — 정의와 역할

등가 직교이방성 판은 $A$(6개)와 $D$(6개)로 완전히 정의되며 "두께"는 이론상 불필요하다. 두께가 필요한 곳은 세 가지이고, 각각 다른 값이 나온다.

| 기호 | 정의 | 용도 | 기본값 |
|---|---|---|---|
| $t_b$ | $\sqrt{12D_{11}/A_{11}}$ (= $\sqrt{12D_{22}/A_{22}}$ = $\sqrt{12D_{12}/A_{12}}$, 채택식에서 세 값이 동일) | 균질 단일층으로 면내·굽힘(11, 12, 22)을 동시에 재현할 때의 두께. 판 처짐·좌굴 판정용 "판 두께" 감각치 | **3.734 mm** |
| $t_s$ | $\sqrt{12D_{66}/A_{66}} = h\lambda$ | 면내 전단·비틀림을 동시에 재현하는 두께. $t_b$와 다르므로 균질 단일층 불가의 근거 | 0.757 mm |
| $t_m$ | $h\lambda = h\,l/c$ | 단위 투영면적당 질량·전개 면적(무게, 열용량, 밀도 입력) | 0.757 mm |
| $h$ | 원판 두께 | 응력 복원(쉘 구성식), 부식 여유·공차 반영은 05 §2.5 | 0.600 mm |

채택식에서 $D_{11}/A_{11}=D_{12}/A_{12}=D_{22}/A_{22}=h^2(12J_2+h^2J_1)/(12\lambda h^2)$ … 임을 `equivalent_plate._selfcheck`가 확인한다(세 값의 편차 $<10^{-9}$). 따라서 국부축에서는 **$t_b$ 하나 + 공학상수 4개 + 전단 1개**로 11/12/22 성분을 정확히 맞출 수 있고, 남는 불일치는 66 성분 하나다.

### 1.1 균질 직교이방성 단일층으로 넣을 때 (`homogeneous_orthotropic`)

$t_e=t_b$, $E_x=(A_{11}A_{22}-A_{12}^2)/(A_{22}t_e)$, $E_y=(\cdot)/(A_{11}t_e)$, $\nu_{xy}=A_{12}/A_{22}$, 전단은 두 선택지:

| 항목 | 값 | $A_{66}$ 잔차 | $D_{66}$ 잔차 |
|---|---|---|---|
| $E_x$ | 696 MPa | | |
| $E_y$ | 39 142 MPa | | |
| $\nu_{xy}$ | 0.0053 | | |
| $G_{xy}$ (A 기준 $A_{66}/t_e$) | 9 450 MPa | 0 | **+2 331 %** |
| $G_{xy}$ (D 기준 $12D_{66}/t_e^3$) | 389 MPa | **−95.9 %** | 0 |

→ 국부축에서도 면내 전단과 비틀림 중 하나는 24배 틀린다. 쉐브론 판축으로 $\mp\beta$ 회전하면(07) 상황이 더 나빠진다: $\beta=60°$에서 $\sqrt{12\bar D_{ii}/\bar A_{ii}}$가 11: 3.27, 22: 2.12, 66: 3.27, 16: 4.51, 26: 2.69 mm로 모두 달라 **어떤 두께로도 균질 단일층이 성립하지 않는다.** 균질 단일층은 국부축 간이 검토(예: 1방향 주름 스트립의 처짐 확인)에만 쓴다.

## 2. FE 입력 방식 결정

| 방식 | 재현 범위 | 판정 |
|---|---|---|
| (i) **preintegrated general shell section** — ANSYS `SECTYPE,,GENS` + `SSPA/SSPB/SSPD/SSPE/SSPM`, SHELL181/281 | $A, B, D$ 12+6개 성분을 그대로 입력. 16/26 연성항 포함 | **채택** (1차) |
| (ii) 균질 직교이방성 단일층 ($t_b$ + 공학상수) | 11/12/22만 (국부축) | 간이 검토·비교용 |
| (iii) 3층 샌드위치(강성 외피 + 연성 코어)로 $A_{66}, D_{66}$ 맞춤 | 국부축 6성분은 맞출 수 있으나 회전 후 16/26은 층 배향으로만 가능 | GENS 미지원 코드용 대안, 필요 시 |

**(i) 입력 절차** (`equivalent_plate.chevron_gens_deck`)

1. 국부축 채택식 $A, D$ → `transform.zone_stiffness`로 우측 영역($\theta=-\beta$)·좌측 영역($\theta=+\beta$) 판축 강성 $\bar A,\bar D$(16/26 포함) 생성.
2. 영역별 절(section) 2개: `SSPA,A11,A21,A31,A22,A32,A33`, `SSPD,D11,…,D33` — 인덱스 3 = xy. 판축 $\bar A_{16}\to$ A31, $\bar A_{26}\to$ A32, $D$ 동일. `SSPB` = 0(대칭 단면; 비대칭 단면이면 06 §6의 $B_{ij}$ 입력).
3. 요소 좌표계는 두 영역 모두 판축 $(x_p, y_p)$에 맞춘다(ESYS 회전 불필요). apex 선은 요소 경계와 일치(07 §5).
4. 대안 `mode='local_axes'`: 국부축 직교이방 $A, D$ 한 절 + 영역별 ESYS를 $\mp\beta$ 회전. 성분 수가 적어 입력 검증이 쉽다. 두 방식의 결과 일치가 E-1 첫 검사 항목.
5. 횡전단 `SSPE,E11,E21,E22`: Kirchhoff 등가모델에는 없는 량. 1차 추정 $E_{11}=E_{22}=\tfrac56\,G\,h\lambda$ (기본값 46 850 N/mm), $E_{21}=0$. 판 두께 대비 스팬이 크므로 결과에 영향이 없어야 하며 E단계에서 ×0.1/×10 감도로 확인한다.
6. 밀도 `SSPM,DENS`: 단위 두께 가정 밀도 = 단위 투영면적당 질량 $\rho h\lambda$ (기본값 $7.9\times10^{-9}\times0.757=5.98\times10^{-9}$ tonne/mm²).
7. 단면 오프셋: 등가 판 중립면 = 주름 중앙면($z=0$). 판 사이 간격(형상 04 양식 $H$, 판 간 접촉) 모델링 시 노드를 중앙면에 두면 오프셋 불필요.

**명령 인자 순서 근거**: ANSYS 매뉴얼 미러가 차단되어 `ansys-mapdl-core 0.74.1` 명령 레퍼런스 docstring(PyPI 배포본)에서 확인했다(2026-09-09). SSPA/SSPD 순서 `11,21,31,22,32,33`, SSPE `11,21,22`, SSPM `DENS[,T]`(단위 두께 가정), `SECTYPE,secid,GENS` 뒤에 Subtype/REFINEKEY 없음.

**확인 필요(E-1)**: GENS 일반화 변형률의 비틀림 성분이 공학 비틀림 $\kappa_{xy}$(=$2\partial^2 w/\partial x\partial y$, 본 프로젝트 $D$ 관례)인지 텐서 성분인지 매뉴얼 원문으로 미확정. 단일 요소 순수 비틀림 시험(모서리 변위 부과 → $M_{xy}$ 반력)으로 $D_{33}$ 규약을 확인한 뒤 필요하면 33 항에 계수 2 또는 1/2를 적용한다. 회전 후 16/26 항은 같은 규약을 따르므로 함께 확인된다.

**요소 조건수**: 국부축 $A_{22}/A_{11}=56$, 판축에서는 최대/최소 대각 비 ≈ 8로 완화된다(07 §4.1). GENS 입력은 층상 쉘의 두께 적분을 거치지 않으므로 국부축 비율 자체는 문제 없고, 얇은 판 스팬 비($t_b$/스팬)가 조건수를 지배한다.

## 3. VAM 응력 복원 구현 (`stress_recovery.ye_recovery`)

### 3.1 식 (대칭 단면, Ye Eq. 14–17, 20; 프로젝트 기호)

거시 변형률 $\mathbf m=\{\epsilon_{xx},\epsilon_{yy},\gamma_{xy},\kappa_{xx},\kappa_{yy},\kappa_{xy}\}$ ($\gamma_{xy}=2\epsilon_{xy}$, $\kappa_{\alpha\beta}=w_{,\alpha\beta}$). 단면 위치 $x$에서 $a=1+z'^2$, $\sqrt a=ds/dx$, $z''$, 쉘 곡률 $\kappa=z''/a^{3/2}$, 두께 보정 $c_h=1+h^2\kappa^2/48$.

$$
\mathcal C=-\Big(\frac{12J_2}{h^2}+J_1\Big),\quad
c_1=-\frac{\epsilon_{xx}+\nu\epsilon_{yy}}{\mathcal C},\quad
c_4=\frac{\kappa_{xx}+\nu\kappa_{yy}}{\lambda},\quad
\alpha_1=\Big\langle\frac{\sqrt a}{c_h}\Big\rangle^{-1},\quad c_2=\alpha_1\gamma_{xy}
$$
$$
\begin{aligned}
\gamma^0_{11}&=c_1\sqrt a-\nu a(\epsilon_{yy}+z\kappa_{yy}) &
\rho^0_{11}&=a\Big(\frac{12c_1z}{h^2}+c_4\Big)-\nu\sqrt a\,\kappa_{yy}\\
\gamma^0_{22}&=\epsilon_{yy}+z\kappa_{yy} &
\rho^0_{22}&=\kappa_{yy}/\sqrt a\\
2\gamma^0_{12}&=\frac{1}{c_h}\Big(\sqrt a\,c_2-\frac{h^2z''\kappa_{xy}}{12a}\Big) &
2\rho^0_{12}&=\frac{1}{c_h}\Big(-2\sqrt a\,\kappa_{xy}+\frac{z''c_2}{2a}\Big)
\end{aligned}
$$

물리 성분: $\varepsilon_s=\gamma^0_{11}/a$, $\varepsilon_y=\gamma^0_{22}$, $\gamma_{sy}=2\gamma^0_{12}/\sqrt a$, $\kappa_s=\rho^0_{11}/a$, $\kappa_y=\rho^0_{22}$, $2\kappa_{sy}=2\rho^0_{12}/\sqrt a$. 표면 변형률 $e(\zeta)=\varepsilon-\zeta\kappa$, $\zeta=\pm h/2$(법선 $n$은 산에서 $+z$; 외측 = $+h/2$). 평면응력 $\sigma_s,\sigma_y,\tau_{sy}$ → von Mises. 국부 합력 $N_s, N_y, M_s, M_y$도 함께 반환.

**부호 확정**: Ye 인쇄본 Eq. (14)의 $\rho^0_{11}$ 마지막 항과 $\rho^0_{22}$ 부호는 OCR·전사 불확실(02 §2.6). 본 구현은 (a) 평판 극한에서 $\kappa_s\to\kappa_{xx}$, $\kappa_y\to\kappa_{yy}$로 환원되고, (b) 복원 변형률로 계산한 단위셀 에너지가 $\tfrac12\mathbf m^T[A,D]\mathbf m$과 일치하는 부호 조합으로 확정했다(§3.2 [A], [D]). 잘못된 부호에서는 [D]의 `kxx+kyy` 교차항이 실패한다.

### 3.2 검증 결과 (`python calc/stress_recovery.py`, 전부 통과)

| 검사 | 내용 | 결과 |
|---|---|---|
| [A] 평판 극한 | $H\to0$ 사인형에서 6개 복원 변형률 = 입력 거시 변형률 | 상대오차 $<10^{-13}$ |
| [B] Case 1 ($\epsilon_{xx}$) | 능선 $N_s=A_{11}\epsilon_{xx}$, $M_s=N_xf$, 내측 $\sigma_s/(N_x/h)=K_t=1+6f/h=17.0$ (Xia Case 1, Briassoulis B2) | 1e-5 |
| [C] Case 5 ($\kappa_{yy}$) | 능선 $N_y=Ehf\kappa_{yy}=185.3$ N/mm (Briassoulis Eq. 12); Xia 가정 203.6은 불일치 | 1e-5 |
| [D] 에너지 | 6개 단일 성분 + 3개 조합(`exx+eyy`, `kxx+kyy`, 6성분 동시)에서 셀 에너지 = 거시 에너지(`stiffness.ye_full`) | $\le5\times10^{-6}$ |
| [E] Briassoulis Fig. 6 경향 | 순수 $\kappa_{yy}$에서 능선 $|\sigma_y|_{max}/(6M_y^*/t^2)$: $f/t$=0.25→1.18, 0.5→0.97, 1→0.54, 2→0.22, 4→0.07 | 두꺼운 판에서 >1, 얇은 판에서 급감 — Fig. 6 경향 재현 |
| [F] 왕복 | 합력 → 채택식 역행렬 → 복원 → 최대 von Mises가 직접 입력과 동일 | 0 |

[B], [C]의 1e-5는 수치 적분 표본(800점)에 따른 것이고 표본 수를 늘리면 감소한다.

### 3.3 기본값 판 단위 합력 응답 (국부축, 1 단위씩; 최대 von Mises와 위치)

| 합력 | max σ_vM [MPa] | 위치 | 능선 외/내 $\sigma_s$ | 능선 외/내 $\sigma_y$ |
|---|---|---|---|---|
| $N_x=10$ N/mm | 253 | 산 내측 = 골 외측 | −250 / +283 | −79 / +81 |
| $N_y=100$ N/mm | 132 | 균일 ($N_y/(h\lambda)$) | 0 | 132 / 132 |
| $N_{xy}=50$ N/mm | 159 | 능선 ($\tau=N_{xy}\alpha_1/h$ 성) | 0 | 0 |
| $M_x=10$ N·mm/mm | 150 | 산·골 | ∓167 | −54 / +44 |
| $M_y=100$ N·mm/mm | 216 | 산·골 내측 | 0 | 148 / 216 |
| $M_{xy}=20$ N·mm/mm | 471 | 능선 | 0 | 0 |

읽는 법: 능선 가로 인장 $N_x$는 $K_t=17$로 지배적 응력집중을 만든다(10 N/mm에 283 MPa). 능선 방향 굽힘 $M_y$는 $\sigma_y$ 막 성분($Ehf\kappa_{yy}$)이 능선에 집중되고 판 굽힘 성분은 얇은 원판 굽힘이라 작다. 비틀림 $M_{xy}$는 $D_{66}$이 작아 같은 모멘트에서 큰 곡률을 만들어 응력이 크다 — 쉐브론 회전 후 16/26 연성으로 $M_{xy}$가 항상 생기므로(07 §4.1) 두께 사이징에서 무시할 수 없다.

## 4. 합력 → 최대응력 절차 (1단계 수식 계산의 표준 경로)

1. 판축 하중(압력차, 포트 하중, 열응력)으로부터 판축 합력 $\{N_p, M_p\}$ — 폐형식(등가 판 이론) 또는 등가 판 FE.
2. 영역별 국부축 합력: `transform.resultants_to_local` ($\theta=\mp\beta$).
3. 국부축 거시 변형률: 채택식 $A^{-1}N$, $D^{-1}M$ (`recovery_from_resultants`).
4. VAM 복원 → 단면 위치·양면 $\sigma_s,\sigma_y,\tau_{sy}$ → von Mises 최대값과 위치, 막/굽힘 분리($N_s,M_s$).
5. 판정: 막 성분(1차 일반 막응력)과 막+굽힘(1차 국부 막 + 굽힘)을 각각 코드 허용치와 비교. 피로는 능선 표면 응력 진폭 사용.
6. apex·테두리·포트 보정계수(D단계) 곱.

## 5. 한계·방침

- **$t/R_{min}$**: 기본값 $t/R=0.6/1.5=0.4$. 쉘 이론(VAM 포함) 유효범위 $t/R\lesssim0.2$ 밖. 복원식은 $c_h=1+h^2\kappa^2/48$(기본값 능선에서 1.003)로 1차 보정만 있다. **방침**: 강성은 채택식 유지(Briassoulis FE ±2%로 확인됨), 능선 응력은 D-단계 솔리드 단위셀 FE와 대조해 $t/R$별 보정계수를 정하고 그 전까지는 복원값을 그대로 쓴다(두꺼운 원호는 곡률 내측 응력이 복원값보다 큰 방향이므로 보수적이지 않음을 명시).
- **비틀림 규약(E-1)**: §2의 GENS $D_{33}$/16/26 규약 확인 전까지 APDL 출력에는 "verify" 주석이 붙는다.
- **좌굴·접촉**: 복원식은 선형 응력만 준다(08 §6).
- **비대칭 단면**: $B_{ij}\ne0$이면 복원식에 $\mathcal B$, $\alpha_2$ 항 추가 필요(현재 미구현, 06 §6과 함께).

## 6. 다음 단계 연결

- **C-1**: YAML 입력 → `geometry` → `stiffness.adopted` → `transform` → `equivalent_plate`(등가 두께·APDL) + `stress_recovery`(합력 → 최대응력) 보고서 모듈.
- **C-3**: $K_t$, $A_{22}/A_{11}$, §3.3 단위 응답의 $f/t$, $H/p$, $\beta$ 의존성.
- **D**: 단위셀 6케이스(비틀림 포함) 쉘 FE로 채택식·복원식 검증, 솔리드 FE로 $t/R=0.4$ 보정.
- **E-1**: GENS 절 단일 요소 시험(6 케이스, 비틀림 규약), `plate_axes` vs `local_axes` 일치, SSPE 감도.
