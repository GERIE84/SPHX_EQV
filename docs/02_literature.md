# 02. 문헌조사

주름형상 박판의 등가 직교이방성 모델에 관한 문헌을 정리한다. 본 프로젝트의 대상은 쉐브론 판이지만, 등가모델의 이론적 뼈대는 일방향 주름판(사인/원호/사다리꼴 단면) 연구에서 출발하므로, 먼저 일방향 주름판 모델을 정리하고 쉐브론 확장에 필요한 점을 별도로 기록한다.

## 1. 기준 문헌 — Lang & Su (2022)

**서지사항**
- Kun Lang, Mingzhou Su (Xi'an University of Architecture & Technology, School of Civil Engineering)
- *Equivalent orthotropic model for corrugated plates based on simplified constitutive relation*, Heliyon 8 (2022) e11264. DOI 10.1016/j.heliyon.2022.e11264
- 공개본: PMC9638743 / ScienceDirect S240584402202552X / SSRN 프리프린트 DOI 10.2139/ssrn.4046743 (2022-03). CC BY-NC-ND 4.0.
- 원문 사본: `docs/ref/Lang_Su_2022_Heliyon_e11264.pdf` (12쪽, 2022-11-03 corrected proof)
- 키워드: corrugated plate, equivalent orthotropic model, simplified constitutive relation, equivalent bending stiffnesses, natural frequency prediction.

**요지**
- 기존 대표체적요소(RVE) 기반 등가 직교이방성 모델(Xia et al. 2012)은 Dirichlet(일반화 변형률) 경계조건을 사용하므로 Voigt 상한과 같은 성격을 가지며, 특히 **주름 능선 방향 굽힘강성 $D_{22}$를 과대평가**한다.
- Kirchhoff 박판 가정으로 원판(박판) 구성관계를 단순화한 뒤 Xia et al.의 RVE 절차(등가 에너지법 + 등가 힘법)를 그대로 적용하여 새로운 등가 강성식을 유도한다. 결과적으로 $D_{22}$에서 계수 $1/(1-\nu^2)$가 제거되며, 이는 변분점근법(VAM) 기반 Ye et al. (2014)의 단순화식 및 Kolpakov & Kolpakov (2018, 2020)와 일치한다.
- 검증: 사인형·사다리꼴 주름판의 등가 강성 비교(Table 5, 6), 4변 단순지지 등분포하중 중앙 처짐(Table 7), 고유진동수(Table 8), 원호-접선(arc-and-tangent) 주름판 단면 특성치와 규격값 비교(Table 9), 스웨덴 Enköping 파형강판 아치 암거 실물시험 재해석(5.3절).

### 1.1 좌표계 및 형상 파라미터 정의 (Fig. 1, 2, 6)

| 기호 | 정의 | 비고 |
|---|---|---|
| $x$ | **corrugated direction** — 파형이 반복되는 방향(주름 능선에 수직). 판이 유연한 방향 | Fig. 1 |
| $y$ | **transverse direction** — 주름 능선(골/산)이 뻗은 방향. 판이 강한 방향 | Fig. 1 |
| $z$ | 판 중앙면 법선 방향 | |
| $s, n$ | 국부 곡선좌표: $s$는 단면 곡선을 따르는 호 길이 좌표, $n$은 원판 두께 방향 법선. $\theta$는 $s$축과 $x$축 사이 각 ($\cos\theta = dx/ds$) | Fig. 1 |
| $h$ | 원판(박판) 두께 | 등가 판 두께로도 그대로 사용됨 (§1.5 참조) |
| $f$ | 주름 반높이(half-height). 파고 $=2f$ | Fig. 2 |
| $c$ | 반주기(half-period). 피치 $=2c$ | Fig. 2 |
| $l$ | 반주기의 전개 호 길이(expanded half-length). 한 주기의 전개 길이 $=2l$ | |
| $\alpha$ | 사다리꼴 주름의 경사면 각도 / 원호-접선 주름의 접선 각도 | Fig. 2b, 6 |
| $I_1$ | $I_1 = \displaystyle\int_0^{2l}\left(\frac{dx}{ds}\right)^2 ds = \int_0^{2l}\cos^2\theta\, ds$ | 한 주기 전개길이에 대한 적분 |
| $I_2$ | $I_2 = \displaystyle\int_0^{2l} z^2\, ds$ | 단면의 2차 모멘트 관련 ($I_u = hI_2/2c$) |

> **주의 (용어)**: 이 논문은 "corrugated direction"을 파형 진행 방향(능선에 수직), "transverse direction"을 능선 방향으로 쓴다. 따라서 논문의 "equivalent transverse bending stiffness"는 $D_{22}$, 즉 **능선 방향 축에 대한 굽힘(강한 방향)** 강성이다. 본 프로젝트 문서에서 "주름 방향"이라는 표현은 혼동을 일으키므로, 형상 파라미터 정의서(03_next_steps B1)에서 "능선 방향(y, 강)" / "파형 진행 방향(x, 유연)"으로 통일한다.

파형별 $l, I_1, I_2$:

- **사인형** $z(x) = f\sin(\pi x/c)$ — Eq. (14)
  $$l=\int_0^{c}\sqrt{1+z'^2}\,dx,\quad I_1=\int_0^{2c}\cos^2\theta\sqrt{1+z'^2}\,dx=\int_0^{2c}\frac{dx}{\sqrt{1+z'^2}},\quad I_2=\int_0^{2c}z^2\sqrt{1+z'^2}\,dx$$
  (해석해 없음, 수치적분)
- **사다리꼴** (경사각 $\alpha$, Fig. 2b) — Eq. (15)
  $$l=\frac{2f}{\sin\alpha}+c-\frac{2f}{\tan\alpha},\quad I_1=\frac{4f\cos^2\alpha}{\sin\alpha}+2c-\frac{4f}{\tan\alpha},\quad I_2=\frac{4f^3}{3\sin\alpha}+2f^2\left(c-\frac{2f}{\tan\alpha}\right)$$
- **원호-접선(arc-and-tangent)** (내측 곡률반경 $r$, 접선부 길이 $T_L$, 접선각 $\alpha$, $R=r+h/2$) — Eq. (26), (27)
  $$z(x)=\begin{cases}x\tan\alpha & 0\le x\le X_1\\ f-R+\sqrt{R^2-(X_2-x)^2} & X_1\le x\le X_2\end{cases},\quad X_1=\tfrac{T_L}{2}\cos\alpha,\ X_2=X_1+R\sin\alpha$$
  $$l=T_L+2R\alpha,\quad I_1=R(2\alpha+\sin 2\alpha)+2T_L\cos^2\alpha,\quad I_2=4\int_{S_1}^{S_2}\left(f-R+R\sqrt{1-\sin^2\tfrac{T_L/2-s+R\alpha}{R}}\right)^2ds+\frac{T_L^3}{6}\sin^2\alpha$$
  ($S_1=T_L/2,\ S_2=S_1+R\alpha$; $I_2$ 첫 항은 수치적분)

### 1.2 등가 직교이방성 판 구성관계 (Eq. 1, 2)

등가 판 좌표계에서 연성 강성(B 행렬)은 0이고,
$$\{N_x,N_y,N_{xy},M_x,M_y,M_{xy}\}^T=\mathrm{diag}\!\left(\begin{bmatrix}\bar A_{11}&\bar A_{12}&0\\ \bar A_{12}&\bar A_{22}&0\\0&0&\bar A_{66}\end{bmatrix},\begin{bmatrix}\bar D_{11}&\bar D_{12}&0\\ \bar D_{12}&\bar D_{22}&0\\0&0&\bar D_{66}\end{bmatrix}\right)\{\varepsilon_x,\varepsilon_y,\varepsilon_{xy},\kappa_x,\kappa_y,\kappa_{xy}\}^T$$
(논문은 등가 판 물성에 윗줄 $\bar{\ }$를 붙여 원판 물성과 구분한다. 이하 표에서는 윗줄을 생략한다.)

원판(등방성 박판)의 국부좌표 구성관계 — Eq. (3):
$A_{11}=A_{22}=\frac{Eh}{1-\nu^2},\ A_{12}=\frac{\nu Eh}{1-\nu^2},\ A_{66}=\frac{Eh}{2(1+\nu)},\ D_{11}=D_{22}=\frac{Eh^3}{12(1-\nu^2)},\ D_{12}=\frac{\nu Eh^3}{12(1-\nu^2)},\ D_{66}=\frac{Eh^3}{24(1+\nu)}$

### 1.3 유도 절차 — RVE 등가 에너지법 / 등가 힘법 (Xia et al. 2012 절차, Table 1)

- RVE = 한 주기($2c$) × 폭 $b$. 변형에너지 $U=\tfrac12\iint \mathbf N^T\mathbf S\mathbf N\,ds\,dy$ (Eq. 4), $\mathbf S=\mathbf K^{-1}$ (Eq. 5).
- 등가 판 에너지 $\bar U=\tfrac12(2c)b\{\bar\varepsilon,\bar\kappa\}^T[\bar A,\bar D]\{\bar\varepsilon,\bar\kappa\}$ (Eq. 6)와 같다고 두고, 6개의 단위 일반화 변형률 경계조건 $\{1,0,0,0,0,0\}^T$ … $\{0,0,0,0,0,1\}^T$을 순서대로 부과한다.
- 등가 에너지법: 대각항 $\bar A_{11}=2U/(2cb)$ 등. 등가 힘법: 비대각항 $\bar A_{12}=\bar N_y$ (조건 $\{1,0,\dots\}$일 때) 등.

Xia et al. (2012) 결과 — Eq. (7):

| 성분 | Xia et al. (2012) |
|---|---|
| $\bar A_{11}$ | $\dfrac{2c}{I_1/A_{11}+I_2/D_{11}}=\dfrac{1}{1-\nu^2}\dfrac{2cEh^3}{I_1h^2+12I_2}$ |
| $\bar A_{12}$ | $\dfrac{\nu}{1-\nu^2}\dfrac{2cEh^3}{I_1h^2+12I_2}$ |
| $\bar A_{22}$ | $\dfrac{\nu^2}{1-\nu^2}\dfrac{2cEh^3}{I_1h^2+12I_2}+\dfrac{lEh}{c}$ |
| $\bar A_{66}$ | $\dfrac{c}{l}A_{66}=\dfrac{cEh}{2l(1+\nu)}$ |
| $\bar D_{11}$ | $\dfrac{c}{l}D_{11}=\dfrac{1}{1-\nu^2}\dfrac{cEh^3}{12l}$ |
| $\bar D_{12}$ | $\dfrac{\nu}{1-\nu^2}\dfrac{cEh^3}{12l}$ |
| $\bar D_{22}$ | $\dfrac{1}{2c}\left[I_2A_{22}+I_1D_{22}\right]=\dfrac{1}{1-\nu^2}\dfrac{12I_2Eh+I_1Eh^3}{24c}$ |
| $\bar D_{66}$ | $\dfrac{l}{c}D_{66}=\dfrac{1}{1+\nu}\dfrac{lEh^3}{24c}$ |

### 1.4 구성관계 단순화 (본 논문의 핵심, §2.2, Eq. 8–13)

Kirchhoff 가정: (1) 두께 방향 수직변형률 $\varepsilon_n$ 무시, (2) $\tau_{ns},\tau_{yn},\sigma_n$은 $\tau_{sy},\sigma_s,\sigma_y$보다 작아 이들이 만드는 변형률은 무시(평형에는 필요). 이에 따라 3D Hooke 법칙을 평면응력 3식으로 축약(Eq. 8)한 뒤, 각 일반화 변형률 경계조건에서 국부 변형률 상태를 따져 **추가로 단순화**한다.

| 경계조건 | 국부 변형률 상태 | 단순화 결과 | 수정 원판 강성 |
|---|---|---|---|
| $\{1,0,0,0,0,0\}^T$ | $\varepsilon_y=0,\gamma_{sy}=0,\kappa_y=0,\kappa_{sy}=0$ → $\sigma_y,\tau_{sy}$가 만드는 변형률 무시 | $\varepsilon_s=\sigma_s/E=N_s/(Eh)$ (1축 Hooke) — Eq. (9) | $\tilde A_{11}=Eh,\ \tilde A_{12}=\nu\tilde A_{11}=\nu Eh$ |
| $\{0,1,0,0,0,0\}^T$ | $\varepsilon_s\approx0$으로 가정 | $\varepsilon_y=\sigma_y/E=N_y/(Eh)$ — Eq. (10) | $\tilde A_{22}=Eh$ |
| $\{0,0,0,0,1,0\}^T$ | $\{\varepsilon_s,\varepsilon_y,\gamma_{sy},\kappa_s,\kappa_y,\kappa_{sy}\}=\{0,z,0,dx/ds,0,0\}$ — Eq. (11) | $M_y=\dfrac{Eh^3}{12}\kappa_y$ — Eq. (12) | $\tilde D_{22}=\dfrac{Eh^3}{12}$ |

단순화된 원판 구성관계 — Eq. (13):
$$\mathbf K_{\text{simp}}=\mathrm{diag}\!\left(\begin{bmatrix}Eh&\nu Eh&0\\ \nu Eh&Eh&0\\0&0&\frac{Eh}{2(1+\nu)}\end{bmatrix},\begin{bmatrix}\frac{Eh^3}{12(1-\nu^2)}&\frac{\nu Eh^3}{12(1-\nu^2)}&0\\ \frac{\nu Eh^3}{12(1-\nu^2)}&\frac{Eh^3}{12}&0\\0&0&\frac{Eh^3}{24(1+\nu)}\end{bmatrix}\right)$$
즉 $A_{11},A_{12},A_{22}$에서 $1/(1-\nu^2)$가 제거되고 $D_{22}$가 $Eh^3/12$로 바뀐다. $A_{66},D_{11},D_{12},D_{66}$은 그대로다. 이 $\mathbf K_{\text{simp}}$를 §1.3 절차에 넣어 등가 강성을 얻는다.

### 1.5 등가 강성식 비교 — Table 2 (Xia 2012 / 본 논문 / Ye 2014)

| 성분 | Xia et al. (2012) [31] | **Lang & Su (2022)** | Ye et al. (2014) [35] — 대칭 주름 | Ye et al. (2014) — 얕은 주름 |
|---|---|---|---|---|
| $\bar A_{11}$ | $\dfrac{1}{1-\nu^2}\dfrac{2cEh^3}{I_1h^2+12I_2}$ | $\dfrac{2cEh^3}{I_1h^2+12I_2(1-\nu^2)}$ | $\dfrac{1}{1-\nu^2}\dfrac{2cEh^3}{12I_2}$ | $\dfrac{1}{1-\nu^2}\dfrac{I_1Eh}{2c}$ |
| $\bar A_{12}$ | $\dfrac{\nu}{1-\nu^2}\dfrac{2cEh^3}{I_1h^2+12I_2}$ | $\dfrac{2c\nu Eh^3}{I_1h^2+12I_2(1-\nu^2)}$ | $\dfrac{\nu}{1-\nu^2}\dfrac{2cEh^3}{12I_2}$ | $\dfrac{\nu}{1-\nu^2}\dfrac{I_1Eh}{2c}$ |
| $\bar A_{22}$ | $\dfrac{\nu^2}{1-\nu^2}\dfrac{2cEh^3}{I_1h^2+12I_2}+\dfrac{lEh}{c}$ | $\dfrac{2cE\nu^2h^3}{I_1h^2+12I_2(1-\nu^2)}+(1-\nu^2)\dfrac{lEh}{c}$ | $\dfrac{Ehl}{c}$ | $\dfrac{\nu^2}{1-\nu^2}\dfrac{I_1Eh}{2c}+\dfrac{lEh}{c}$ |
| $\bar A_{66}$ | $\dfrac{cEh}{2l(1+\nu)}$ | 동일 | 동일 | 동일 |
| $\bar D_{11}$ | $\dfrac{1}{1-\nu^2}\dfrac{cEh^3}{12l}$ | 동일 | 동일 | 동일 |
| $\bar D_{12}$ | $\dfrac{\nu}{1-\nu^2}\dfrac{cEh^3}{12l}$ | 동일 | 동일 | 동일 |
| $\bar D_{22}$ | $\dfrac{1}{1-\nu^2}\dfrac{12I_2Eh+I_1Eh^3}{24c}$ | $\dfrac{12I_2Eh+I_1Eh^3}{24c}=\dfrac{I_2Eh}{2c}+\dfrac{I_1Eh^3}{24c}$ | $\dfrac{I_2Eh}{2c}$ | $\dfrac{I_1Eh^3}{24c}+\dfrac{\nu^2}{1-\nu^2}\dfrac{cEh^3}{12l}$ |
| $\bar D_{66}$ | $\dfrac{1}{1+\nu}\dfrac{lEh^3}{24c}$ | 동일 | 동일 | 동일 |

**$\bar D_{22}$(능선 방향 굽힘) 표현식 비교 — Table 3** (논문 기호로 통일):

| 출처 | 논문 기호 표현 |
|---|---|
| Kolpakov & Kolpakov (2018) [36], (2020) [37], Ye et al. (2014) [35] | $\nu^2\bar D_{11}+\dfrac{I_2Eh}{2c}+\dfrac{I_1Eh^3}{24c}$ (세 식 동일) |
| Xia et al. (2012) [31] | $\dfrac{1}{1-\nu^2}\left(\dfrac{I_2Eh}{2c}+\dfrac{I_1Eh^3}{24c}\right)$ |
| **Lang & Su (2022)** | $\dfrac{I_2Eh}{2c}+\dfrac{I_1Eh^3}{24c}$ |

- 저자 해석: [35, 36, 37]의 제1항($\nu^2\bar D_{11}$, $h^3$ 차수)은 파형 형상 정보가 없고 작으므로 무시하면 본 논문 식과 같아진다. Xia 식은 $1/(1-\nu^2)>1$ 계수 때문에 항상 크다 → Dirichlet 경계조건 결과가 상한이라는 주장의 근거.
- 실용적으로 $\bar D_{22}\approx E\,I_u$ ($I_u=hI_2/2c$: 단위폭당 단면 2차 모멘트) + 판 자체 굽힘 기여 $I_1Eh^3/24c$.

### 1.6 등가 판 탄성상수 변환 및 등가 밀도 (Eq. 16, 17)

FE 입력용 직교이방성 탄성상수는 **원판 두께 $h$를 등가 판 두께로 그대로 사용**하고 굽힘강성 $\bar D$만 맞추도록 환산한다:
$$E_x=\frac{12(\bar D_{11}\bar D_{22}-\bar D_{12}^2)}{h^3\bar D_{22}},\quad E_y=\frac{12(\bar D_{11}\bar D_{22}-\bar D_{12}^2)}{h^3\bar D_{11}},\quad \nu_{xy}=\frac{\bar D_{12}}{\bar D_{22}},\quad G_{xy}=\frac{12\bar D_{66}}{h^3},\quad G_{zx}=\frac{E_x}{2(1+\nu)},\quad G_{yz}=\frac{E_y}{2(1+\nu_{xy})}$$
$$\rho_{eq}\,2cbh=\rho\,2lbh\ \Rightarrow\ \rho_{eq}=\rho\,\frac{l}{c}$$

> **본 프로젝트 관점 주의**: 이 변환은 굽힘(out-of-plane) 하중 전용이다. 두께를 $h$로 유지하면 면내 강성 $E_yh$는 $\bar A_{22}$와 일치하지 않는다(Table 11에서 $E_\theta=1.2\times10^8$ MPa처럼 비물리적으로 큰 값이 나오는 이유). 면내·굽힘을 동시에 만족시키려면 (i) 등가 두께 $t_{eq}$를 별도로 정의하거나 (ii) ANSYS 층상 쉘/사용자 정의 단면(SECTYPE,,GENS 또는 SHELL181 preintegrated section)으로 $\bar A,\bar D$를 직접 입력해야 한다. → 03_next_steps B4, E1에서 다룬다.

### 1.7 검증 결과 요약

**Table 4 검증용 판 파라미터**

| 파형 | $E$ (GPa) | $\nu$ | $\rho$ (kg/m³) | $c$ (m) | $f$ (m) | $h$ (m) | $\alpha$ |
|---|---|---|---|---|---|---|---|
| 사인형 | 30 | 0.2 | 7830 | 0.32 | 0.11 | 0.005 | – |
| 사다리꼴 | 21 | 0.3 | 7850 | 0.0508 | 0.0127 | 0.00635 | 45° |

**Table 5 사인형 등가 강성** (기준 = VAPAS / Ye et al.)

| 성분 | Briassoulis | Xia et al. | VAPAS [48] | Ye et al. | **Lang & Su** |
|---|---|---|---|---|---|
| $\bar A_{11}$ (N/m) | 39639 | 47613 | 48152 | 47613 | 47612 |
| $\bar A_{12}$ (N/m) | 7928 | 9523 | 9630 | 9523 | 9522 |
| $\bar A_{22}$ (MN/m) | 187.08 | 187.08 | 186.92 | 187.08 | 179.60 |
| $\bar A_{66}$ (MN/m) | 62.500 | 50.113 | 50.097 | 50.113 | 50.113 |
| $\bar D_{11}$ (N·m) | 261.004 | 261.004 | 263.972 | 261.004 | 261.004 |
| $\bar D_{12}$ (N·m) | 52.20 | 52.20 | 52.95 | 52.20 | 52.20 |
| $\bar D_{22}$ (N·m) | 907830 | 1068260 | 1022874 | 1025540 | 1025530 |
| $\bar D_{66}$ (N·m) | 260.42 | 162.39 | 163.38 | 162.39 | 162.39 |

**Table 6 사다리꼴 등가 강성** (기준 = MSG-TW)

| 성분 | Samanta & Mukhopadhyay | MSG-TW [40] | Xia et al. | **Lang & Su** |
|---|---|---|---|---|
| $\bar A_{11}$ (MN/m) | 4.150 | 4.052 | 4.052 | 4.042 |
| $\bar A_{12}$ (MN/m) | 1.245 | 1.216 | 1.216 | 1.213 |
| $\bar A_{22}$ (MN/m) | 176.888 | 161.332 | 161.332 | 146.844 |
| $\bar A_{66}$ (MN/m) | 42.489 | 42.489 | 42.489 | 42.489 |
| $\bar D_{11}$ (N·m) | 371.205 | 407.913 | 407.917 | 407.917 |
| $\bar D_{12}$ (N·m) | 0 | 122.375 | 122.375 | 122.375 |
| $\bar D_{22}$ (N·m) | 31647 | 16242 | 17809 | 16206 |
| $\bar D_{66}$ (N·m) | 208.032 | 208.033 | 208.032 | 208.032 |

- Lang & Su의 $\bar D_{22}$는 Ye/VAPAS/MSG-TW와 0.2% 이내로 일치. Xia는 사인형 +4.2%, 사다리꼴 +9.6% 과대.
- 반면 Lang & Su의 **$\bar A_{22}$(능선 방향 면내 강성)는 과소평가**: 사인형 −3.7%(Ye 대비), 사다리꼴 −9.0%(MSG-TW 대비). 저자도 "면내 강성 정확도는 충분히 높지 않으나 굽힘강성 기반 탄성상수 산정에는 사용 가능"이라고 인정.
- **Table 5·6의 Xia/Lang & Su 열 전체를 `calc/verify_langsu_tables.py`로 재계산하여 논문값을 소수점까지 재현함** (2026-09-08). 위 수식 전사에 오류가 없음을 확인.

**Table 7 — 4변 단순지지, 등분포 1000 N/m², 중앙 처짐 (사다리꼴, ABAQUS S4R 대비)**: $N$(주름 수)=8~20, $\gamma=a/b$=0.5~2. Xia 모델은 항상 처짐 과소(−4~−14%). Lang & Su는 $N=8$에서 −5.5~+2.2%, $N\ge12$에서 ±2.6% 이내, $N=20$에서 ±1% 이내. → RVE 모델이므로 판 크기 ≫ 단일 파형일 때 정확도 급상승.

**Table 8 — 고유진동수 (사다리꼴, $N=20$, 4변 단순지지)**: $\omega_{mn}=\pi^2\sqrt{[D_1(m/a)^4+D_2(n/b)^4+2H(m/a)^2(n/b)^2]/(\rho_{eq}h)}$, $D_1=\bar D_{11}, D_2=\bar D_{22}, H=\bar D_{12}+2\bar D_{66}$ (Eq. 25). Lang & Su 오차: 저차 모드 |0.5%| 이내, $\omega_{17}$(고차, 능선 방향 7반파)에서 +4.8% (Xia +9.0%).

**Table 9 — 원호-접선 주름판(150×50 mm, GB/T 34567-2017) 단면 특성치**: $A_u=hl/c$, $I_u=hI_2/2c$, $W_u=I_u/(f+h/2)$ (Eq. 28–30). 규격값 대비 $A_u$ 오차 0%, $I_u$는 $h$=2~10 mm에서 0.15~5.76%.

**5.3 실물 암거 (Enköping, 200×55×3 mm 주름, 스팬 6 m, PLAXIS 3D + HSs 토질모델)**: 크라운 처짐 측정 65 mm vs 본 모델 64 / 공학법[46] 60 mm; 추력 −120 vs 112 / 106 kN/m; 모멘트 10.3 kNm/m는 양쪽 모두 근접.

### 1.8 본 프로젝트 관점의 시사점

1. **채택 수식**: 굽힘강성 $\bar D_{11},\bar D_{12},\bar D_{22},\bar D_{66}$은 Lang & Su(=Ye 단순화식)를 기준으로 삼는다. 면내 강성은 Lang & Su가 $\bar A_{22}$를 과소평가하므로 **Xia et al. (2012) 또는 Ye et al. (2014) 식을 사용**하는 것이 타당하다(두 식은 $\bar A_{11},\bar A_{12},\bar A_{66}$에서 실질적으로 같다).
2. **필요한 형상 입력은 $c, l, I_1, I_2, h$ 다섯 개**뿐이다. 임의 단면(사인/원호/사다리꼴/원호-접선)은 $l, I_1, I_2$ 수치적분으로 처리 가능 → 계산 시트의 입력 구조를 이렇게 잡는다.
3. **쉐브론 확장**: 논문은 일방향 주름만 다룬다. 쉐브론 판은 능선이 판 축에 대해 $\pm\beta$로 기울어진 두 영역의 조합이므로, 위 $\bar A,\bar D$를 $\pm\beta$ 회전 변환한 뒤 조합하는 절차(03_next_steps B3)가 별도로 필요하다. 회전 후에는 $\bar A_{16},\bar A_{26},\bar D_{16},\bar D_{26}$이 생기며, $\pm\beta$ 영역이 같은 비율이면 판 평균에서는 상쇄된다(국부적으로는 남음).
4. **응력 평가**: 논문은 강성(변위·진동수) 검증만 수행했고 응력 환산은 없다. 두께 결정용 응력은 별도 상세 FE로 환산계수를 산정해야 한다(01_plan 3.1-4).
5. **등가 두께**: Eq. (16)처럼 원판 두께를 유지하면 면내 강성이 맞지 않는다. 본 프로젝트 2단계(조립체 FE)에서는 면내·굽힘 동시 만족이 필요하므로 preintegrated shell section 방식 또는 등가 두께·등가 탄성계수 2변수 방식을 검토한다.
6. **적용 한계**: RVE 기반이므로 판 크기 ≫ 주름 피치일 때 유효. 전열판은 피치 대비 판 크기가 충분히 크므로 문제 없으나, 가스켓/실링 부근 국부 영역에는 부적합.

## 2. 관련 문헌

Lang & Su (2022) 서론에서 정리한 계보와 본 프로젝트에 필요한 확인 사항. 원문 확보 여부를 표시한다.

| 구분 | 문헌 | 확인할 내용 | 원문 |
|---|---|---|---|
| 고전 | Huber (1923), Seydel (1931), Timoshenko & Woinowsky-Krieger (1959) | 최초 등가 강성식, 사인/사다리꼴 굽힘·비틀림 강성 | 미확보 |
| 일방향 주름판 등가 강성 | Briassoulis (1986), *Equivalent orthotropic properties of corrugated sheets*, Comput. Struct. 23(2) 129–138 | 사인형 주름의 면내·굽힘 등가 강성 고전식. Table 5에서 $\bar A_{11}$ −17%, $\bar D_{22}$ −11%, $\bar D_{66}$ +60% 편차 | 미확보 (수치는 Table 5로 확보) |
| 일방향 주름판 등가 강성 | Samanta & Mukhopadhyay (1999), Eng. Struct. 21(3) 277–287 | 사다리꼴 주름의 인장+굽힘 강성. Table 6에서 $\bar D_{22}$ +95%, $\bar D_{12}=0$ | 미확보 (수치는 Table 6으로 확보) |
| RVE 균질화 | **Xia, Friswell, Saavedra Flores (2012)**, *Equivalent models of corrugated panels*, IJSS 49(13) 1453–1462 | 임의 단면 일반식(Eq. 7), 유도 절차 §2.3. 면내 강성의 기준식 | 미확보 — **우선 확보 필요** (유도 절차 상세, 검증 FE) |
| VAM 균질화 | **Ye, Berdichevsky, Yu (2014)**, *An equivalent classical plate model of corrugated structures*, IJSS 51 2073–2083 | 얕은/깊은 주름 모두 유효한 완전식, 연성 강성 $B$ 최초 제시. 현재 "가장 타당한" 식으로 평가됨 | 미확보 — **우선 확보 필요** (일반식은 Lang & Su 기호로 표현 불가) |
| 차원축소 | Kolpakov & Kolpakov (2018 arXiv:1811.01718; 2020 IJES 154 103327) | 곡선보 문제로 축소한 최단순 형태, 임의 파형(비대칭 포함) 확장 | 미확보 (arXiv는 공개) |
| MSG | Deo & Yu (2021), IJSS 208–209 262–271 (MSG-TW) | 사다리꼴 기준값(Table 6) 출처 | 미확보 |
| 등가판 검증 | Aoki & Maysenhölder (2017), IJSS 108 11–23 | 사인/사다리꼴 자유판 고유진동수 실험·FE, 등가판 모델 고차 모드 한계 | 미확보 |
| 원호 주름 해석해 | Kress & Winkler (2010, 2011), Compos. Struct. 92, 93 | 원호 단면 주름의 해석적 하중응답 — SPHX 사인/원호 단면과 관련 | 미확보 |
| 쉐브론 판 (열교환기) | 판형 열교환기 쉐브론 판의 구조 강성/좌굴 연구 | 쉐브론 각도의 영향, 교차 적층 접촉점 거동 | 문헌 탐색 필요 |
| 응력 평가 | 주름판 응력집중/국부 굽힘 연구 | 등가 공칭응력 → 실제 최대응력 환산 계수 | 문헌 탐색 필요 |

## 3. 정리 방향

1. 일방향 주름판의 등가 직교이방성 강성식은 §1.5 표(Xia / Lang & Su / Ye)로 1차 정리 완료. Briassoulis, Samanta & Mukhopadhyay 식은 원문 확보 후 같은 기호로 추가한다.
2. 쉐브론 판은 "능선 방향이 $\pm\beta$로 교차하는 두 영역"으로 보고, 각 영역의 등가 강성을 판 좌표계로 변환한 뒤 조합하는 방식과, 쉐브론 단위 셀을 직접 균질화하는 방식을 비교한다.
3. 응력 환산은 문헌값이 부족할 것으로 예상되므로, 단일 판 상세 FE로 직접 계수를 산정하는 것을 기본으로 한다.
