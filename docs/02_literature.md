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

### 1.3 유도 절차 — RVE 등가 에너지법 / 등가 힘법 (Xia et al. 2012 절차, Table 1; 케이스별 국부 내력 가정은 §3.2 참조)

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
4. **응력 평가**: 논문은 강성(변위·진동수) 검증만 수행했고 응력 환산은 없다. 두께 결정용 응력은 Briassoulis (1986) 응력집중식(§4.1)과 Ye (2014) 복원식(§2.6)으로 1차 추정하고 상세 FE로 보정한다(01_plan 3.1-4).
5. **등가 두께**: Eq. (16)처럼 원판 두께를 유지하면 면내 강성이 맞지 않는다. 본 프로젝트 2단계(조립체 FE)에서는 면내·굽힘 동시 만족이 필요하므로 preintegrated shell section 방식 또는 등가 두께·등가 탄성계수 2변수 방식을 검토한다.
6. **적용 한계**: RVE 기반이므로 판 크기 ≫ 주름 피치일 때 유효. 전열판은 피치 대비 판 크기가 충분히 크므로 문제 없으나, 가스켓/실링 부근 국부 영역에는 부적합.

## 2. Ye, Berdichevsky & Yu (2014) — 변분점근법(VAM) 등가 고전판 모델

**서지사항**
- Zheng Ye (Utah State Univ.), Victor L. Berdichevsky (Wayne State Univ.), Wenbin Yu (Purdue Univ.)
- *An equivalent classical plate model of corrugated structures*, Int. J. Solids Struct. 51 (2014) 2073–2083. DOI 10.1016/j.ijsolstr.2014.02.025
- 원문 확보: 2026-09-08 (저작권 문서이므로 저장소에는 넣지 않음. 프로젝트 공유 폴더에 보관)
- Lang & Su (2022)가 "현재 가장 타당한 등가 강성식"으로 평가한 문헌. 얕은/깊은 주름 모두 유효하고, **인장–굽힘 연성 강성 $B_{ij}$와 국부 변형률 복원식을 최초로 제시**.

### 2.1 기호 및 Lang & Su 기호와의 대응

| Ye 기호 | 정의 | Lang & Su 대응 |
|---|---|---|
| $x$ | 파형 진행 방향(corrugation direction), $y$는 능선 방향 | 동일 |
| $\varepsilon$ | 한 주기의 투영 길이 | $2c$ |
| $X=x/\varepsilon\in[-\tfrac12,\tfrac12]$ | 무차원 셀 좌표 | – |
| $\langle f\rangle=\int_{-1/2}^{1/2}f\,dX$ | 셀 평균 (= $\frac{1}{2c}\int_0^{2c}f\,dx$) | – |
| $x_3=\varepsilon\phi(X)$, $\langle x_3\rangle=0$ | 주름 단면 형상(중앙면), 평균 0으로 원점 조정 | $z(x)$ |
| $\varphi=d\phi/dX=dx_3/dx$, $a=1+\varphi^2$ | 기울기, $\sqrt a=ds/dx$ | $\cos\theta=1/\sqrt a$ |
| $\varphi'=d\varphi/dX$ | 기울기의 $X$ 도함수(곡률 관련) | – |
| $T$ | 주름 높이(중앙면 기준 rise) | $f$ |
| $h$ | 원판 두께 | $h$ |
| $S=\varepsilon\langle\sqrt a\rangle$ | 한 주기 전개 길이 — Eq. (4) | $2l$ |
| $I_y=h\varepsilon^2\langle\phi^2\sqrt a\rangle$ | 단위폭당 단면 2차 모멘트 — Eq. (4) | $hI_2/(2c)$ |
| $\mu=E/[2(1+\nu)]$ | 전단탄성계수 | $G$ |
| $\langle\sqrt a\rangle$ | | $l/c$ |
| $\langle 1/\sqrt a\rangle$ | | $I_1/(2c)$ |
| $\varepsilon^2\langle\phi^2\sqrt a\rangle$ | | $I_2/(2c)$ |

구성관계 Eq. (1): $\{N_{xx},N_{yy},N_{xy},M_{xx},M_{yy},M_{xy}\}$ ↔ $\{\epsilon_{xx},\epsilon_{yy},2\epsilon_{xy},\kappa_{xx},\kappa_{yy},2\kappa_{xy}\}$, 강성행렬 $[A\ B;\ B\ D]$에서 $A_{16}=A_{26}=B_{16}=B_{26}=D_{16}=D_{26}=0$. 변형률 정의 Eq. (104): $\epsilon_{xx}=v_{1,1},\ \epsilon_{yy}=v_{2,2},\ 2\epsilon_{xy}=v_{1,2}+v_{2,1},\ \kappa_{xx}=v_{3,11},\ \kappa_{yy}=v_{3,22},\ \kappa_{xy}=v_{3,12}$ ($v_i$: 등가 판 변위 = 셀 평균 변위). **전단 변형률·곡률은 공학 변형률($2\epsilon_{xy}, 2\kappa_{xy}$) 기준**이므로 $A_{66}, D_{66}$은 $N_{xy}=A_{66}\cdot 2\epsilon_{xy}$ 형태다 (Lang & Su Eq. (1)의 $\varepsilon_{xy},\kappa_{xy}$도 같은 값이 나오므로 동일 관례로 판단).

### 2.2 일반식 (임의 단면, 비대칭 포함) — Eq. (11)–(13)

보조량 — Eq. (12), (13):
$$\mathcal B=\frac{\langle\sqrt a\,\phi\rangle}{\langle\sqrt a\rangle},\qquad
\mathcal C=-12\langle\varphi\mathcal A\rangle\frac{\varepsilon^2}{h^2}-\left\langle\frac{1}{\sqrt a}\right\rangle,\qquad
\alpha_1=1\Big/\left\langle\frac{\sqrt a}{1+\frac{\varphi'^2h^2}{48\varepsilon^2a^3}}\right\rangle,\qquad
\alpha_2=\alpha_1\left\langle\frac{\frac{h^2\varphi'}{12\varepsilon a}}{1+\frac{\varphi'^2h^2}{48\varepsilon^2a^3}}\right\rangle$$
$$\mathcal A(X)=-\int_0^X\sqrt a\,\phi(Y)\,dY+\mathcal B\int_0^X\sqrt a\,dY$$

등가 강성 — Eq. (11):
$$\begin{aligned}
A_{11}&=\frac{E}{1-\nu^2}\frac{12\varepsilon^2\langle\varphi\mathcal A\rangle}{h\mathcal C^2}+\frac{Eh}{1-\nu^2}\left\langle\frac1{\sqrt a}\right\rangle\frac1{\mathcal C^2}, & A_{12}&=\nu A_{11}, & A_{22}&=Eh\langle\sqrt a\rangle+\nu^2A_{11}, & A_{66}&=\mu h\alpha_1,\\
B_{11}&=\frac{E}{1-\nu^2}\frac{12\varepsilon^3\langle\varphi\mathcal A\rangle}{h\mathcal C^2}\mathcal B+\frac{Eh}{1-\nu^2}\left\langle\frac1{\sqrt a}\right\rangle\frac{\mathcal B\varepsilon}{\mathcal C^2}\ (=\varepsilon\mathcal B A_{11}), & B_{12}&=\nu B_{11}, & B_{22}&=Eh\varepsilon\langle\sqrt a\,\phi\rangle+\nu^2B_{11}, & B_{66}&=\mu h\alpha_2,\\
D_{11}&=\frac{Eh^3}{12(1-\nu^2)}\left(\frac{12^2\varepsilon^4\mathcal B^2}{h^4\mathcal C^2}\langle\varphi\mathcal A\rangle+\frac1{\langle\sqrt a\rangle}\right)+\frac{Eh}{1-\nu^2}\frac{\varepsilon^2\mathcal B^2}{\mathcal C^2}\left\langle\frac1{\sqrt a}\right\rangle, & D_{12}&=\nu D_{11},\\
D_{22}&=Eh\varepsilon^2\langle\phi^2\sqrt a\rangle+\frac{Eh^3}{12}\left\langle\frac1{\sqrt a}\right\rangle+\nu^2D_{11}, &
D_{66}&=\frac{\mu h}{4}\left\langle\frac{\sqrt a}{3}h^2-\frac1{\sqrt a}\frac{\frac{h^4\varphi'^2}{12^2\varepsilon^2a^2}-a\alpha_2^2}{1+\frac{\varphi'^2h^2}{48\varepsilon^2a^3}}\right\rangle.
\end{aligned}$$

- 평판으로 퇴화($\varphi=\phi=0$, $\sqrt a=1$, $\mathcal C=-1$, $\alpha_1=1$, $\alpha_2=0$)시 고전 등방 판 강성으로 정확히 환원됨 — Eq. (21). 기존 식 중 Seydel(굽힘)·모든 면내 강성식은 이 검사를 통과하지 못한다.
- **연성 강성 $B_{ij}$는 일반적으로 0이 아니다.** $\phi(X)$가 기함수(대칭 주름: $\phi(-X)=-\phi(X)$)일 때만 $\mathcal B=0,\ \alpha_2=0$이 되어 모두 소멸 — Eq. (18).

### 2.3 대칭 주름 — Eq. (19) (본 프로젝트 사인/원호 단면에 해당)

$$\begin{aligned}
A_{11}&=\frac{E}{1-\nu^2}\frac{12\varepsilon^2\langle\varphi\mathcal A\rangle}{h\mathcal C^2}+\frac{Eh}{1-\nu^2}\left\langle\frac1{\sqrt a}\right\rangle\frac1{\mathcal C^2}, & A_{12}&=\nu A_{11}, & A_{22}&=Eh\langle\sqrt a\rangle+\nu^2A_{11}, & A_{66}&=\mu h\alpha_1,\\
B_{11}&=B_{12}=B_{22}=B_{66}=0,\\
D_{11}&=\frac{Eh^3}{12(1-\nu^2)}\frac1{\langle\sqrt a\rangle}, & D_{12}&=\nu D_{11}, & D_{22}&=Eh\varepsilon^2\langle\phi^2\sqrt a\rangle+\frac{Eh^3}{12}\left\langle\frac1{\sqrt a}\right\rangle+\nu^2D_{11}, & D_{66}&=\frac{\mu h}4\left\langle\frac{\sqrt a}3h^2-\frac1{\sqrt a}\frac{\frac{h^4\varphi'^2}{12^2\varepsilon^2a^2}}{1+\frac{\varphi'^2h^2}{48\varepsilon^2a^3}}\right\rangle.
\end{aligned}$$
(대칭이면 $\mathcal A(X)=-\int_0^X\sqrt a\,\phi\,dY$이고, 부분적분으로 $\langle\varphi\mathcal A\rangle=\langle\phi^2\sqrt a\rangle$ — 사인형 수치검증에서 두 값 일치 확인.)

### 2.4 점근 단순화

**얕은 주름** ($\phi\sim\delta\ll1$, $h/\varepsilon$ 임의) — Eq. (22), (23): $\mathcal B\sim\delta$, $\mathcal C=-\langle1/\sqrt a\rangle$, $\alpha_1=1/\langle\sqrt a\rangle$
$$A_{11}=\frac{Eh}{1-\nu^2}\frac1{\langle1/\sqrt a\rangle},\ A_{22}=Eh\langle\sqrt a\rangle+\nu^2A_{11},\ A_{66}=\frac{\mu h}{\langle\sqrt a\rangle},\ B_{11}=\frac{Eh}{1-\nu^2}\frac{\mathcal B\varepsilon}{\langle1/\sqrt a\rangle},\ B_{22}=Eh\varepsilon\langle\sqrt a\phi\rangle+\nu^2B_{11},\ B_{66}=\mu h\alpha_2,$$
$$D_{11}=\frac{Eh^3}{12(1-\nu^2)}\frac1{\langle\sqrt a\rangle},\ D_{12}=\nu D_{11},\ D_{22}=\frac{Eh^3}{12}\left\langle\frac1{\sqrt a}\right\rangle+\nu^2D_{11},\ D_{66}=\frac{\mu h^2}{12}\langle\sqrt a\rangle\ \text{(원문 표기; }h^3\text{의 오식으로 보임)}$$
→ Lang & Su Table 2 "Ye — shallow" 열과 일치 ($\frac{Eh}{1-\nu^2}\frac{1}{\langle1/\sqrt a\rangle}=\frac{1}{1-\nu^2}\frac{2cEh}{I_1}$… 단, Lang & Su는 $\frac{I_1Eh}{2c}$로 전사했음. $\langle1/\sqrt a\rangle=I_1/2c$이므로 정확한 전사는 $\frac{2cEh}{(1-\nu^2)I_1}$. **Lang & Su Table 2의 Ye-shallow $A_{11}, A_{12}, A_{22}$ 열은 역수가 뒤바뀐 오식**으로 판단됨. 사인형 수치: $2cEh/((1-\nu^2)I_1)=1.914\times10^8$ N/m vs $I_1Eh/(2c(1-\nu^2))=1.276\times10^8$; 어느 쪽도 얕은 주름 가정이 성립하지 않는 $f/c=0.34$ 사례라 47,613과 크게 다름.)

**얇은 판** ($h/\varepsilon\ll1$, 대부분의 주름판) — Eq. (24), (25): $\mathcal C\approx-12\langle\varphi\mathcal A\rangle\varepsilon^2/h^2$, $\alpha_1\approx1/\langle\sqrt a\rangle$
$$A_{11}=\frac{Eh^3}{12(1-\nu^2)\varepsilon^2\langle\varphi\mathcal A\rangle},\ A_{22}=Eh\langle\sqrt a\rangle,\ A_{66}=\frac{\mu h}{\langle\sqrt a\rangle},\ B_{11}=\frac{Eh^3\mathcal B}{12(1-\nu^2)\varepsilon\langle\varphi\mathcal A\rangle},\ B_{22}=Eh\varepsilon\langle\phi\sqrt a\rangle,\ B_{66}=\frac{\mu h^3}{12\varepsilon}\frac{\langle\varphi'/a\rangle}{\langle\sqrt a\rangle},$$
$$D_{11}=\frac{Eh^3}{12(1-\nu^2)}\left(\frac{\mathcal B^2}{\langle\varphi\mathcal A\rangle}+\frac1{\langle\sqrt a\rangle}\right),\ D_{12}=\nu D_{11},\ D_{22}=Eh\varepsilon^2\langle\phi^2\sqrt a\rangle,\ D_{66}=\frac{\mu h^3}{12}\langle\sqrt a\rangle$$
차수 관계 — Eq. (26)–(28): $A_{11}\sim A_{12}\sim(h/\varepsilon)^2A_{22}\sim(h/\varepsilon)^2A_{66}$; $D_{11}\sim D_{12}\sim(h/\varepsilon)^2D_{22}\sim D_{66}$; $B_{11}\sim B_{12}\sim(h/\varepsilon)^2B_{22}\sim B_{66}$.

**얇은 판 + 대칭 주름** — Eq. (29) (= Lang & Su Table 2 "Ye — symmetric" 열):
$$A_{11}=\frac{Eh^3}{12(1-\nu^2)\varepsilon^2\langle\phi^2\sqrt a\rangle}=\frac{1}{1-\nu^2}\frac{2cEh^3}{12I_2},\quad A_{22}=Eh\langle\sqrt a\rangle=\frac{Ehl}{c},\quad A_{66}=\frac{\mu h}{\langle\sqrt a\rangle}=\frac{cEh}{2l(1+\nu)},$$
$$D_{11}=\frac{Eh^3}{12(1-\nu^2)\langle\sqrt a\rangle}=\frac{cEh^3}{12l(1-\nu^2)},\quad D_{12}=\nu D_{11},\quad D_{22}=Eh\varepsilon^2\langle\phi^2\sqrt a\rangle=\frac{I_2Eh}{2c},\quad D_{66}=\frac{\mu h^3}{12}\langle\sqrt a\rangle=\frac{lEh^3}{24c(1+\nu)}$$

### 2.5 기존 식 판정 (§2.3 Discussion)

| 모델 | 굽힘 | 면내 |
|---|---|---|
| Seydel (1931) Eq. (5) | $D_{12}=0$만 틀림. $D_{11},D_{22}=EI_y,D_{66}$은 얇은 판 극한에서 정확 | – |
| 1960–70년대 통용식 Eq. (8) | – | $T^2:=\langle x_3^2\sqrt a\rangle/2$로 정의하면 정확 |
| Briassoulis (1986) Eq. (6), (9) | $D_{11},D_{12}$ 정확. **$D_{22},D_{66}$ 틀림** | **$A_{11},A_{66}$ 틀림** |
| Xia et al. (2012) Eq. (7), (10) | **$D_{22}$만 틀림** ($1/(1-\nu^2)$ 과대) | $A_{11}$은 고차항 무시 시 정확, $A_{22}$ 근사적으로 정확 |
| Briassoulis 국부응력 복원식 | 강성 절반이 틀리므로 복원식도 신뢰 불가 | |

### 2.6 국부 변형률 복원식 — Eq. (14)–(17), (20) ★ 응력 평가에 직접 활용 가능

등가 판 해석으로 얻은 $\epsilon_{yy},\kappa_{yy},\kappa_{xy}$와 $v_{i,jk}$로부터 원래 주름 쉘의 국부 막 변형률 $\gamma^0_{\alpha\beta}$·굽힘 변형률 $\rho^0_{\alpha\beta}$(국부 곡선좌표 1=$s$방향, 2=$y$)를 복원한다:
$$\gamma^0_{11}=c_1\sqrt a-\nu a(\epsilon_{yy}+x_3\kappa_{yy}),\qquad 2\gamma^0_{12}=\frac{\sqrt a\,c_2-\frac{h^2\varphi'\kappa_{xy}}{12\varepsilon a}}{1+\frac{\varphi'^2h^2}{48\varepsilon^2a^3}},\qquad \gamma^0_{22}=\epsilon_{yy}+x_3\kappa_{yy},$$
$$\rho^0_{11}=a\left(c_1\frac{12x_3}{h^2}+c_4\right)+\nu\sqrt a\,\kappa_{yy},\qquad 2\rho^0_{12}=\frac{-2\sqrt a\,\kappa_{xy}+\frac{\varphi'}{2\varepsilon a}c_2}{1+\frac{\varphi'^2h^2}{48\varepsilon^2a^3}},\qquad \rho^0_{22}=-\frac1{\sqrt a}\kappa_{yy}.$$
$$c_1=\frac{\varepsilon\mathcal B(v_{3,11}+\nu v_{3,22})-(v_{1,1}+\nu v_{2,2})}{\mathcal C},\quad c_2=\alpha_1(v_{1,2}+v_{2,1})-\alpha_2v_{3,12},\quad c_4=\frac{v_{3,11}+\nu v_{3,22}}{\langle\sqrt a\rangle}-\frac{12}{h^2}c_1\frac{\langle x_3\sqrt a\rangle}{\langle\sqrt a\rangle}$$
대칭 주름 — Eq. (20): $c_1=-\dfrac{v_{1,1}+\nu v_{2,2}}{\mathcal C},\quad c_2=\alpha_1(v_{1,2}+v_{2,1}),\quad c_4=\dfrac{v_{3,11}+\nu v_{3,22}}{\langle\sqrt a\rangle}$.

국부 응력합력은 쉘 변형에너지 밀도 Eq. (43) $U=\mu h[\sigma(a^{\alpha\beta}\gamma_{\alpha\beta})^2+a^{\alpha\beta}a^{\gamma\delta}\gamma_{\alpha\gamma}\gamma_{\beta\delta}]+\frac{\mu h^3}{12}[\sigma(a^{\alpha\beta}\rho_{\alpha\beta})^2+a^{\alpha\beta}a^{\gamma\delta}\rho_{\alpha\gamma}\rho_{\beta\delta}]$, $\sigma=\nu/(1-\nu)$의 구성관계로 얻고, 이어서 3D 응력으로 복원한다.

> **본 프로젝트 시사점**: 01_plan 3.1-4의 "공칭응력 × 형상 응력집중 계수" 접근을 **수식 기반 국부 응력 복원**으로 대체·보완할 수 있다. 주름 단면 내 위치 $X$별로 $\gamma^0,\rho^0$을 계산해 최대 표면 응력 $\sigma_s=\frac{E}{1-\nu^2}[(\gamma^0_{11}\pm\frac h2\rho^0_{11})+\nu(\gamma^0_{22}\pm\frac h2\rho^0_{22})]$ 등을 구하면 두께 결정용 응력을 상세 FE 없이 1차 추정 가능. 상세 FE는 이 복원식의 검증 용도로 위치가 바뀐다. (03_next_steps B4, D3에 반영)

### 2.7 수치 검증 — Table 1 (사인형 $\varepsilon=0.64$ m, $T=0.11$ m, $h=0.005$ m, $E=30$ GPa, $\nu=0.2$; Lang & Su Table 4와 동일 판)

| 성분 | Eqs. (5),(8) 고전식 | Xia et al. (2012) | VAPAS (Lee & Yu 2011) | **Ye (2014)** |
|---|---|---|---|---|
| $A_{11}$ (N/m) | 53805 | 47613 | 48152 | 47613 |
| $A_{12}$ (N/m) | 10761 | 9523 | 9630 | 9523 |
| $A_{22}$ (N/m) | 1.8708e8 | 1.8708e8 | 1.8692e8 | 1.8708e8 |
| $A_{66}$ (N/m) | 5.0113e7 | 5.0113e7 | 5.0097e7 | 5.0113e7 |
| $D_{11}$ (N·m) | 261.004 | 261.004 | 263.972 | 261.004 |
| $D_{12}$ (N·m) | 52.20 | 52.20 | 52.95 | 52.20 |
| $D_{22}$ (N·m) | 1025270 | 1068260 | 1022874 | 1025540 |
| $D_{66}$ (N·m) | 162.39 | 162.39 | 163.38 | 162.39 |

- **`calc/verify_ye_tables.py`로 Eq. (19) 전체식(⟨φ𝒜⟩, 𝒞, α₁ 포함)과 Eq. (5), (8), (6), (9)를 구현하여 위 표의 'Present'·'Eqs.(5),(8)' 열 및 Lang & Su Table 5의 Briassoulis 열을 재현함** (2026-09-08). Briassoulis $D_{66}$은 Ye 전사식 $Eh^3/(24(1+\nu))$로 130.21, Lang & Su Table 5는 260.42 — **원문 확인 결과(§4.1) Ye가 옳고 Lang & Su는 비틀림 관례($B_{xy}=2D_{66}$)를 변환하지 않은 오류**.
- 11개 주름 정사각 사인형 판의 등분포 압력 해석을 ANSYS 쉘 상세모델과 비교해 일치 확인(Ye 2013 박사논문에 수록).

**Table 2 — 비대칭 지수-사인 단면** $\phi(X)=\eta[e^{\sin2\pi X}-\langle e^{\sin2\pi X}\rangle]$, $\eta=0.1$, $\varepsilon=1$ m, $h=0.005$ m: $B_{11}=204.26$ N, $B_{12}=40.85$ N, $B_{22}=794841$ N (VAPAS 225.98 / 42.64 / 817802). 비대칭 주름에서는 $B_{22}$가 무시할 수 없는 크기. 쉐브론 판 단면이 대칭(사인/원호)이면 해당 없음.

### 2.8 Xia (2012) 식의 2차 문헌 전사 불일치 — **해결됨 (2026-09-08, §3.3 참조)**

Ye Eq. (10)이 전사한 Xia의 $A_{22}=\nu^2A_{11}+\frac S\varepsilon Eh\left(\frac1{1-\nu^2}-\frac{1-\nu^2}{4(1+\nu)^2}\right)$는 $\nu=0.2$에서 괄호항이 0.875가 되어 Ye Table 1의 Xia 열($1.8708\times10^8$)을 재현하지 못한다($1.64\times10^8$). Lang & Su Eq. (7)의 전사 $\frac{\nu^2}{1-\nu^2}\frac{2cEh^3}{I_1h^2+12I_2}+\frac{lEh}{c}$는 정확히 재현한다. 또 Ye의 Xia $A_{11}=\frac{Eh^3}{12(1-\nu^2)}\frac{1}{\langle1/\sqrt a\rangle h^2/12+I_y/h}$는 Lang & Su 전사와 동치다. → Xia 원문(§3.3 Table 2)으로 확인: **Lang & Su 전사가 원문과 일치하고 Ye Eq. (10)의 전사가 오류**다. 원식은 Ye 자신의 $A_{22}$와 동일하므로 Xia 면내 강성은 정확하다.

## 3. Xia, Friswell & Saavedra Flores (2012) — RVE 균질화 등가 직교이방성 모델

**서지사항**
- Yi Xia, Michael I. Friswell (Swansea Univ.), Erick I. Saavedra Flores (Swansea Univ. / Univ. de Santiago de Chile)
- *Equivalent models of corrugated panels*, Int. J. Solids Struct. 49(13) (2012) 1453–1462. DOI 10.1016/j.ijsolstr.2012.02.023
- 원문 확보: 2026-09-08 (저작권 문서이므로 저장소에는 넣지 않음)
- Lang & Su (2022)가 유도 절차를 그대로 차용한 모체 논문. 모핑 항공기 스킨용 주름 복합재 적층판이 동기이며, **원판이 직교이방성(복합재 적층)이어도 적용 가능**한 일반형으로 제시됨.

### 3.1 좌표계·형상 (Fig. 1, Eq. 1–4)

- 전역 $xyz$: $x$ = 파형 진행 방향(corrugation direction, 유연), $y$ = 능선 방향(transverse, 강), $z$ 법선. Lang & Su와 동일.
- 국부 $(s,n,y)$: $s$ = $xz$ 평면 내 접선 방향(호 길이), $n$ = 법선. $\cos\theta=dx/ds,\ \sin\theta=dz/ds$.
- 원판 직교이방성 주축이 $s, y$와 일치한다고 가정(그래야 강성행렬의 비대각 대부분이 0).
- 형상 입력은 $c$(반주기), $l$(반주기 전개길이), $I_1=\int_0^{2l}(dx/ds)^2ds$, $I_2=\int_0^{2l}z^2ds$ 네 개뿐 — Eq. (19). Lang & Su §1.1과 동일 정의.
- 등가 판은 연성 강성 $\mathbf B$를 무시한 직교이방성 Kirchhoff 판 — Eq. (5), (6).

### 3.2 유도 절차 — 6개 일반화 변형률 경계조건과 국부 내력 가정 (§2.3)

RVE(한 주기 $2c$ × 폭 $b$)에 Table 1의 단위 일반화 변형률을 하나씩 부과하고, 국부 내력 분포를 평형·대칭 논거로 가정한 뒤 (i) 등가 에너지법(대각항: $\bar A_{11}=2U/(2cb)$ 등) 또는 (ii) 등가 힘법(비대각항: 국부 내력의 주기 평균)으로 등가 강성을 정한다. **각 케이스의 국부 내력 가정이 곧 국부 응력 복원식**이 되므로 아래에 기록한다.

| Case | 부과 변형률 | 국부 변형률 가정 | 국부 내력 분포 | 결과 |
|---|---|---|---|---|
| 1 | $\bar\epsilon_x=1$ | $\epsilon_y=\gamma_{sy}=\kappa_y=\kappa_{sy}=0$ | $N_s=\bar N_x\frac{dx}{ds}$ (Eq. 11), $N_y=\frac{A_{12}}{A_{11}}N_s$ (12), $M_s=\bar N_x z$ (13), $M_y=\frac{D_{12}}{D_{11}}M_s$ (14) | $\bar A_{11}=\dfrac{2c}{I_1/A_{11}+I_2/D_{11}}$ (21), $\bar A_{12}=\dfrac{A_{12}}{A_{11}}\bar A_{11}$ (23) |
| 2 | $\bar\epsilon_y=1$ | $\gamma_{sy}=\kappa_y=\kappa_{sy}=0$ | $N_s=\bar A_{12}\frac{dx}{ds}$ (24), $\epsilon_s=\frac{1}{A_{11}}(\bar A_{12}\frac{dx}{ds}-A_{12})$ (25), $N_y=\frac{A_{12}\bar A_{12}}{A_{11}}\frac{dx}{ds}+\frac{A_{11}A_{22}-A_{12}^2}{A_{11}}$ (26) | $\bar A_{22}=\dfrac{A_{12}}{A_{11}}\bar A_{12}+\dfrac lc\dfrac{A_{11}A_{22}-A_{12}^2}{A_{11}}$ (28) |
| 3 | $\bar\gamma_{xy}=1$ | – | $N_{sy}=\bar N_{xy}=\bar A_{66}$ 일정 | $\bar A_{66}=\dfrac cl A_{66}$ (31) |
| 4 | $\bar\kappa_x=1$ | $\epsilon_s=\epsilon_y=\gamma_{sy}=\kappa_y=\kappa_{sy}=0$ | $M_s=\bar M_x$ (32), $M_y=\frac{D_{12}}{D_{11}}M_s$ (33) | $\bar D_{11}=\dfrac cl D_{11}$ (36), $\bar D_{12}=\dfrac{D_{12}}{D_{11}}\bar D_{11}$ (37) |
| 5 | $\bar\kappa_y=1$ | $\epsilon_s=0,\ \epsilon_y=z,\ \gamma_{sy}=0,\ \kappa_s=0,\ \kappa_y=\frac{dx}{ds},\ \kappa_{sy}=0$ | $N_y=A_{22}z$, $M_y=D_{22}\frac{dx}{ds}$ | $\bar D_{22}=\dfrac{1}{2c}[A_{22}I_2+D_{22}I_1]$ (40) |
| 6 | $\bar\kappa_{xy}=1$ | $\kappa_{sy}=1$, 나머지 0 | $M_{sy}=D_{66}$ 일정 | $\bar D_{66}=\dfrac lc D_{66}$ (43) |

- Case 1의 에너지 적분에서 $S_{11}+2\frac{A_{12}}{A_{11}}S_{12}+\frac{A_{12}^2}{A_{11}^2}S_{22}=\frac1{A_{11}}$, 굽힘도 동형으로 $\frac1{D_{11}}$ (Eq. 16, 17) → $U=\frac12 b\bar A_{11}^2\left[\frac{I_1}{A_{11}}+\frac{I_2}{D_{11}}\right]$ (18).
- **등가 판의 Poisson 비는 원판과 동일** ($\bar A_{12}/\bar A_{11}=A_{12}/A_{11}$, $\bar D_{12}/\bar D_{11}=D_{12}/D_{11}$).
- Lang & Su는 이 절차에서 원판 강성 $A_{11},A_{12},A_{22},D_{22}$만 $\tilde A_{11}=Eh,\ \tilde A_{12}=\nu Eh,\ \tilde A_{22}=Eh,\ \tilde D_{22}=Eh^3/12$로 바꿔 넣은 것이다(§1.4). Case 5의 $\epsilon_y=z$ 가정에서 $\sigma_s$가 생기지 않는다고 보면 $A_{22}\to Eh$가 되어 $D_{22}$의 $1/(1-\nu^2)$가 사라진다.

### 3.3 일반식 — Table 2 (임의 단면, 임의 직교이방성 원판)

| 성분 | Xia et al. (2012) | 등방성 원판 대입 ($A_{11}=\frac{Eh}{1-\nu^2}$ 등) |
|---|---|---|
| $\bar A_{11}$ | $\dfrac{2c}{I_1/A_{11}+I_2/D_{11}}$ | $\dfrac1{1-\nu^2}\dfrac{2cEh^3}{I_1h^2+12I_2}$ |
| $\bar A_{12}$ | $\dfrac{A_{12}}{A_{11}}\bar A_{11}$ | $\nu\bar A_{11}$ |
| $\bar A_{22}$ | $\dfrac{A_{12}}{A_{11}}\bar A_{12}+\dfrac lc\dfrac{A_{11}A_{22}-A_{12}^2}{A_{11}}$ | $\nu^2\bar A_{11}+\dfrac{lEh}{c}$ |
| $\bar A_{66}$ | $\dfrac cl A_{66}$ | $\dfrac{cEh}{2l(1+\nu)}$ |
| $\bar D_{11}$ | $\dfrac cl D_{11}$ | $\dfrac1{1-\nu^2}\dfrac{cEh^3}{12l}$ |
| $\bar D_{12}$ | $\dfrac{D_{12}}{D_{11}}\bar D_{11}$ | $\nu\bar D_{11}$ |
| $\bar D_{22}$ | $\dfrac1{2c}[I_2A_{22}+I_1D_{22}]$ | $\dfrac1{1-\nu^2}\dfrac{12I_2Eh+I_1Eh^3}{24c}$ |
| $\bar D_{66}$ | $\dfrac lc D_{66}$ | $\dfrac{lEh^3}{24c(1+\nu)}$ |

**§2.8의 의문 해소**: 원문 Table 2의 $\bar A_{22}$는 Lang & Su Eq. (7) 전사 $\frac{\nu^2}{1-\nu^2}\frac{2cEh^3}{I_1h^2+12I_2}+\frac{lEh}{c}$와 정확히 일치한다. Ye (2014) Eq. (10)의 전사 $\nu^2A_{11}+\frac S\varepsilon Eh\left(\frac1{1-\nu^2}-\frac{1-\nu^2}{4(1+\nu)^2}\right)$는 **오류**다(사다리꼴 예제에서 $1.556\times10^8$ vs 원문 $1.613\times10^8$, −3.6%). Ye의 "A22 is approximately correct" 판정은 이 잘못된 전사에 근거했으나, 원식 $\nu^2\bar A_{11}+lEh/c$는 Ye 자신의 $Eh\langle\sqrt a\rangle+\nu^2A_{11}$과 **동일**하므로 결론(정확)은 그대로다.

### 3.4 단면별 폐형식 — Table 3, 4

**사다리꼴** (경사각 $\alpha$, Table 3): $l=\frac{2f}{\sin\alpha}+c-\frac{2f}{\tan\alpha}$, $I_1=\frac{4f\cos^2\alpha}{\sin\alpha}+2c-\frac{4f}{\tan\alpha}$, $I_2=\frac{4f^3}{3\sin\alpha}+2f^2\left(c-\frac{2f}{\tan\alpha}\right)$ (Lang & Su Eq. 15와 동일).

**원호(round) 주름** (Table 4; 반원 반경 $R$ + 수직 직선 $L$, 반주기 폭 $c=2R$, $l=\pi R+2L$): $I_1=\pi R$, $I_2=\frac{4L^3}{3}+2\pi L^2R+8LR^2+\pi R^3$. → `calc/verify_xia_tables.py`에서 수치적분과 일치 확인. **SPHX 전열판 단면(원호+직선)에 가장 가까운 폐형식**이므로 계산 시트의 기본 단면형 후보.

Table 3·4에 병기된 기존 식 (Xia의 전사):

| 성분 | Samanta & Mukhopadhyay (1999), 등방 $E,\nu,t$ | Yokozeki et al. (2006), 원호 주름 |
|---|---|---|
| $\bar A_{11}$ | $\dfrac{2c}{I_2}\dfrac{Et^3}{12}$ | $\dfrac{4RD_{11}}{I_2}$ |
| $\bar A_{12}$ | $\nu\bar A_{11}$ | – |
| $\bar A_{22}$ | $\dfrac lc Et$ | $\dfrac{l}{2R}A_{22}$ |
| $\bar A_{66}$ | $\dfrac cl\dfrac{Et}{2(1+\nu)}$ | – |
| $\bar D_{11}$ | $\dfrac cl\dfrac{Et^3}{12}$ | $\dfrac{2R}{l}D_{11}$ |
| $\bar D_{12}$ | $0$ | – |
| $\bar D_{22}$ | $\dfrac{Et}{2c}I_2$ | $\left[I_2+\dfrac{(3\pi R+8L)}{12}t^2\right]\dfrac{A_{22}}{4R}$ |
| $\bar D_{66}$ | $\dfrac lc\dfrac{Et^3}{6(1+\nu)}$ | – |

> **주의 — S&M 열의 계수 불일치**: Table 3에 인쇄된 S&M 수식을 사다리꼴 예제에 대입하면 Table 5의 S&M 수치와 $\bar A_{11},\bar A_{12},\bar A_{22}$는 $(1-\nu^2)$배(=0.91), $\bar D_{22}$는 0.5배, $\bar D_{66}$은 4배 차이가 난다(`verify_xia_tables.py`). Table 5 수치는 $\bar A_{11}=\frac{2c}{I_2}\frac{Et^3}{12(1-\nu^2)}$, $\bar A_{22}=\frac lc\frac{Et}{1-\nu^2}$, $\bar D_{22}=\frac{Et}{c}I_2$, $\bar D_{66}=\frac lc\frac{Et^3}{24(1+\nu)}$에 해당한다. S&M 원문의 변형률·곡률 정의(공학 전단, $I_2$ 적분 구간) 차이로 추정되며, Lang & Su Table 6의 S&M 열은 Xia Table 5 수치를 그대로 옮긴 것이다. S&M 원문 확인 전까지는 **Table 5 수치 기준 형태**를 쓴다.

### 3.5 검증

**Table 5 — 사다리꼴** ($E=21$ GPa, $\nu=0.3$, $c=0.0508$ m, $f=0.0127$ m, $t=0.00635$ m, $\alpha=45°$, $b=1.016$ m). FEM = ANSYS SHELL63 단위셀(4131 노드, 4000 요소)에 Table 1의 일반화 변형률 경계조건을 부과하고 경계 반력·모멘트로 등가 강성 산출.

| 성분 | Xia | S&M | 단위셀 FEM | Xia vs FEM |
|---|---|---|---|---|
| $\bar A_{11}$ (MN/m) | 4.052 | 4.150 | 4.051 | 0.01% |
| $\bar A_{12}$ (MN/m) | 1.216 | 1.245 | 1.215 | 0.01% |
| $\bar A_{22}$ (MN/m) | 161.332 | 176.888 | 163.910 | −1.57% |
| $\bar A_{66}$ (MN/m) | 42.489 | 42.489 | 42.797 | −0.72% |
| $\bar D_{11}$ (N·m) | 407.917 | 371.205 | 406.512 | 0.35% |
| $\bar D_{12}$ (N·m) | 122.375 | 0 | 121.954 | 0.35% |
| $\bar D_{22}$ (kN·m) | 17.809 | 31.647 | 17.805 | 0.02% |
| $\bar D_{66}$ (N·m) | 208.032 | 208.032 | 207.96 | 0.04% |

- **중요한 해석**: 단위셀 FEM은 Xia와 같은 Dirichlet 경계조건을 쓰므로 $\bar D_{22}$가 Xia 식(17.809)과 일치한다. 즉 단위셀 검증은 Xia $\bar D_{22}$의 "상한" 성격을 드러내지 못한다. 전체 판 해석(Lang & Su Table 7: Xia 모델 처짐 −4~−14% 과소)과 VAM(Ye: 16.2 kN·m)이 이를 드러낸다. → 본 프로젝트 D단계 상세 FE 검증은 **단위셀 + 전체 판** 두 수준으로 해야 한다.
- Fig. 3, 4: Case 1·5의 국부 내력 $N_s, N_y, M_s, M_y$ 분포가 FEM과 매우 근접 → §3.2의 국부 내력 가정이 타당함을 뜻하며, 이는 국부 응력 추정에 그대로 쓸 수 있다.
- 10주기 사다리꼴 판(1.016 m × 1.016 m, 4변 고정, 69 kN/m², SHELL63 71,120 요소) vs 등가판: 처짐 근접, 상세모델이 약간 강함(파형 방향 고정단의 국부 보강 효과를 등가판이 못 잡음 — Fig. 6).

**Table 7 — 원호 주름 복합재** (AS4/3501-6 $[0/90]_s$, $E_1=148$, $E_2=10.5$, $G_{12}=5.61$ GPa, $\nu_{12}=0.3$; $R=L=3$ mm, $b=15$ mm; SHELL181 1344 요소): Xia vs FEM 오차 $\bar A_{11}$ 0.12%, $\bar A_{22}$ 0.00%, $\bar A_{66}$ 0.78%, $\bar D_{11}$ 0.29%, $\bar D_{22}$ 0.02%, $\bar D_{66}$ 0.02%. Yokozeki 식과도 근접. Fig. 7, 8: $w=4R$, $L$ 변화에 따른 8개 강성의 파라메트릭 곡선(1~5 mm).

### 3.6 본 프로젝트 관점의 시사점

1. **원판 직교이방성 허용**: 압연 방향 이방성이나 클래드 판을 쓰는 경우에도 원판 $A_{ij}, D_{ij}$만 바꿔 넣으면 된다. 계산 모듈 입력을 등방 $E,\nu,t$가 아니라 원판 ABD로 받도록 설계하면 확장성이 생긴다.
2. **국부 응력 1차 추정식**: Case 1 가정 $N_s=\bar N_x\cos\theta,\ M_s=\bar N_x z$에서 파형 진행 방향 인장 $\bar N_x$에 대한 원판 표면 응력은 $\sigma_s=\frac{\bar N_x\cos\theta}{t}\pm\frac{6\bar N_x z}{t^2}$, 최대는 $|z|=f$인 산·골에서 $\approx\bar N_x(\frac{1}{t}+\frac{6f}{t^2})$. Case 5 가정에서 능선 방향 굽힘 $\bar M_y$에 대해 $N_y=A_{22}z\bar\kappa_y,\ M_y=D_{22}\frac{dx}{ds}\bar\kappa_y$, $\bar\kappa_y=\bar M_y/\bar D_{22}$. Ye §2.6 복원식과 상호 검증할 것.
3. **원호 주름 폐형식**(Table 4)이 SPHX 단면에 가장 가깝다. 실제 판 단면(원호+직선 또는 사인)에 맞춰 $l, I_1, I_2$만 바꾸면 되므로 계산 시트에 단면 라이브러리(사인/사다리꼴/원호-직선/원호-접선)를 두고 선택하게 한다.
4. **검증 전략**: 단위셀 FE(Dirichlet BC)는 Xia식과 일치할 뿐 정확도를 판정하지 못한다. 전체 판 굽힘 해석과 비교해야 $\bar D_{22}$ 오차가 드러난다.

## 4. 고전 등가 강성식 통합 비교표 (Ye 2014 §2.1 전사, Ye 기호)

Briassoulis·Samanta & Mukhopadhyay 원문 미확보 상태에서, Ye (2014)가 전사한 형태를 2차 출처로 기록한다. Samanta & Mukhopadhyay (1999) 수식은 Xia (2012) Table 3 전사로 §3.4에 별도 기록했다(Ye는 Easley 인용식이라 "불완전·부정확"으로 분류).

| 성분 | Seydel (1931) Eq. (5) / 통용식 Eq. (8) | Briassoulis (1986) Eq. (6), (9) [사인형 가정] | Xia et al. (2012) Eq. (7), (10) | Ye (2014) 얇은판·대칭 Eq. (29) |
|---|---|---|---|---|
| $A_{11}$ | $\dfrac{Eh^3}{6(1-\nu^2)T^2}$ | $\dfrac{Eh^3}{h^2+6(1-\nu^2)T^2\left(\frac{S^2}{\varepsilon^2}-\frac{S}{2\pi\varepsilon}\sin\frac{2\pi S}{\varepsilon}\right)}$ | $\dfrac{Eh^3}{12(1-\nu^2)}\dfrac{1}{\langle1/\sqrt a\rangle\frac{h^2}{12}+\frac{I_y}{h}}$ | $\dfrac{Eh^3}{12(1-\nu^2)\varepsilon^2\langle\phi^2\sqrt a\rangle}$ |
| $A_{12}$ | $\nu A_{11}$ | $\nu A_{11}$ | $\nu A_{11}$ | $\nu A_{11}$ |
| $A_{22}$ | $\dfrac S\varepsilon Eh$ | $\dfrac S\varepsilon Eh$ | $\nu^2A_{11}+\dfrac S\varepsilon Eh(\cdots)$ — §2.8 참조 | $Eh\langle\sqrt a\rangle$ |
| $A_{66}$ | $\dfrac\varepsilon S\dfrac{Eh}{2(1+\nu)}$ | $\dfrac{Eh}{2(1+\nu)}$ | $\dfrac\varepsilon S\dfrac{Eh}{2(1+\nu)}$ | $\dfrac{\mu h}{\langle\sqrt a\rangle}$ |
| $D_{11}$ | $\dfrac\varepsilon S\dfrac{Eh^3}{12(1-\nu^2)}$ | 동일 | 동일 | $\dfrac{Eh^3}{12(1-\nu^2)\langle\sqrt a\rangle}$ (동일) |
| $D_{12}$ | $0$ | $\nu D_{11}$ | $\nu D_{11}$ | $\nu D_{11}$ |
| $D_{22}$ | $EI_y$ | $\dfrac{EhT^2}{2}+\dfrac{Eh^3}{12(1-\nu^2)}$ | $\dfrac{EI_y}{1-\nu^2}+\left\langle\dfrac1{\sqrt a}\right\rangle\dfrac{Eh^3}{12(1-\nu^2)}$ | $Eh\varepsilon^2\langle\phi^2\sqrt a\rangle=EI_y$ |
| $D_{66}$ | $\dfrac S\varepsilon\dfrac{Eh^3}{24(1+\nu)}$ | $\dfrac{Eh^3}{24(1+\nu)}$ | $\dfrac S\varepsilon\dfrac{Eh^3}{24(1+\nu)}$ | $\dfrac{\mu h^3}{12}\langle\sqrt a\rangle$ (동일) |

사인형 검증판($\varepsilon=0.64$, $T=0.11$, $h=0.005$, $E=30$ GPa, $\nu=0.2$) 기준 오차 (기준 = VAPAS):

| 성분 | Seydel/통용식 | Briassoulis | Xia | Lang & Su | Ye |
|---|---|---|---|---|---|
| $A_{11}$ | +11.7% | −17.7% | −1.1% | −1.1% | −1.1% |
| $A_{22}$ | +0.1% | +0.1% | +0.1% | −3.9% | +0.1% |
| $A_{66}$ | 0.0% | +24.8% | 0.0% | 0.0% | 0.0% |
| $D_{11}$ | −1.1% | −1.1% | −1.1% | −1.1% | −1.1% |
| $D_{22}$ | +0.2% | −11.2% | +4.4% | +0.3% | +0.3% |
| $D_{66}$ | −0.6% | −20% (관례 변환 후; L&S표의 +59%는 관례 오류) | −0.6% | −0.6% | −0.6% |

**채택 권고 (2026-09-08 기준, 2026-09-09 Briassoulis 원문 확인 후 유지)**: 굽힘 $D_{11},D_{12},D_{22},D_{66}$ 및 면내 $A_{11},A_{12},A_{66}$은 Ye (2014) Eq. (19) 또는 그 얇은판 극한 Eq. (29)(= Lang & Su)로 계산. $A_{22}$는 Ye Eq. (19) $Eh\langle\sqrt a\rangle+\nu^2A_{11}$ 사용(Lang & Su의 $(1-\nu^2)$ 계수 제외). 필요 형상 입력은 $\langle\sqrt a\rangle,\langle1/\sqrt a\rangle,\langle\phi^2\sqrt a\rangle$ 세 적분(= $l, I_1, I_2$)이며, 정밀 $A_{11}, D_{66}$에는 $\langle\varphi\mathcal A\rangle, \alpha_1$ 추가.

### 4.1 Briassoulis (1986) 원문 확인 (2026-09-09)

**서지사항**: D. Briassoulis (Univ. of Illinois, Agricultural Engineering), *Equivalent orthotropic properties of corrugated sheets*, Computers & Structures 23(2) (1986) 129–138. 스캔본 10쪽. 저작권 문서이므로 저장소에는 넣지 않음. 재계산: `calc/verify_briassoulis_tables.py`.

**내용 요약**: 기존(Davies, Easley, El-Atrouzy 등 1960–70년대) 등가 강성식을 검토하고, 9절점 Lagrangian 쉘 FE로 **한 주기 단위셀에 일정 변형률·곡률 상태를 부과**해 실제 강성을 구한 뒤(Table 1, 2), 불일치하는 항을 Castigliano 제2정리로 수정(Table 3 "Present"). 부록 B에서 능선부 국부 응력집중식을 유도.

**벤치마크 형상** (§4.1): 원호+접선 표준 주름, $c=2$ in, $f=0.21875$ in, $l=2.046$ in, $\alpha=17.54°$, $R=2.0208$ in, $t=0.25$ in, $E=30\times10^6$ psi, $\mu=0.3$. 얕고 두꺼운 주름: $f/c=0.109$, $l/c=1.029$, $f/t=0.875$, $t/R=0.124$. `geometry.arc_tangent(p=4, H=0.4375, R=2.0208)`로 $l=2.058$, $\alpha=17.7°$ 재현(공표값과 0.6% 차이, 반올림).

**기호 대응**: 그의 $D_x, D_\mu, D_y, D_{xy}$ = 면내 $A_{11}, A_{12}, A_{22}, A_{66}$; $B_x, B_\mu, B_y, B_{xy}$ = 굽힘 $D_{11}, D_{12}, D_{22}$, 그리고 비틀림. $x$ = 파형 진행 방향(우리와 동일). $\mu_1=\mu\,cE_x/(lE)$, $\mu_2=\mu$.

**$D_{66}$ 2배 문제 해결**: 그의 비틀림 강성은 $M_{xy}=B_{xy}\,w_{,xy}$ 관례이고 Eq. (14) $B_{xy}=Et^3/[12(1+\mu)]$는 **등방 평판의 비틀림 강성 그 자체**다. Ye/Xia/본 프로젝트 관례($M_{xy}=D_{66}\cdot2\kappa_{xy}$)로는 $D_{66}=B_{xy}/2=Et^3/[24(1+\nu)]$. 따라서 **Ye Eq. (6)의 전사가 옳고, Lang & Su Table 5의 260.42는 관례를 변환하지 않은 오류**(130.21이 맞음). 어느 쪽이든 Briassoulis 식에는 $\lambda=l/c$ 계수가 없다는 점은 같다.

**독립 FE 벤치마크 — 우리 모델과의 비교** (kips/in; $D_{66}$은 $2D_{66}=B_{xy}$로 환산):

| 성분 | Briassoulis FE | Xia | Lang & Su | Ye 얇은판 | **채택식** | 비고 |
|---|---|---|---|---|---|---|
| $A_{11}$ | 1440 | 1460 (1.014) | 1436 (0.997) | 1765 (1.225) | **1460 (1.014)** | Briassoulis Eq. (8): 1421 |
| $A_{12}$ | 439 | 438 (0.998) | 431 (0.981) | 529 (1.206) | **438 (0.998)** | |
| $A_{22}$ | 7709 | 7849 (1.018) | 7153 (0.928) | 7718 (1.001) | **7849 (1.018)** | |
| $A_{66}$ | 2800 | 2803 (1.001) | 동일 | 동일 | **2803 (1.001)** | |
| $D_{11}$ | 41.5 | 41.7 (1.005) | 동일 | 동일 | **41.7 (1.005)** | |
| $D_{12}$ | 12.5 | 12.5 (1.001) | 동일 | 동일 | **12.5 (1.001)** | 기존식은 0 |
| $D_{22}$ | 216.5 | 242.2 (1.119) | 220.4 (1.018) | 182.5 (0.843) | **224.2 (1.035)** | Briassoulis Eq. (11): 222.4 |
| $2D_{66}$ | 30.0 | 30.9 (1.031) | 동일 | 동일 | **30.9 (1.031)** | 그의 FE는 $\lambda$ 계수 없음을 지지 |

- 채택식은 면내·$D_{11}$·$D_{12}$에서 ±2%, $D_{22}$에서 +3.5%로 1986년 FE(한 주기 4×4 요소의 성긴 격자)와 일치. Xia $D_{22}$는 +12%로 여기서도 상한 과대가 드러남. 얇은판 선두항은 $f/t\approx0.9$인 두꺼운 주름에서 ±20% 오차 — **$h^2J_1$ 항을 생략하면 안 되는 경우**의 실례.
- $D_{66}$: 그의 FE는 평판값(λ 계수 없음)을 주고 VAPAS(Ye Table 1)는 λ 계수를 지지한다. 얕은 주름이라 차이가 3%에 불과해 그의 FE로는 판별이 어렵다. PHE 비율($\lambda=1.26$)에서는 26% 차이가 나므로 **D단계 상세 FE에서 비틀림 케이스를 반드시 포함**해 확정한다. 채택식은 VAPAS·Ye를 따라 λ 계수를 유지.
- 그의 A(면내) 결론(Eq. 8, 10)에 대한 Ye의 "wrong" 판정은 깊은 주름 점근에서의 판정이며, 그의 얕은 형상에서는 $A_{11}$ 1421 vs FE 1440(−1.3%)으로 실용상 문제가 없다. 다만 $A_{66}=Et/(2(1+\mu))$(λ 없음)는 $A_{66}$ FE 2800 vs 2885(+3%)로 그의 FE와도 어긋난다.

**응력집중식 (Appendix B)** — 본 프로젝트 응력 평가에 직접 사용:
- 일정 $\varepsilon_x$(파형 가로 인장) 상태: 국부 모멘트 $M_x=N_x^*z$ (Eq. 6), 응력 $\sigma_x=\dfrac{N_x^*}{t}\left(1\pm\dfrac{6z}{t}\right)$ (B1). 능선($z=f$) 내측 섬유에서 최대 $\sigma_{x,\max}=\dfrac{N_x^*}{t}\left(1+\dfrac{6f}{t}\right)$ (B2) → **$K_t=1+6f/t$**. 외측 섬유는 $f/t>1/3$이면 반대 부호 $\dfrac{N_x^*}{t}\left(1-\dfrac{6f}{t}\right)$ (B3).
- 일정 $\kappa_y$(능선 방향 굽힘) 상태: 균일 모멘트 $M$과 $z$에 비례하는 축력 $N_y=\dfrac{M_y^*}{I/A-f^2/2}f\sin\dfrac{\pi x}{c}$ (A6), $\sigma_y=\dfrac{6M_y^*}{t^2}\dfrac{\mp\frac{z}{t}2(1-\mu^2)\mp1}{(f/t)^26(1-\mu^2)+1}$ (B5). 얇은 판에서는 오히려 감소 계수, 두꺼운 판에서 최대 1.134배(Fig. 6).
- Xia §3.2 Case 1 가정($N_s=\bar N_x\cos\theta$, $M_s=\bar N_xz$)과 동일한 물리이며 능선($\theta=0$)에서 정확히 일치. **PHE 기본값($f=1.6$, $t=0.6$)에서 $K_t=17$** — 파형 가로 방향 막력이 조금만 있어도 능선 내측에 큰 굽힘응력이 생긴다는 뜻이며, $A_{11}$이 작은 이유와 같은 물리다. 두께 결정 시 이 항이 지배할 가능성이 크다.

## 5. 관련 문헌

Lang & Su (2022) 서론에서 정리한 계보와 본 프로젝트에 필요한 확인 사항. 원문 확보 여부를 표시한다.

| 구분 | 문헌 | 확인할 내용 | 원문 |
|---|---|---|---|
| 고전 | Huber (1923), Seydel (1931), Timoshenko & Woinowsky-Krieger (1959) | 최초 등가 강성식, 사인/사다리꼴 굽힘·비틀림 강성 | 미확보 — Seydel 굽힘식·1960년대 면내식은 Ye (2014) Eq. (5), (8)로 확보 (§3) |
| 일방향 주름판 등가 강성 | Briassoulis (1986), *Equivalent orthotropic properties of corrugated sheets*, Comput. Struct. 23(2) 129–138 | 사인형 주름의 면내·굽힘 등가 강성 고전식 + 능선 응력집중식. Table 5에서 $\bar A_{11}$ −17%, $\bar D_{22}$ −11% 편차 | **확보 (2026-09-09)** — §4.1. $D_{66}$ 관례 문제 해결, 독립 FE 벤치마크 및 응력집중식 $K_t=1+6f/t$ 확보 |
| 일방향 주름판 등가 강성 | Samanta & Mukhopadhyay (1999), Eng. Struct. 21(3) 277–287 | 사다리꼴 주름의 인장+굽힘 강성(Easley 1969/1975 인용). Table 6에서 $\bar D_{22}$ +95%, $\bar D_{12}=0$ | 미확보 — 수식은 Xia (2012) Table 3 전사로 확보 (§3.4). 인쇄 수식과 수치 사이 계수 불일치 있음 → 원문 확인은 우선순위 낮음 |
| RVE 균질화 | **Xia, Friswell, Saavedra Flores (2012)**, *Equivalent models of corrugated panels*, IJSS 49(13) 1453–1462 | 임의 단면 일반식(Eq. 7), 유도 절차 §2.3. 면내 강성의 기준식 | **확보 (2026-09-08)** — §3에 전사·검증 완료. $A_{22}$ 원식 확정, 직교이방성 원판 일반형, 원호 주름 폐형식 확보 |
| VAM 균질화 | **Ye, Berdichevsky, Yu (2014)**, *An equivalent classical plate model of corrugated structures*, IJSS 51 2073–2083 | 얕은/깊은 주름 모두 유효한 완전식, 연성 강성 $B$ 최초 제시. 현재 "가장 타당한" 식으로 평가됨 | **확보 (2026-09-08)** — §2에 전사·검증 완료. 국부 변형률 복원식 포함 |
| 차원축소 | Kolpakov & Kolpakov (2018 arXiv:1811.01718; 2020 IJES 154 103327) | 곡선보 문제로 축소한 최단순 형태, 임의 파형(비대칭 포함) 확장 | 미확보 (arXiv는 공개) |
| MSG | Deo & Yu (2021), IJSS 208–209 262–271 (MSG-TW) | 사다리꼴 기준값(Table 6) 출처 | 미확보 |
| 등가판 검증 | Aoki & Maysenhölder (2017), IJSS 108 11–23 | 사인/사다리꼴 자유판 고유진동수 실험·FE, 등가판 모델 고차 모드 한계 | 미확보 |
| 원호 주름 해석해 | Kress & Winkler (2010, 2011), Compos. Struct. 92, 93 | 원호 단면 주름의 해석적 하중응답 — SPHX 사인/원호 단면과 관련 | 미확보 |
| 원호 주름 등가식 | Yokozeki et al. (2006), Composites A 37 1578–1586 | 원호+직선 주름 복합재의 등가 강성 실험·해석 | 미확보 — 수식은 Xia (2012) Table 4 전사로 확보 (§3.4) |
| 쉐브론 판 (열교환기) | 판형 열교환기 쉐브론 판의 구조 강성/좌굴 연구 | 쉐브론 각도의 영향, 교차 적층 접촉점 거동 | 문헌 탐색 필요 |
| 응력 평가 | 주름판 응력집중/국부 굽힘 연구 | 등가 공칭응력 → 실제 최대응력 환산 계수 | 문헌 탐색 필요 |

## 6. 정리 방향

1. 일방향 주름판의 등가 직교이방성 강성식은 §1.5(Lang & Su 기호), §3.3(Xia 일반형), §4(Ye 기호 통합표)로 정리 완료. 채택식은 §4 말미의 권고를 따른다. 3대 기준 문헌(Lang & Su, Ye, Xia) 원문 모두 확보·전사·수치검증됨.
2. 쉐브론 판은 "능선 방향이 $\pm\beta$로 교차하는 두 영역"으로 보고, 각 영역의 등가 강성을 판 좌표계로 변환한 뒤 조합하는 방식과, 쉐브론 단위 셀을 직접 균질화하는 방식을 비교한다.
3. 응력 환산은 Briassoulis (1986) Appendix B의 능선 응력집중식(§4.1, `calc/stress_recovery.py`), Xia (2012) §3.2의 케이스별 국부 내력 가정, Ye (2014) §2.6의 국부 변형률 복원식을 1차 수단으로 삼고(셋을 상호 검증), 단일 판 상세 FE는 복원식 검증 및 쉐브론 꺾임부 보정계수 산정에 사용한다. 상세 FE 검증은 단위셀(일정 변형률·곡률 부과)과 전체 판 두 수준으로 하며, 비틀림 케이스로 $D_{66}$의 $\lambda$ 계수를 확정한다(§3.5, §4.1).
