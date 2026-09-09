# 05. 쉐브론 전열판 형상 파라미터 정의서 (B-1)

등가모델(수식·엑셀·FE)에서 사용하는 **모든 형상 파라미터의 정의, 기호, 단위, 측정 기준, 파생 관계, 유효 범위**를 확정한다. 이 문서의 기호가 이후 `calc/`, `fea/`, 보고서의 표준이다. 값 자체는 `04_plate_dimension_form.md`(입력 양식)에 기입하고, 이 문서는 정의만 다룬다. 그림은 `fig/plate_dimensions.png` (a)~(e) 패널을 참조한다.

구현: `calc/geometry.py` (단면 라이브러리, 파생량 계산, 문헌값 자체검증), 입력 템플릿 `calc/plate_input_template.yaml`.

---

## 1. 좌표계와 방향 용어

세 개의 좌표계를 구분한다. 문헌·업계 용어가 서로 반대인 경우가 많으므로 **이 문서의 용어만 사용**한다.

| 좌표계 | 축 | 정의 | 그림 |
|---|---|---|---|
| **국부(단위셀) 좌표** $(x, y, z)$ | $x$ | 파형 진행 방향 = **능선에 수직**. 주름을 가로지르는 방향. 면내·굽힘 모두 **유연** | (b), (c) |
| | $y$ | **능선 방향**. 골·산이 뻗은 방향. **강함** | (b), (c) |
| | $z$ | 판 중앙면 법선 | (b) |
| **판(전역) 좌표** $(x_p, y_p, z_p)$ | $y_p$ | 유동축 = 두 포트 중심을 잇는 방향 | (a), (c) |
| | $x_p$ | 판면 내 $y_p$에 수직 | (a), (c) |
| | $z_p$ | 판 법선 ($=z$) | |
| **단면 곡선좌표** $(s, n, y)$ | $s$ | 단면 곡선을 따르는 호 길이. $\cos\theta = dx/ds$, $\sin\theta = dz/ds$ | (b) |
| | $n$ | 원판 두께 방향 법선 | |

- 쉐브론 각 $\beta$ = 능선($y$)과 유동축($y_p$) 사이 각, 판면 내 측정. apex 선 우측 영역은 $+\beta$, 좌측 영역은 $-\beta$ (그림 (c)).
- 국부→판 좌표 변환은 $z$축 회전 $\pm\beta$. 변환식은 B-3(`06_stiffness_transformation.md`, 예정)에서 다룬다.
- 문헌 용어 대응: Lang & Su·Xia의 "corrugated direction" = 본 문서 $x$, "transverse direction" = $y$. Ye의 $x_1 = x$, $x_2 = y$. 따라서 **문헌의 $D_{22}$("transverse bending stiffness")는 능선 방향 굽힘(강한 방향)**이다.

## 2. 1차 파라미터 (입력)

### 2.1 단면 (그림 (b)) — 단위셀 정의에 필수

| 기호 | 명칭 | 단위 | 측정 기준 | 문헌 기호 |
|---|---|---|---|---|
| `profile` | 단면 형식 | – | `sinusoidal` / `arc_tangent` / `trapezoidal` / `round` | – |
| $p$ | 피치 | mm | 능선에 **수직** 방향, 산–산 거리 | $2c$ (Lang & Su, Xia), $\varepsilon$ (Ye) |
| $H$ | 주름 깊이 | mm | **중앙면** 산–골 높이차. 외면 기준 전체 높이는 $H_o = H + t$ | $2f$, $2T$ (Ye) |
| $t$ | 판 두께 | mm | 성형 전 소재 두께(공칭). 성형 후 최소 두께는 $t_{min}$으로 별도 | $h$ |
| $R_c$ | 산 원호 반경 | mm | 중앙면 기준. `arc_tangent`, `round` | $R = r + h/2$ (Lang & Su) |
| $R_v$ | 골 원호 반경 | mm | 중앙면 기준. 미입력 시 $R_v = R_c$ | – |
| $\alpha$ | 플랭크 각 | ° | 경사부(접선부)와 $x$축(수평) 사이 각 | $\alpha$ |
| $T_L$ | 접선(직선) 길이 | mm | 산 원호 끝 ~ 골 원호 시작 | $T_L$ |
| $w_c, w_v$ | 산·골 평탄폭 | mm | `trapezoidal`만. 대칭이면 $w_c = w_v = c - H/\tan\alpha$ | – |

**입력 조합 규칙** (과잉 지정 금지, 하나의 조합만 준다):

| 형식 | 필수 | 선택/파생 |
|---|---|---|
| `sinusoidal` | $p, H, t$ | – |
| `arc_tangent` 모드 A | $p, H, t, R_c$ (+$R_v$) | $\alpha, T_L$ 파생 |
| `arc_tangent` 모드 B | $p, H, t, \alpha$ | $R_c = R_v = R$, $T_L$ 파생 |
| `trapezoidal` | $p, H, t, \alpha$ | $w$ 파생 |
| `round` | $R, L, t$ | $p = 4R$, $H = 2R + 2L$ 파생 (= `arc_tangent`의 $\alpha = 90°$ 특수해) |

`arc_tangent`의 기하 구속 (반주기: 산 원호 → 접선 → 골 원호):
$$\text{수평: } (R_c + R_v)\sin\alpha + T_L\cos\alpha = \frac p2,\qquad \text{수직: } (R_c + R_v)(1-\cos\alpha) + T_L\sin\alpha = H$$
두 식에서 미지수 2개를 정한다. $T_L \ge 0$ 조건에서 $\sin\alpha \le p/[2(R_c+R_v)]$, $1-\cos\alpha \le H/(R_c+R_v)$ — 이 범위를 벗어나면 "반경이 너무 큼"으로 입력 오류.

### 2.2 쉐브론 패턴 (그림 (a), (c), (e2))

| 기호 | 명칭 | 단위 | 정의 |
|---|---|---|---|
| $\beta$ | 쉐브론 각 | ° | §1 정의. 유동축 기준. 도면이 수평축 기준 $\beta'$를 주면 $\beta = 90° - \beta'$ |
| $n_{apex}$ | apex 선 수 | – | 단일 V = 1, W 패턴 = 2 이상 |
| $s_{apex}$ | apex 선 간격 | mm | $n_{apex} \ge 2$일 때 |
| `apex_join` | apex 연결 방식 | – | `continuous`(산–산 연속) / `offset`(반피치 어긋남, 산–골) |
| $r_{apex}$ | apex 꺾임부 필렛 반경 | mm | 또는 |
| $w_{apex}$ | apex 평탄 띠 폭 | mm | 있으면. 없으면 0 |
| `stack_rotation` | 인접 판 배열 | – | `rot180`(±β 교차, 표준) / `same` |

### 2.3 판 외형·테두리·포트 (그림 (a), (e1))

| 기호 | 명칭 | 단위 | 정의 |
|---|---|---|---|
| `plate_shape` | 판 외형 | – | `circular` / `rectangular` |
| $D_p$ ($L_p \times W_p$) | 판 외경 (길이×폭) | mm | 테두리 포함 |
| $D_e$ ($L_e \times W_e$) | 유효 주름 영역 | mm | 주름이 완전히 형성된 영역의 외곽 |
| $b_m$ | 외주 평탄 테두리 폭 | mm | 원형이면 $b_m = (D_p - D_e)/2 - w_{ro}$ |
| $w_{ro}$ | 런아웃 폭 | mm | 주름 깊이가 $H \to 0$으로 감소하는 전이 구간 |
| $h_{rim}$ | 테두리 오프셋 | mm | 테두리 평면과 주름 중앙면의 높이차 (부호: 채널 중심 쪽 +) |
| `rim_joint` | 외주 접합 | – | `laser_weld` / `gasket` / `braze`; 어느 판 쌍이 외주에서 접합되는지 |
| $w_{weld}$ | 용접 비드 폭 | mm | |
| $N_{port}, d_{port}, e_{port}$ | 포트 수·지름·중심거리 | – / mm / mm | 포트 간 거리 $= 2e_{port}$ |
| $b_{port}, h_{port}$ | 포트 테두리 폭·오프셋 | mm | 외주와 동일 정의 |
| `port_joint` | 포트 접합 | – | 외주와 동일 선택지 |
| $d_{flat,port}$ | 포트 주위 무주름 영역 지름 | mm | 있으면 |

### 2.4 적층 (그림 (d))

| 기호 | 명칭 | 단위 | 정의 |
|---|---|---|---|
| $N_p$ | 판 장수 | – | |
| $s_p$ | 판 피치 | mm | 인접 판 중앙면 간 거리. 산–산 접촉이면 $s_p = H + t$ |
| $b_{ch}$ | 채널 간극 | mm | 유로 높이. 고온측/저온측 각각 |
| `contact` | 접촉 조건 | – | `point` / `line`, 용접·브레이징 여부 |

### 2.5 재료·두께 보정

| 기호 | 명칭 | 단위 | 정의 |
|---|---|---|---|
| $E, \nu$ | 탄성계수·Poisson 비 | MPa, – | 설계온도 값 |
| $t_{min}$ | 성형 후 최소 두께 | mm | 위치 명시. 응력 평가 시 두께로 사용 (보수적) |
| $c_a$ | 부식 여유 | mm | 응력 평가 시 $t - c_a$ |
| $\Delta p, \Delta H, \Delta t$ | 제작 공차 | mm | 파라메트릭 범위 설정용 |

## 3. 파생 형상량 (계산)

모든 등가 강성식은 단면을 **네 개의 스칼라 $c, l, I_1, I_2$** (+$h$)로만 참조한다. 이것이 단면 형식과 등가모델 사이의 인터페이스다.

| 기호 | 정의 | 물리 의미 |
|---|---|---|
| $c = p/2$ | 반주기 | |
| $f = H/2$ | 반높이 | |
| $l$ | 반주기 전개 호 길이 $= \int_0^{c}\sqrt{1+z'^2}\,dx$ | 한 주기 전개길이 $2l$; $l/c$ = 전개 비(면적 확대 비) |
| $I_1 = \displaystyle\int_0^{2l}\left(\frac{dx}{ds}\right)^2 ds = \int_0^{2c}\frac{dx}{\sqrt{1+z'^2}}$ | 한 주기 | 판 자체 굽힘의 방향 여현 가중 길이 |
| $I_2 = \displaystyle\int_0^{2l} z^2\,ds$ | 한 주기, 중앙면 기준 ($\langle z\rangle = 0$) | 단면 2차 모멘트 관련: $I_u = hI_2/(2c)$ |

Ye (2014) 셀 평균과의 관계 ($\varepsilon = 2c$):
$$\langle\sqrt a\rangle = \frac lc,\qquad \left\langle\frac1{\sqrt a}\right\rangle = \frac{I_1}{2c},\qquad \varepsilon^2\langle\phi^2\sqrt a\rangle = \frac{I_2}{2c},\qquad \langle\varphi\mathcal A\rangle = \langle\phi^2\sqrt a\rangle \text{ (대칭 단면)}$$

단위폭당 단면 특성 (Lang & Su Eq. 28–30): $A_u = hl/c$, $I_u = hI_2/(2c)$, $W_u = I_u/(f + h/2)$.

### 3.1 단면별 폐형식

| 형식 | $l$ | $I_1$ | $I_2$ |
|---|---|---|---|
| `sinusoidal` | 수치적분 | 수치적분 | 수치적분 |
| `arc_tangent` (일반, $R_c \ne R_v$ 허용) | $(R_c+R_v)\alpha + T_L$ | $(R_c+R_v)\left(\alpha + \tfrac12\sin2\alpha\right) + 2T_L\cos^2\alpha$ | $2\,[I_{2,c} + I_{2,f} + I_{2,v}]$ (아래) |
| `trapezoidal` | $\dfrac{2f}{\sin\alpha} + w$ | $\dfrac{4f\cos^2\alpha}{\sin\alpha} + 2w$ | $\dfrac{4f^3}{3\sin\alpha} + 2f^2 w$,  $w = c - \dfrac{2f}{\tan\alpha}$ |
| `round` ($c = 2R$) | $\pi R + 2L$ | $\pi R$ | $\dfrac{4L^3}{3} + 2\pi L^2R + 8LR^2 + \pi R^3$ |

`arc_tangent`의 $I_2$ 성분 ($q = \alpha/2 + \sin2\alpha/4$, $a_c = f - R_c$, $a_v = f - R_v$):
$$I_{2,c} = R_c\left(a_c^2\alpha + 2a_cR_c\sin\alpha + R_c^2 q\right),\quad I_{2,v} = R_v\left(a_v^2\alpha + 2a_vR_v\sin\alpha + R_v^2 q\right),\quad I_{2,f} = \frac{z_1^3 - z_2^3}{3\sin\alpha}$$
$$z_1 = f - R_c(1-\cos\alpha),\qquad z_2 = -f + R_v(1-\cos\alpha)$$
$R_c = R_v$이면 Lang & Su Eq. (27)과 일치한다($I_1 = R(2\alpha + \sin2\alpha) + 2T_L\cos^2\alpha$). Lang & Su는 $I_2$ 첫 항을 수치적분으로 두었으나 위 폐형식으로 대체한다.

### 3.2 검증 (`python3 calc/geometry.py`)

| 검증 | 기준 | 결과 |
|---|---|---|
| 사인형 $l, I_1, I_2$ | Lang & Su Table 4 사인형 판 → `verify_langsu_tables.py` 수치 | 1e-4 이내 일치 |
| 사다리꼴 폐형식 | Lang & Su Eq. (15) | 일치 |
| `round` 폐형식 vs `arc_tangent`($\alpha = 90°$) | Xia Table 4 | 완전 일치 |
| `arc_tangent` $\alpha, T_L$ 역산 | Lang & Su Table 10 (200×55×3 암거판, $r = 53$, $T_L = 32.3$, $\alpha_0 = 45.187°$) | $\alpha$ 일치, $T_L$ 32.17 vs 32.3 (공표값 자체의 반올림 불일치 0.4%) |
| 비대칭 `arc_tangent` ($R_c \ne R_v$) 폐형식 | 중앙면 좌표 수치적분 | $l, I_1, I_2$ 1e-6 이내 |

## 4. 무차원 군과 유효 범위

| 무차원 군 | 정의 | 의미 / 유효 범위 |
|---|---|---|
| $H/p$ | 깊이비 | 주름 "깊음" 정도. Ye 얕은 주름 근사($\phi \ll 1$)는 $H/p \lesssim 0.1$에서만 유효. PHE는 보통 0.2~0.5 → **깊은 주름 식(Ye Eq. 19 / Lang & Su)** 필수 |
| $t/p$, $t/H$ | 박판비 | Ye 얇은 판 극한($h/\varepsilon \ll 1$)의 근거. PHE는 $t/p \sim 0.05$~0.1 → 선두항 식(Eq. 29 = Lang & Su) 오차 확인 필요 → C-3 파라메트릭에서 Eq. 19 전체식과 비교 |
| $t/R_{min}$ | 쉘 유효성 | 고전 쉘 이론(Ye §3) 가정: $t \ll R_{min}$. $t/R_{min} > 0.2$이면 두께 방향 응력 분포가 비선형 → 상세 FE에서 솔리드 요소 검증 권고 |
| $l/c$ | 전개비 | 면적 확대 비. $A_{22}, A_{66}, D_{11}, D_{66}$이 직접 의존 |
| $D_e/p$ (또는 $L_e/p$) | RVE 유효성 | 판 크기 / 주름 피치. Lang & Su Table 7: $N \ge 12$에서 처짐 오차 3% 이내, $N \ge 20$에서 1% 이내. **$D_e/p \ge 20$ 권고** |
| $\beta$ | 쉐브론 각 | 0° = 능선이 유동축과 평행(일방향 주름), 90° = 능선이 유동축에 수직(워시보드). 통상 25°~65° |

## 5. 기본값 세트 (값 미확정 시 파라메트릭용, **가정**)

실측값이 오기 전 계산 모듈 개발·검증에 쓰는 **가정값**이다. 일반 PHE/SPHX 판의 통상 범위이며 문헌 검증값이 아니다. 실제 값이 확정되면 이 절은 삭제한다.

| 파라미터 | 기본값 | 파라메트릭 범위 |
|---|---|---|
| `profile` | `arc_tangent` | + `sinusoidal` 비교 |
| $p$ | 9.0 mm | 6 ~ 14 mm |
| $H$ | 3.2 mm | 2 ~ 5 mm |
| $t$ | 0.6 mm | 0.4 ~ 1.0 mm |
| $R_c = R_v$ | 1.5 mm | 1.0 ~ 3.0 mm |
| $\beta$ | 60° | 25° ~ 65° |
| $E, \nu$ | 193 GPa, 0.3 (316L 상온) | 설계온도 값으로 교체 |

기본값의 파생량 (`calc/geometry.py` 출력): $\alpha = 44.3°$, $T_L = 3.36$ mm, $l = 5.680$ mm ($l/c = 1.262$), $I_1 = 7.261$ mm, $I_2 = 12.98$ mm³, $A_u = 0.757$ mm²/mm, $I_u = 0.865$ mm⁴/mm. 무차원: $H/p = 0.36$, $t/p = 0.067$, $t/R_{min} = 0.40$ (→ 두께 방향 응력 비선형 가능성, §4 참조).

## 6. 입력 파일 규약

- 형식: YAML (`calc/plate_input_template.yaml`), 단위 **mm–N–MPa**, 각도 **도(°)**.
- 키 이름은 이 문서의 기호를 ASCII로 옮긴 것 (`p, H, t, R_c, R_v, alpha_deg, T_L, beta_deg, D_p, D_e, b_m, w_ro, h_rim, ...`).
- 단면은 §2.1 입력 조합 규칙 중 하나만 채운다. 나머지 키는 비우거나 삭제.
- `geometry.from_dict()`가 §2.1 규칙으로 파생량을 계산하고, 규칙 위반(과잉 지정, 반경 과대)은 오류로 반환한다.

## 7. 이 정의서에서 확정한 결정사항

1. 모든 단면 치수는 중앙면 기준, 피치는 능선 수직 방향 — 도면이 다르면 입력 단계에서 변환한다.
2. 등가모델은 단면을 $c, l, I_1, I_2, h$로만 본다. 단면 형식 추가는 `geometry.py`에 폐형식 또는 수치적분 함수를 추가하는 것으로 끝난다.
3. 쉐브론 각은 유동축 기준 $\beta$ 하나로 통일하고, apex 좌우를 $\pm\beta$로 구분한다.
4. `round` 단면은 `arc_tangent`의 특수해로 취급한다(별도 식을 유지하되 검증용).
5. 비대칭 단면($R_c \ne R_v$)을 허용하되, 이 경우 Ye (2014)의 연성 강성 $B_{ij} \ne 0$ 가능성을 B-2에서 검토한다.
6. RVE 유효성 기준 $D_e/p \ge 20$을 채택하고, 미달 시 전체 판 FE 비교로 오차를 별도 산정한다.
