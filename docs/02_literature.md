# 02. 문헌조사

주름형상 박판의 등가 직교이방성 모델에 관한 문헌을 정리한다. 본 프로젝트의 대상은 쉐브론 판이지만, 등가모델의 이론적 뼈대는 일방향 주름판(사인/원호/사다리꼴 단면) 연구에서 출발하므로, 먼저 일방향 주름판 모델을 정리하고 쉐브론 확장에 필요한 점을 별도로 기록한다.

## 1. 기준 문헌 — Lang & Su (2022)

**서지사항**
- Kun Lang, Mingzhou Su (Xi'an University of Architecture & Technology, School of Civil Engineering)
- *Equivalent orthotropic model for corrugated plates based on simplified constitutive relation*, Heliyon, 2022.
- 공개본: PMC9638743 / ScienceDirect S240584402202552X / SSRN 프리프린트 DOI 10.2139/ssrn.4046743 (2022-03).
- 키워드: corrugated plate, equivalent orthotropic model, simplified constitutive relation, equivalent bending stiffnesses, natural frequency prediction.

**요지 (초록 기준)**
- 기존 주름판 등가 직교이방성 모델은 Dirichlet 경계조건을 전제로 유도되어 실제 적용 시 정확도가 떨어지는 경우가 많다는 점을 문제로 제기한다.
- Kirchhoff 가정에 기초한 **단순화된 구성관계(simplified constitutive relation)** 로부터 등가 굽힘강성(equivalent bending stiffnesses)을 유도하여 개선된 등가 직교이방성 모델을 제안한다.
- 정적 응답 해석과 동적 응답 시험으로 검증하였으며, 저차·고차 고유진동수 예측 정확도가 향상되었다고 보고한다. 검증 사례에는 파형강판 아치 암거(corrugated metal pipe arch culvert) 실물 시험이 포함된다.

**본 프로젝트 관점의 시사점**
- 적용 대상이 토목 구조물(파형강판 암거)이므로 주름 단면은 일방향 원호/사인형이다. 쉐브론 판에 그대로 쓸 수 없고, 주름 방향 등가 강성을 구한 뒤 쉐브론 각도 β에 대한 좌표 변환과 V자 꺾임부의 영향을 추가로 다루어야 한다.
- 논문의 핵심은 "등가 굽힘강성 유도 방식"이므로, 본 프로젝트에서는 (1) 주름 방향(x)과 직각 방향(y)의 등가 D_x, D_y, D_xy 유도 절차와 (2) 면내 등가 강성(E_x, E_y, G_xy)의 취급 방식을 우선 확인한다.
- 검증이 고유진동수(강성 대표값) 중심이므로, 본 프로젝트의 목적인 **응력(두께 결정)** 평가에는 등가 강성 외에 형상 응력집중 환산이 별도로 필요하다.

**확인 필요 항목 (원문 열람 후 기입)**
- 본문의 등가 강성 수식(기호, 적분 형태) 및 단면 형상 파라미터 정의.
- 기존 모델(예: Samanta & Mukhopadhyay, Briassoulis, Xia et al.)과의 비교 오차 수치.
- 면내 강성 유도 포함 여부 및 막-굽힘 연성 취급.

## 2. 관련 문헌 (조사 예정)

| 구분 | 문헌 | 확인할 내용 |
|---|---|---|
| 일방향 주름판 등가 강성 | Briassoulis (1986), *Equivalent orthotropic properties of corrugated sheets* | 사인형 주름의 면내·굽힘 등가 강성 고전식 |
| 일방향 주름판 등가 강성 | Samanta & Mukhopadhyay (1999) | 사다리꼴 주름의 등가 강성 |
| 일방향 주름판 등가 강성 | Xia, Friswell, Saavedra Flores (2012), *Equivalent models of corrugated panels* | 일반 단면 주름에 대한 균질화 기반 등가 모델, FEA 검증 |
| 주름판 균질화 | Ye, Berdichevsky, Yu (2014) 등 균질화(MSG/VAM) 접근 | 임의 단면 주름의 등가 판 강성 |
| 쉐브론 판 (열교환기) | 판형 열교환기 쉐브론 판의 구조 강성/좌굴 연구 | 쉐브론 각도의 영향, 교차 적층 접촉점 거동 |
| 응력 평가 | 주름판 응력집중/국부 굽힘 연구 | 등가 공칭응력 → 실제 최대응력 환산 계수 |

## 3. 정리 방향

1. 일방향 주름판의 등가 직교이방성 강성식을 단면 형상별(사인/원호/사다리꼴)로 표로 정리한다.
2. 쉐브론 판은 "주름 방향이 ±β로 교차하는 두 영역"으로 보고, 각 영역의 등가 강성을 판 좌표계로 변환한 뒤 조합하는 방식과, 쉐브론 단위 셀을 직접 균질화하는 방식을 비교한다.
3. 응력 환산은 문헌값이 부족할 것으로 예상되므로, 단일 판 상세 FE로 직접 계수를 산정하는 것을 기본으로 한다.
