"""엑셀 계산 시트(calc/SPHX_EQV_calc.xlsx) Input 시트에 삽입하는 입력 위치 안내 그림 3장.

  xl_fig1_section.png  단면 입력: p, H (mid / outer), t, R_c, R_v, α, T_L, 산·골, +z/−z 면, t_min·CA → t_calc
  xl_fig2_plan.png     평면 입력: β (flow 기준 / horizontal 기준), 우·좌 영역 (θ = ∓β), 피치 기준 (normal / axis), 판축·국부축
  xl_fig3_loads.png    하중 입력: 국부축·판축 합력 N, M 의 방향과 부호
  xl_fig4_stress_points.png  응력 출력 위치 (산·골 × +z/−z 면, Stress 시트 열 대응)
라벨은 시트 항목 기호와 동일하게 쓴다 (한글 폰트가 없어 영문·기호). 한글 설명은 시트 셀에 있다.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Arc, Polygon, FancyBboxPatch, Rectangle
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "mathtext.fontset": "dejavusans"})
BLUE, RED, GREY, GREEN, ORANGE, PURPLE = "#1f4e79", "#c0392b", "#7f8c8d", "#1e8449", "#d68910", "#7d3c98"
YEL = "#fff3b0"


def tag(ax, x, y, text, color, fs=11, ha="center", va="center", bold=True):
    """노란 배경 라벨 = 시트 입력 항목."""
    ax.text(x, y, text, ha=ha, va=va, color=color, fontsize=fs, fontweight="bold" if bold else "normal",
            bbox=dict(boxstyle="round,pad=0.25", fc=YEL, ec=color, lw=0.8))


def dim_h(ax, x1, x2, y, text, color=RED, off=0.05, fs=11, tagged=True):
    ax.annotate("", (x1, y), (x2, y), arrowprops=dict(arrowstyle="<->", color=color, lw=1.2, shrinkA=0, shrinkB=0))
    (tag if tagged else lambda *a, **k: ax.text(*a, **k))(ax, (x1 + x2) / 2, y + off, text, color, fs=fs, va="bottom") if tagged else \
        ax.text((x1 + x2) / 2, y + off, text, ha="center", va="bottom", color=color, fontsize=fs)


def dim_v(ax, x, y1, y2, text, color=RED, off=0.06, fs=11, ha="left", tagged=True):
    ax.annotate("", (x, y1), (x, y2), arrowprops=dict(arrowstyle="<->", color=color, lw=1.2, shrinkA=0, shrinkB=0))
    if tagged:
        tag(ax, x + off, (y1 + y2) / 2, text, color, fs=fs, ha=ha)
    else:
        ax.text(x + off, (y1 + y2) / 2, text, ha=ha, va="center", color=color, fontsize=fs)


def arc_tangent_profile(p, H, R, n=600):
    f = lambda a: (p / 2 - 2 * R * np.sin(a)) / np.cos(a) - (H - 2 * R * (1 - np.cos(a))) / np.sin(a)
    alpha = brentq(f, 1e-3, np.pi / 2 - 1e-3)
    TL = (p / 2 - 2 * R * np.sin(alpha)) / np.cos(alpha)
    th = np.linspace(0, alpha, n // 4)
    xa, za = R * np.sin(th), H / 2 - R + R * np.cos(th)
    s = np.linspace(0, TL, n // 4)
    xf, zf = xa[-1] + s * np.cos(alpha), za[-1] - s * np.sin(alpha)
    xq, zq = np.concatenate([xa, xf[1:]]), np.concatenate([za, zf[1:]])
    xh = np.concatenate([xq, p / 2 - xq[::-1][1:]]); zh = np.concatenate([zq, -zq[::-1][1:]])
    return np.concatenate([xh, p - xh[::-1][1:]]), np.concatenate([zh, zh[::-1][1:]]), alpha, TL


def offset_curve(x, z, d):
    tx, tz = np.gradient(x), np.gradient(z); nrm = np.hypot(tx, tz)
    return x - d * tz / nrm, z + d * tx / nrm


# ============================================================================= Fig 1: section
def fig_section(path):
    fig, ax = plt.subplots(figsize=(11, 6.4))
    p, H, R, t = 2.0, 1.0, 0.36, 0.13
    x1, z1, alpha, TL = arc_tangent_profile(p, H, R)
    X = np.concatenate([x1 - p, x1[1:], x1[1:] + p]); Z = np.concatenate([z1, z1[1:], z1[1:]])
    xo, zo = offset_curve(X, Z, t / 2); xi, zi = offset_curve(X, Z, -t / 2)
    ax.fill(np.concatenate([xo, xi[::-1]]), np.concatenate([zo, zi[::-1]]), color="#d6e4f0", ec="none")
    ax.plot(xo, zo, "k", lw=1.2); ax.plot(xi, zi, "k", lw=1.2)
    ax.plot(X, Z, color=BLUE, lw=1.0, ls="--")
    ax.text(-p - 0.25, -H / 2 - 0.22, "dashed = mid-surface", color=BLUE, fontsize=8, va="top", ha="left")
    ax.axhline(0, color=GREY, lw=0.6, ls=":")
    # p (crest to crest, normal to ridge)
    yp = H / 2 + t / 2 + 0.40
    ax.annotate("", (0, yp), (p, yp), arrowprops=dict(arrowstyle="<->", color=RED, lw=1.2, shrinkA=0, shrinkB=0))
    tag(ax, p / 2, yp + 0.07, r"$p$", RED, fs=13, va="bottom")
    for xx in (0, p):
        ax.plot([xx, xx], [H / 2 + t / 2, yp], color=RED, lw=0.6, ls=":")
    # H (mid-surface) on the left, H_o (outer) on the right
    xd = -p - 0.35
    ax.annotate("", (xd, -H / 2), (xd, H / 2), arrowprops=dict(arrowstyle="<->", color=RED, lw=1.2, shrinkA=0, shrinkB=0))
    tag(ax, xd - 0.08, 0.0, r"$H$", RED, fs=13, ha="right")
    ax.text(xd - 0.08, -0.32, "(mid)", color=RED, fontsize=8.5, ha="right")
    ax.plot([xd, -p], [H / 2, H / 2], color=RED, lw=0.6, ls=":"); ax.plot([xd, -p / 2], [-H / 2, -H / 2], color=RED, lw=0.6, ls=":")
    xd2 = 2 * p + 0.35
    ax.annotate("", (xd2, -H / 2 - t / 2), (xd2, H / 2 + t / 2), arrowprops=dict(arrowstyle="<->", color=RED, lw=1.2, shrinkA=0, shrinkB=0))
    tag(ax, xd2 + 0.08, 0.0, r"$H_o$", RED, fs=13, ha="left")
    ax.text(xd2 + 0.08, -0.32, "(outer)\n$= H + t$", color=RED, fontsize=8.5, ha="left", va="top")
    ax.plot([2 * p, xd2], [H / 2 + t / 2, H / 2 + t / 2], color=RED, lw=0.6, ls=":"); ax.plot([1.5 * p, xd2], [-H / 2 - t / 2, -H / 2 - t / 2], color=RED, lw=0.6, ls=":")
    # t at a flank of the third period
    xm = p / 4 + p; nx, nz = np.sin(alpha), np.cos(alpha)
    ax.annotate("", (xm + nx * t / 2, nz * t / 2), (xm - nx * t / 2, -nz * t / 2), arrowprops=dict(arrowstyle="<->", color=RED, lw=1.2, shrinkA=0, shrinkB=0))
    ax.plot([xm + nx * t / 2, xm + 0.42], [nz * t / 2, 0.30], color=RED, lw=0.6)
    tag(ax, xm + 0.45, 0.30, r"$t$", RED, fs=13, ha="left")
    # R_c (crest of period 2), R_v (valley of period 1)
    cc = (0, H / 2 - R)
    ax.plot([cc[0], cc[0] + R * np.sin(0.7)], [cc[1], cc[1] + R * np.cos(0.7)], color=GREEN, lw=1.2); ax.plot(*cc, "o", color=GREEN, ms=3)
    tag(ax, cc[0] + 0.30, cc[1] - 0.20, r"$R_c$", GREEN, fs=12, ha="left")
    cv = (-p / 2, -H / 2 + R)
    ax.plot([cv[0], cv[0] - R * np.sin(0.7)], [cv[1], cv[1] - R * np.cos(0.7)], color=GREEN, lw=1.2); ax.plot(*cv, "o", color=GREEN, ms=3)
    tag(ax, cv[0] - 0.30, cv[1] + 0.22, r"$R_v$", GREEN, fs=12, ha="right")
    # alpha at the flank of period 3 (right side), T_L at flank of period 2
    xi3 = p / 4 + p + p / 2   # descending flank of period 3? use inflection at x = 3p/4 + p (ascending)
    xi0 = p / 4 + p * 0.0 + p / 2 + p / 4   # x = p (crest) ... choose inflection at x = 5p/4 (descending flank of period 2)
    xi0 = 5 * p / 4
    ax.plot([xi0 - 0.45, xi0 + 0.55], [0, 0], color=GREEN, lw=0.8, ls=":")
    ax.add_patch(Arc((xi0, 0), 0.7, 0.7, angle=0, theta1=360 - np.degrees(alpha), theta2=360, color=GREEN, lw=1.3))
    tag(ax, xi0 + 0.62, -0.30, r"$\alpha$", GREEN, fs=13, ha="left")
    xs0, zs0 = R * np.sin(alpha) - p, H / 2 - R + R * np.cos(alpha)      # flank of period 1
    xe0, ze0 = xs0 + TL * np.cos(alpha), zs0 - TL * np.sin(alpha)
    off = 0.24
    ax.annotate("", (xs0 - off * nx, zs0 - off * nz), (xe0 - off * nx, ze0 - off * nz), arrowprops=dict(arrowstyle="<->", color=GREEN, lw=1.1, shrinkA=0, shrinkB=0))
    ax.text((xs0 + xe0) / 2 - 0.50, (zs0 + ze0) / 2 - 0.05, r"$T_L$", color=GREEN, fontsize=11, ha="right")
    # crest / valley, faces at crest of period 3
    ax.text(-p, H / 2 + t / 2 + 0.08, "crest", fontsize=9.5, ha="center", color=GREY, fontweight="bold")
    ax.text(p / 2, -H / 2 - t / 2 - 0.10, "valley", fontsize=9.5, ha="center", color=GREY, fontweight="bold", va="top")
    ax.annotate("+z face", (2 * p, H / 2 + t / 2), (2 * p - 0.9, H / 2 + 0.55), fontsize=9, color=PURPLE, arrowprops=dict(arrowstyle="->", color=PURPLE, lw=0.9))
    ax.annotate("−z face", (2 * p, H / 2 - t / 2), (2 * p - 0.75, H / 2 - 0.62), fontsize=9, color=PURPLE, arrowprops=dict(arrowstyle="->", color=PURPLE, lw=0.9), va="top")
    # local axes (bottom-left)
    ox, oy = -p - 0.9, -H / 2 - 0.75
    ax.annotate("", (ox + 0.9, oy), (ox, oy), arrowprops=dict(arrowstyle="->", color="k", lw=1.5))
    ax.annotate("", (ox, oy + 0.75), (ox, oy), arrowprops=dict(arrowstyle="->", color="k", lw=1.5))
    ax.text(ox + 0.95, oy, r"$x$ ($\perp$ ridge)", fontsize=10, va="center")
    ax.text(ox + 0.05, oy + 0.75, r"$z$", fontsize=11, va="bottom")
    ax.text(ox, oy - 0.18, r"$y$ = along ridge (into page)", fontsize=8.5, color=GREY, va="top")
    # notes
    ax.text(-p - 1.25, H / 2 + 1.30, "Section input (Input sheet)  —  section taken perpendicular to the ridges", fontsize=11.5, fontweight="bold")
    ax.text(-p - 1.25, H / 2 + 1.05, r"mode A: give $p, H, t, R_c$ (and $R_v$ if different) $\rightarrow$ $\alpha, T_L$ computed.     mode B: give $p, H, t, \alpha$ $\rightarrow$ $R_c = R_v$ computed.", fontsize=9.5, color=GREY)
    ax.text(-p - 1.25, H / 2 + 0.85, r"$R_c, R_v$ are mid-surface radii (outer radius $- t/2$).   dims_reference: mid $\rightarrow$ enter $H$;  outer $\rightarrow$ enter $H_o$ (sheet uses $H = H_o - t$).", fontsize=9.5, color=GREY)
    bx, by = p / 2 + 0.2, -H / 2 - 0.55
    ax.add_patch(FancyBboxPatch((bx, by - 0.66), 3.55, 0.66, boxstyle="round,pad=0.04", fc="#fdfefe", ec=GREY, lw=0.8))
    ax.text(bx + 0.1, by - 0.08, "Thickness used in the calculation", fontsize=9, fontweight="bold", va="top")
    ax.text(bx + 0.1, by - 0.30, r"$t_{calc}$ = ( $t_{min}$ if given, else $t$ ) $-$ CA", fontsize=9.5, va="top")
    ax.text(bx + 0.1, by - 0.50, r"$t$ nominal;  $t_{min}$ min. after forming (ridge thinning);  CA corrosion allowance", fontsize=8, color=GREY, va="top")
    ax.set_xlim(-p - 1.3, 2 * p + 1.4); ax.set_ylim(-H / 2 - 1.35, H / 2 + 1.45); ax.set_aspect("equal"); ax.axis("off")
    fig.tight_layout(pad=0.2); fig.savefig(path, dpi=130); plt.close(fig)


# ============================================================================= Fig 2: plan / chevron
def fig_plan(path):
    fig, ax = plt.subplots(figsize=(11, 6.0))
    beta = np.radians(60); cot = 1 / np.tan(beta)
    W, Hh = 6.0, 3.4; xc = W / 2
    ax.add_patch(Polygon([[0, 0], [W, 0], [W, Hh], [0, Hh]], closed=True, fc="#f4f6f7", ec="k", lw=1.4))
    xs = np.linspace(0, W, 400)
    for k in np.arange(-W / 2 * cot - 0.3, Hh + 0.3, 0.32):
        ys = k + np.abs(xs - xc) * cot; m = (ys >= 0) & (ys <= Hh)
        ax.plot(np.where(m, xs, np.nan), np.where(m, ys, np.nan), color=BLUE, lw=0.9)
    ax.plot([xc, xc], [0, Hh], color=GREY, ls="-.", lw=1.0)
    ax.text(xc + 0.06, 0.10, "apex line", color=GREY, fontsize=9, bbox=dict(fc="white", ec="none", pad=1))
    tag(ax, W * 0.25, Hh + 0.18, r"zone = left  ($\theta=+\beta$)", BLUE, fs=10)
    tag(ax, W * 0.75, Hh + 0.18, r"zone = right  ($\theta=-\beta$)", BLUE, fs=10)
    # flow axis
    ax.annotate("", (xc, Hh + 0.75), (xc, Hh + 0.02), arrowprops=dict(arrowstyle="->", color="k", lw=1.6))
    ax.text(xc + 0.1, Hh + 0.62, r"flow axis $y_p$ (port to port)", fontsize=9.5)

    def local_axes(px, py, sgn):
        ey = np.array([sgn * np.sin(beta), np.cos(beta)]); ex = np.array([sgn * np.cos(beta), -np.sin(beta)])
        ax.plot([px, px], [py, py + 1.1], color="k", lw=1.0, ls=":"); ax.text(px + 0.04, py + 1.12, r"$y_p$", fontsize=9, va="bottom")
        ax.annotate("", (px + 1.0 * ey[0], py + 1.0 * ey[1]), (px, py), arrowprops=dict(arrowstyle="->", color=GREEN, lw=2.2))
        ax.annotate("", (px + 0.8 * ex[0], py + 0.8 * ex[1]), (px, py), arrowprops=dict(arrowstyle="->", color=RED, lw=2.2))
        ax.text(px + 1.1 * ey[0], py + 1.1 * ey[1], r"$y$ ($\parallel$ ridge)", color=GREEN, fontsize=9.5, ha="left" if sgn > 0 else "right", va="center", bbox=dict(fc="white", ec="none", pad=0.5))
        ax.text(px + 0.95 * ex[0], py + 0.95 * ex[1], r"$x$ ($\perp$ ridge)", color=RED, fontsize=9.5, ha="left" if sgn > 0 else "right", va="top", bbox=dict(fc="white", ec="none", pad=0.5))
        th = np.degrees(beta); t1, t2 = (90 - th, 90) if sgn > 0 else (90, 90 + th)
        ax.add_patch(Arc((px, py), 1.3, 1.3, angle=0, theta1=t1, theta2=t2, color=GREEN, lw=1.5))
        ax.plot(px, py, "ko", ms=3)
        return ex, ey

    local_axes(xc - 1.7, 1.1, -1)
    ex, ey = local_axes(xc + 1.7, 1.1, +1)
    tag(ax, xc + 1.7 + 0.55, 1.1 + 0.82, r"$\beta$  (beta_reference = flow)", GREEN, fs=10.5, ha="left")
    # beta from horizontal on the right-zone axes
    px, py = xc + 1.7, 1.1
    ax.plot([px, px + 1.15], [py, py], color=ORANGE, lw=1.0, ls=":")
    ax.add_patch(Arc((px, py), 1.0, 1.0, angle=0, theta1=0, theta2=90 - np.degrees(beta), color=ORANGE, lw=1.5))
    tag(ax, px + 1.25, py + 0.05, r"$90°-\beta$  (beta_reference = horizontal)", ORANGE, fs=10, ha="left")
    # pitch normal vs along x_p (bottom-left zone)
    qx, qy = 0.55, 0.55
    e_n = np.array([-np.cos(beta), -np.sin(beta)]) * -1   # normal to left-zone ridge pointing right-up
    # left zone ridges: y = k + (xc - x) cot ; direction (1, -cot) normalized ; normal (cot, 1)/n
    dvec = np.array([1, -cot]); dvec /= np.linalg.norm(dvec); nvec = np.array([cot, 1]); nvec /= np.linalg.norm(nvec)
    spacing = 0.32 * np.sin(beta)          # normal distance between ridge lines (k step 0.32 along y_p)
    # normal pitch arrow between two adjacent ridge lines
    ax.annotate("", (qx + 2 * spacing * nvec[0], qy + 2 * spacing * nvec[1]), (qx, qy), arrowprops=dict(arrowstyle="<->", color=RED, lw=1.6, shrinkA=0, shrinkB=0))
    tag(ax, qx + 2 * spacing * nvec[0] + 0.12, qy + 2 * spacing * nvec[1] + 0.05, r"$p$  (pitch_reference = normal)", RED, fs=10, ha="left")
    # along x_p
    dx_axis = 2 * spacing / np.cos(beta)   # horizontal distance between the same two ridges... along x_p: p / sin(beta)? ridge angle to x_p is 90-beta, so p_axis = p / sin(90-beta) = p/cos(beta)
    ax.annotate("", (qx + dx_axis, qy), (qx, qy), arrowprops=dict(arrowstyle="<->", color=ORANGE, lw=1.6, shrinkA=0, shrinkB=0))
    tag(ax, qx + dx_axis + 0.12, qy - 0.22, r"$p_{axis} = p/\cos\beta$  (pitch_reference = axis)", ORANGE, fs=10, ha="left")
    ax.plot(qx, qy, "ko", ms=3)
    # plate axes
    ax.annotate("", (-0.6, 1.0), (-0.6, -0.45), arrowprops=dict(arrowstyle="->", color="k", lw=1.5))
    ax.annotate("", (0.9, -0.45), (-0.6, -0.45), arrowprops=dict(arrowstyle="->", color="k", lw=1.5))
    ax.text(-0.57, 1.02, r"$y_p$", fontsize=11); ax.text(0.95, -0.45, r"$x_p$", fontsize=11, va="center")
    ax.text(-0.6, Hh + 1.05, "Chevron / plate-axis input (Input sheet)", fontsize=11.5, fontweight="bold")
    ax.text(-0.6, -0.85, r"$\beta$: angle between ridge ($y$) and flow axis ($y_p$). Drawing angle from horizontal $x_p$ $\Rightarrow$ set beta_reference = horizontal (sheet uses $90°-$input).", fontsize=9, color=GREY)
    ax.text(-0.6, -1.08, r"Pitch measured along $x_p$ on the plan drawing $\Rightarrow$ pitch_reference = axis (sheet converts $p = p_{axis}\cos\beta$). Section A–A of Fig. 1 is taken along $x$ ($\perp$ ridge).", fontsize=9, color=GREY)
    ax.text(-0.6, -1.31, r"Load case zone: right ($\theta=-\beta$) / left ($\theta=+\beta$) selects which half the plate-axis resultants are rotated into.", fontsize=9, color=GREY)
    ax.set_xlim(-0.8, W + 3.2); ax.set_ylim(-1.5, Hh + 1.25); ax.set_aspect("equal"); ax.axis("off")
    fig.tight_layout(pad=0.2); fig.savefig(path, dpi=130); plt.close(fig)


# ============================================================================= Fig 3: loads and stress points
def fig_loads(path):
    """국부축 판 요소의 합력 6성분."""
    fig, ax = plt.subplots(figsize=(11, 5.4))
    a = 2.0
    ax.add_patch(Polygon([[0, 0], [a, 0], [a, a], [0, a]], closed=True, fc="#eef3f8", ec="k", lw=1.3))
    for k in np.linspace(0.2, a - 0.2, 6):
        ax.plot([k, k], [0.05, a - 0.05], color=BLUE, lw=0.8, alpha=0.7)
    ax.text(a / 2, a / 2, "ridges ∥ y", color=BLUE, fontsize=9, ha="center", bbox=dict(fc="white", ec="none", pad=0.5))
    ax.annotate("", (a + 1.3, 0), (a, 0), arrowprops=dict(arrowstyle="->", color="k", lw=1.4)); ax.text(a + 1.35, 0, r"$x$ ($\perp$ ridge)", va="center", fontsize=10.5)
    ax.annotate("", (0, a + 1.1), (0, a), arrowprops=dict(arrowstyle="->", color="k", lw=1.4)); ax.text(0.06, a + 1.1, r"$y$ ($\parallel$ ridge)", fontsize=10.5, va="bottom")
    # membrane resultants
    ax.annotate("", (a + 0.6, a / 2 + 0.5), (a, a / 2 + 0.5), arrowprops=dict(arrowstyle="->", color=RED, lw=2)); tag(ax, a + 0.68, a / 2 + 0.5, r"$N_x$", RED, fs=12, ha="left")
    ax.annotate("", (a / 2 - 0.5, a + 0.55), (a / 2 - 0.5, a), arrowprops=dict(arrowstyle="->", color=RED, lw=2)); tag(ax, a / 2 - 0.5, a + 0.68, r"$N_y$", RED, fs=12, va="bottom")
    ax.annotate("", (a, a / 2 - 0.6), (a, a / 2 - 0.05), arrowprops=dict(arrowstyle="->", color=ORANGE, lw=2)); tag(ax, a + 0.14, a / 2 - 0.62, r"$N_{xy}$", ORANGE, fs=12, ha="left")
    ax.annotate("", (a / 2 + 0.6, a), (a / 2 + 0.05, a), arrowprops=dict(arrowstyle="->", color=ORANGE, lw=2))
    # moment resultants (double-headed vectors along edges)
    ax.annotate("", (-0.5, a / 2 + 0.6), (-0.5, a / 2 - 0.6), arrowprops=dict(arrowstyle="-|>,head_width=0.3,head_length=0.6", color=GREEN, lw=2.0))
    ax.plot([-0.5, -0.5], [a / 2 + 0.35, a / 2 + 0.6], color=GREEN, lw=2.0)
    tag(ax, -0.62, a / 2, r"$M_x$", GREEN, fs=12, ha="right")
    ax.annotate("", (a / 2 + 0.6, -0.5), (a / 2 - 0.6, -0.5), arrowprops=dict(arrowstyle="-|>,head_width=0.3,head_length=0.6", color=GREEN, lw=2.0))
    ax.plot([a / 2 + 0.35, a / 2 + 0.6], [-0.5, -0.5], color=GREEN, lw=2.0)
    tag(ax, a / 2, -0.75, r"$M_y$", GREEN, fs=12, va="top")
    tag(ax, a + 0.3, -0.75, r"$M_{xy}$", PURPLE, fs=12, va="top", ha="left")
    # explanatory column on the right
    x0 = a + 2.6
    ax.text(x0, a + 1.25, "Load case input (Input sheet)", fontsize=12, fontweight="bold")
    ax.text(x0, a + 0.95, r"$N_x, N_y, N_{xy}$ [N/mm]  membrane forces per unit width", fontsize=9.5)
    ax.text(x0, a + 0.70, r"$M_x, M_y, M_{xy}$ [N·mm/mm]  moments per unit width", fontsize=9.5)
    ax.text(x0, a + 0.30, r"axes = local : values in ridge axes $(x, y)$ as drawn.", fontsize=9.5, color=GREY)
    ax.text(x0, a + 0.05, r"axes = plate : values in plate axes $(x_p, y_p)$; the sheet", fontsize=9.5, color=GREY)
    ax.text(x0, a - 0.20, r"                    rotates them by $\theta = -\beta$ (zone right) / $+\beta$ (left).", fontsize=9.5, color=GREY)
    ax.text(x0, a - 0.65, r"$M_x$ : bending about $y$ $\rightarrow$ curvature $\kappa_x$ across the ridges (flexible).", fontsize=9.5, color=GREY)
    ax.text(x0, a - 0.90, r"$M_y$ : bending along the ridges $\rightarrow$ $\kappa_y$ (stiff direction).", fontsize=9.5, color=GREY)
    ax.text(x0, a - 1.15, r"$M_{xy}$ : twisting moment, $M_{xy} = D_{66}\cdot 2\kappa_{xy}$.", fontsize=9.5, color=GREY)
    ax.text(x0, a - 1.60, r"Sign: tension positive.  $M > 0$ puts the $+z$ face in compression", fontsize=9.5, color=GREY)
    ax.text(x0, a - 1.85, r"($\kappa = +w_{,\alpha\beta}$ convention, surface strain $e = \varepsilon - \zeta\kappa$).", fontsize=9.5, color=GREY)
    ax.text(x0, a - 2.30, "Resultants are NOT computed here: obtain them from a closed-form plate", fontsize=9, color=GREY)
    ax.text(x0, a - 2.52, "solution or the equivalent-plate FE model and enter them in the yellow cells.", fontsize=9, color=GREY)
    ax.set_xlim(-1.2, a + 9.3); ax.set_ylim(-1.35, a + 1.6); ax.set_aspect("equal"); ax.axis("off")
    fig.tight_layout(pad=0.2); fig.savefig(path, dpi=130); plt.close(fig)


def fig_stress_points(path):
    """응력 출력 위치: 산·골 × +z/−z 면."""
    fig, ax = plt.subplots(figsize=(11, 4.6))
    p, H, R, t = 2.0, 1.0, 0.36, 0.16
    x1, z1, alpha, TL = arc_tangent_profile(p, H, R)
    X = np.concatenate([x1 - p, x1[1:]]); Z = np.concatenate([z1, z1[1:]])          # crest at x=0, valleys at ±p/2
    m = (X >= -p / 2 - 0.01) & (X <= p + 0.01); X, Z = X[m], Z[m]
    xo, zo = offset_curve(X, Z, t / 2); xi, zi = offset_curve(X, Z, -t / 2)
    ax.fill(np.concatenate([xo, xi[::-1]]), np.concatenate([zo, zi[::-1]]), color="#d6e4f0", ec="none")
    ax.plot(xo, zo, "k", lw=1.2); ax.plot(xi, zi, "k", lw=1.2); ax.plot(X, Z, color=BLUE, lw=0.8, ls="--")
    pts = [(0, H / 2 + t / 2, "crest, +z face  (Stress sheet col. B)", PURPLE, (0.45, 0.42), "left"),
           (0, H / 2 - t / 2, "crest, −z face  (col. C)", PURPLE, (-0.3, -1.0), "right"),
           (p / 2, -H / 2 + t / 2, "valley, +z face  (col. D)", ORANGE, (0.55, 0.45), "left"),
           (p / 2, -H / 2 - t / 2, "valley, −z face  (col. E)", ORANGE, (0.55, -0.30), "left")]
    for x, z, lab, col, (dx, dz), ha in pts:
        ax.plot(x, z, "o", color=col, ms=8, mec="k")
        ax.annotate(lab, (x, z), (x + dx, z + dz), fontsize=10, color=col, fontweight="bold", arrowprops=dict(arrowstyle="->", color=col, lw=0.9), va="center", ha=ha)
    ax.annotate("", (-p / 2 - 1.5, 0.65), (-p / 2 - 1.5, -0.25), arrowprops=dict(arrowstyle="->", color="k", lw=1.4)); ax.text(-p / 2 - 1.45, 0.67, r"$+z$", fontsize=11, va="bottom")
    ax.text(-p / 2 - 1.5, -1.0, "+z / −z: the two surfaces of the sheet (+z = crest side of the plate)", fontsize=8.5, color=GREY, va="top")
    x0 = p + 2.2
    ax.text(x0, H / 2 + 0.55, "Stress output points (Stress sheet)", fontsize=12, fontweight="bold")
    ax.text(x0, H / 2 + 0.22, r"$\sigma_s$ (across ridge), $\sigma_y$ (along ridge), $\tau_{sy}$, von Mises", fontsize=9.5, color=GREY)
    ax.text(x0, H / 2 - 0.03, "at crest and valley on both surfaces = 4 points;", fontsize=9.5, color=GREY)
    ax.text(x0, H / 2 - 0.28, "the maximum von Mises and its location are reported.", fontsize=9.5, color=GREY)
    ax.text(x0, H / 2 - 0.70, r"For $N_x > 0$: crest −z face and valley +z face carry", fontsize=9.5, color=GREY)
    ax.text(x0, H / 2 - 0.95, r"the stress concentration $K_t = 1 + 6f/t$ (Briassoulis).", fontsize=9.5, color=GREY)
    ax.text(x0, H / 2 - 1.37, "Flank (mid-height) stresses are available only in the Python module.", fontsize=9, color=GREY)
    ax.set_xlim(-p / 2 - 1.6, p + 7.6); ax.set_ylim(-H / 2 - 1.15, H / 2 + 0.85); ax.set_aspect("equal"); ax.axis("off")
    fig.tight_layout(pad=0.2); fig.savefig(path, dpi=130); plt.close(fig)


if __name__ == "__main__":
    fig_section(os.path.join(HERE, "xl_fig1_section.png"))
    fig_plan(os.path.join(HERE, "xl_fig2_plan.png"))
    fig_loads(os.path.join(HERE, "xl_fig3_loads.png"))
    fig_stress_points(os.path.join(HERE, "xl_fig4_stress_points.png"))
    print("saved 4 figures in", HERE)
