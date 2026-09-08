"""쉐브론 전열판 치수 정의 그림 생성 (docs/fig/plate_dimensions.png / .svg).

패널 (a) 평면도, (b) 능선 직각 단면 A-A, (c) 좌표계·쉐브론 각 정의, (d) 적층 단면, (e) 테두리·apex 상세.
기호 설명(한글)은 docs/04_plate_dimension_form.md 참조.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Arc, FancyArrowPatch, Polygon
from scipy.optimize import brentq

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "mathtext.fontset": "dejavusans"})
BLUE, RED, GREY, GREEN = "#1f4e79", "#c0392b", "#7f8c8d", "#1e8449"


def dim_h(ax, x1, x2, y, text, color=RED, off=0.0, fs=11):
    ax.annotate("", (x1, y), (x2, y), arrowprops=dict(arrowstyle="<->", color=color, lw=1.2, shrinkA=0, shrinkB=0))
    ax.text((x1 + x2) / 2, y + off, text, ha="center", va="bottom", color=color, fontsize=fs)


def dim_v(ax, x, y1, y2, text, color=RED, off=0.0, fs=11, ha="left"):
    ax.annotate("", (x, y1), (x, y2), arrowprops=dict(arrowstyle="<->", color=color, lw=1.2, shrinkA=0, shrinkB=0))
    ax.text(x + off, (y1 + y2) / 2, text, ha=ha, va="center", color=color, fontsize=fs)


# ----------------------------------------------------------------------------- (a) plan view
def panel_plan(ax):
    Rp, Re, rport, eport, beta = 1.0, 0.82, 0.13, 0.55, np.radians(60)
    ax.add_patch(Circle((0, 0), Rp, fc="#f4f6f7", ec="k", lw=1.8))
    ax.add_patch(Circle((0, 0), Re, fc="none", ec=GREY, lw=1.2, ls="--"))
    # chevron ridges  y = y0 + |x| cot(beta)   (V opening upward, apex line x=0)
    cot = 1 / np.tan(beta)
    xs = np.linspace(-Re, Re, 400)
    clip = Circle((0, 0), Re, transform=ax.transData)
    for y0 in np.arange(-2.2, 1.2, 0.11):
        ys = y0 + np.abs(xs) * cot
        m = xs**2 + ys**2 <= Re**2
        ln, = ax.plot(np.where(m, xs, np.nan), np.where(m, ys, np.nan), color=BLUE, lw=0.9)
        ln.set_clip_path(clip)
    for sy in (+1, -1):  # ports
        ax.add_patch(Circle((0, sy * eport), rport, fc="white", ec="k", lw=1.5))
    ax.plot([0, 0], [-Rp, Rp], color=GREY, lw=0.8, ls="-.")
    ax.text(0.03, -0.93, "apex line", color=GREY, fontsize=9, ha="left", bbox=dict(fc="white", ec="none", pad=1))
    # beta annotation (right half ridge vs. flow axis y_p)
    px = 0.45
    py = -0.15 + px * cot
    ax.plot([px, px], [py, py + 0.42], color=GREEN, lw=1.0, ls=":")
    L = 0.42
    ax.plot([px, px + L * np.sin(beta)], [py, py + L * np.cos(beta)], color=GREEN, lw=1.6)
    ax.add_patch(Arc((px, py), 0.5, 0.5, angle=0, theta1=90 - np.degrees(beta), theta2=90, color=GREEN, lw=1.4))
    ax.text(px + 0.09, py + 0.30, r"$\beta$", color=GREEN, fontsize=14)
    # dimensions
    dim_h(ax, -Rp, Rp, -1.14, r"$D_p$ (plate outer dia.)", off=-0.14)
    ax.plot([-Rp, -Rp], [-Rp, -1.14], color=RED, lw=0.6, ls=":"); ax.plot([Rp, Rp], [-Rp, -1.14], color=RED, lw=0.6, ls=":")
    dim_h(ax, -Re, Re, 1.08, r"$D_e$ (effective corrugated dia.)", off=0.02)
    ax.plot([-Re, -Re], [Re * 0.57, 1.08], color=RED, lw=0.6, ls=":"); ax.plot([Re, Re], [Re * 0.57, 1.08], color=RED, lw=0.6, ls=":")
    dim_h(ax, Re, Rp, 0.0, "", color=RED)
    ax.text(Rp + 0.03, 0.0, r"$b_m$", color=RED, fontsize=12, va="center")
    dim_v(ax, -1.18, 0, eport, r"$e_{port}$", off=-0.04, ha="right")
    ax.plot([-1.18, 0], [0, 0], color=RED, lw=0.6, ls=":"); ax.plot([-1.18, 0], [eport, eport], color=RED, lw=0.6, ls=":")
    dim_h(ax, -rport, rport, eport + rport + 0.03, r"$d_{port}$", off=0.01)
    # plate axes
    ax.annotate("", (-1.45, -1.0), (-1.45, -1.55), arrowprops=dict(arrowstyle="->", color="k", lw=1.4))
    ax.annotate("", (-0.9, -1.55), (-1.45, -1.55), arrowprops=dict(arrowstyle="->", color="k", lw=1.4))
    ax.text(-1.42, -0.98, r"$y_p$ (flow / port axis)", fontsize=10)
    ax.text(-0.88, -1.57, r"$x_p$", fontsize=10, va="top")
    # section cut A-A perpendicular to ridge on right half
    cx, cy = 0.42, 0.42 * cot - 0.7
    d = np.array([np.cos(beta), -np.sin(beta)])  # perpendicular to ridge dir (sinβ, cosβ)
    a0, a1 = np.array([cx, cy]) - 0.32 * d, np.array([cx, cy]) + 0.32 * d
    ax.plot([a0[0], a1[0]], [a0[1], a1[1]], color=RED, lw=2.0)
    nrm = np.array([np.sin(beta), np.cos(beta)])  # ridge direction (perpendicular to cut)
    for pnt in (a0, a1):
        ax.text(pnt[0] + 0.12 * nrm[0], pnt[1] + 0.12 * nrm[1], "A", color=RED, fontsize=12, fontweight="bold", ha="center", va="center",
                bbox=dict(fc="white", ec="none", pad=1))
    ax.set_xlim(-1.65, 1.6); ax.set_ylim(-1.8, 1.3); ax.set_aspect("equal"); ax.axis("off")
    ax.set_title("(a) Plan view  —  chevron plate (circular SPHX plate shown)", loc="left", fontsize=12, fontweight="bold")
    ax.text(-0.6, -1.72, "Rectangular plate: use $L_p\\times W_p$ and $L_e\\times W_e$ instead of $D_p$, $D_e$.  $b_m$ = weld/gasket margin.", fontsize=8.5, color=GREY, va="top")


# ----------------------------------------------------------------------------- (b) section A-A
def arc_tangent_profile(p, H, R, n=600):
    """대칭 원호-접선 단면. crest at x=0 (z=+H/2), valley at x=p/2 (z=-H/2). returns x, z (mid-surface), alpha, TL."""
    f = lambda a: (p / 2 - 2 * R * np.sin(a)) / np.cos(a) - (H - 2 * R * (1 - np.cos(a))) / np.sin(a)
    alpha = brentq(f, 1e-3, np.pi / 2 - 1e-3)
    TL = (p / 2 - 2 * R * np.sin(alpha)) / np.cos(alpha)
    # quarter: crest arc (0..alpha) + flank
    th = np.linspace(0, alpha, n // 4)
    xa, za = R * np.sin(th), H / 2 - R + R * np.cos(th)
    s = np.linspace(0, TL, n // 4)
    xf, zf = xa[-1] + s * np.cos(alpha), za[-1] - s * np.sin(alpha)
    xq, zq = np.concatenate([xa, xf[1:]]), np.concatenate([za, zf[1:]])
    # mirror about x=p/4, z=0 (point symmetry) to reach valley
    xh = np.concatenate([xq, p / 2 - xq[::-1][1:]])
    zh = np.concatenate([zq, -zq[::-1][1:]])
    # mirror about x=p/2 to full period
    xp_ = np.concatenate([xh, p - xh[::-1][1:]])
    zp_ = np.concatenate([zh, zh[::-1][1:]])
    return xp_, zp_, alpha, TL


def offset_curve(x, z, d):
    tx, tz = np.gradient(x), np.gradient(z)
    nrm = np.hypot(tx, tz)
    nx, nz = -tz / nrm, tx / nrm
    return x + d * nx, z + d * nz


def panel_section(ax):
    p, H, R, t = 2.0, 1.0, 0.36, 0.12
    x1, z1, alpha, TL = arc_tangent_profile(p, H, R)
    X = np.concatenate([x1 - p, x1[1:], x1[1:] + p])
    Z = np.concatenate([z1, z1[1:], z1[1:]])
    xo, zo = offset_curve(X, Z, t / 2)
    xi, zi = offset_curve(X, Z, -t / 2)
    ax.fill(np.concatenate([xo, xi[::-1]]), np.concatenate([zo, zi[::-1]]), color="#d6e4f0", ec="none")
    ax.plot(xo, zo, "k", lw=1.2); ax.plot(xi, zi, "k", lw=1.2)
    ax.plot(X, Z, color=BLUE, lw=1.0, ls="--", label="mid-surface")
    ax.axhline(0, color=GREY, lw=0.6, ls=":")
    # pitch p (crest to crest)
    dim_h(ax, 0, p, H / 2 + t / 2 + 0.42, r"$p$  (pitch, $=2c$)", off=0.03, fs=12)
    ax.plot([0, 0], [H / 2 + t / 2, H / 2 + t / 2 + 0.42], color=RED, lw=0.6, ls=":")
    ax.plot([p, p], [H / 2 + t / 2, H / 2 + t / 2 + 0.42], color=RED, lw=0.6, ls=":")
    # depth H (mid-surface crest to trough)
    xd = -p - 0.55
    dim_v(ax, xd, -H / 2, H / 2, r"$H$ ($=2f$)", off=-0.06, ha="right", fs=12)
    ax.plot([xd, -p], [H / 2, H / 2], color=RED, lw=0.6, ls=":"); ax.plot([xd, -p / 2], [-H / 2, -H / 2], color=RED, lw=0.6, ls=":")
    # overall height H_o
    xd2 = 2 * p + 0.55
    dim_v(ax, xd2, -H / 2 - t / 2, H / 2 + t / 2, r"$H_o=H+t$", off=0.06, fs=11)
    ax.plot([2 * p, xd2], [H / 2 + t / 2, H / 2 + t / 2], color=RED, lw=0.6, ls=":"); ax.plot([1.5 * p, xd2], [-H / 2 - t / 2, -H / 2 - t / 2], color=RED, lw=0.6, ls=":")
    # thickness t at a flank (normal direction)
    xm = p / 4 + p  # inflection of second period
    nx, nz = np.sin(alpha), np.cos(alpha)
    ax.annotate("", (xm + nx * t / 2, nz * t / 2), (xm - nx * t / 2, -nz * t / 2), arrowprops=dict(arrowstyle="<->", color=RED, lw=1.2, shrinkA=0, shrinkB=0))
    ax.text(xm + 0.16, 0.12, r"$t$", color=RED, fontsize=13)
    # crest radius R_c, valley radius R_v
    cc = (0, H / 2 - R)
    ax.plot([cc[0], cc[0] + R * np.sin(0.6)], [cc[1], cc[1] + R * np.cos(0.6)], color=GREEN, lw=1.2)
    ax.plot(*cc, "o", color=GREEN, ms=3)
    ax.text(cc[0] + 0.30, cc[1] + 0.05, r"$R_c$", color=GREEN, fontsize=12, bbox=dict(fc="white", ec="none", pad=1))
    cv = (p / 2, -H / 2 + R)
    ax.plot([cv[0], cv[0] - R * np.sin(0.6)], [cv[1], cv[1] - R * np.cos(0.6)], color=GREEN, lw=1.2)
    ax.plot(*cv, "o", color=GREEN, ms=3)
    ax.text(cv[0] - 0.55, cv[1] - 0.05, r"$R_v$", color=GREEN, fontsize=12, ha="right", bbox=dict(fc="white", ec="none", pad=1))
    # flank angle alpha and tangent length T_L
    xi0 = p / 4 - p  # inflection in first period
    ax.plot([xi0 - 0.35, xi0 + 0.45], [0, 0], color=GREEN, lw=0.8, ls=":")
    ax.add_patch(Arc((xi0, 0), 0.6, 0.6, angle=0, theta1=360 - np.degrees(alpha), theta2=360, color=GREEN, lw=1.3))
    ax.text(xi0 + 0.34, -0.17, r"$\alpha$", color=GREEN, fontsize=13)
    # T_L along flank of first period (from crest arc end to valley arc start)
    xs0, zs0 = R * np.sin(alpha) - p, H / 2 - R + R * np.cos(alpha)
    xe0, ze0 = xs0 + TL * np.cos(alpha), zs0 - TL * np.sin(alpha)
    off = 0.22
    ax.annotate("", (xs0 - off * nx, zs0 - off * nz), (xe0 - off * nx, ze0 - off * nz), arrowprops=dict(arrowstyle="<->", color=GREEN, lw=1.1, shrinkA=0, shrinkB=0))
    ax.text((xs0 + xe0) / 2 - 0.42, (zs0 + ze0) / 2 - 0.30, r"$T_L$", color=GREEN, fontsize=12)
    # labels
    ax.text(-p - 0.05, H / 2 + 0.12, "crest", fontsize=9, ha="center", color=GREY)
    ax.text(-p / 2, -H / 2 - 0.2, "valley", fontsize=9, ha="center", color=GREY)
    ax.annotate("", (-p, H / 2 + t / 2 + 0.05), (-p, H / 2 + t / 2 + 0.05), arrowprops=dict(arrowstyle="-"))
    ax.text(-p - 1.3, -H / 2 - 0.60, "$H$: depth between mid-surface crest and trough (not overall height).  $t$: sheet thickness.  $p$: pitch measured perpendicular to ridges.", fontsize=8.5, color=GREY)
    ax.text(-p - 1.3, -H / 2 - 0.78, "Profile options: (1) sinusoidal → $p,H,t$   (2) arc + tangent → $p,H,t,R_c,R_v,\\alpha$ (or $T_L$)   (3) trapezoidal → $p,H,t,\\alpha$, flat widths", fontsize=8.5, color=GREY)
    ax.text(-p - 1.3, -H / 2 - 0.96, "Local axes: $x$ = across ridges (this section, flexible), $y$ = along ridges (stiff), $z$ = plate normal", fontsize=8.5, color=GREY)
    ax.set_xlim(-p - 1.35, 2 * p + 1.2); ax.set_ylim(-H / 2 - 1.1, H / 2 + 0.95); ax.set_aspect("equal"); ax.axis("off")
    ax.set_title("(b) Section A–A, perpendicular to ridge (arc + tangent profile)", loc="left", fontsize=12, fontweight="bold")
    ax.legend(loc="upper right", fontsize=9, frameon=False)


# ----------------------------------------------------------------------------- (c) axes / chevron angle
def panel_axes(ax):
    beta = np.radians(60)
    cot = 1 / np.tan(beta)
    W, Hh = 6.0, 3.2                      # one patch of plate n (plan view, V-shaped ridges)
    xc = W / 2
    ax.add_patch(Polygon([[0, 0], [W, 0], [W, Hh], [0, Hh]], closed=True, fc="#f4f6f7", ec="k", lw=1.4))
    xs = np.linspace(0, W, 400)
    for k in np.arange(-W / 2 * cot - 0.3, Hh + 0.3, 0.30):
        ys = k + np.abs(xs - xc) * cot
        m = (ys >= 0) & (ys <= Hh)
        ax.plot(np.where(m, xs, np.nan), np.where(m, ys, np.nan), color=BLUE, lw=0.9)
    ax.plot([xc, xc], [0, Hh], color=GREY, ls="-.", lw=0.9)
    ax.text(xc + 0.06, 0.12, "apex line", color=GREY, fontsize=9, bbox=dict(fc="white", ec="none", pad=1))
    ax.text(W * 0.25, Hh + 0.10, r"left half: ridges at $-\beta$", ha="center", color=BLUE, fontsize=10.5)
    ax.text(W * 0.75, Hh + 0.10, r"right half: ridges at $+\beta$", ha="center", color=BLUE, fontsize=10.5)
    ax.text(0.0, Hh + 0.45, "Plate $n$  (plan-view patch)", fontsize=11, fontweight="bold")
    ax.annotate("", (xc, -0.05), (xc, -0.65), arrowprops=dict(arrowstyle="->", color="k", lw=1.6))
    ax.text(xc + 0.12, -0.42, "flow (port to port) = $y_p$", fontsize=9.5)

    def local_axes(px, py, sgn, lab_beta):
        ey = np.array([sgn * np.sin(beta), np.cos(beta)])          # along ridge
        ex = np.array([sgn * np.cos(beta), -np.sin(beta)])         # across ridge (perpendicular), pointing away from apex and down
        ax.plot([px, px], [py, py + 1.05], color="k", lw=1.0, ls=":")       # y_p reference
        ax.text(px + 0.04, py + 1.07, r"$y_p$", fontsize=9, va="bottom", ha="left")
        ax.annotate("", (px + 1.0 * ey[0], py + 1.0 * ey[1]), (px, py), arrowprops=dict(arrowstyle="->", color=GREEN, lw=2.2))
        ax.annotate("", (px + 0.8 * ex[0], py + 0.8 * ex[1]), (px, py), arrowprops=dict(arrowstyle="->", color=RED, lw=2.2))
        ax.text(px + 1.1 * ey[0], py + 1.1 * ey[1], r"$y$  ($\parallel$ ridge)", color=GREEN, fontsize=9.5,
                ha="left" if sgn > 0 else "right", va="center", bbox=dict(fc="white", ec="none", pad=0.5))
        ax.text(px + 0.95 * ex[0], py + 0.95 * ex[1], r"$x$  ($\perp$ ridge)", color=RED, fontsize=9.5,
                ha="left" if sgn > 0 else "right", va="top", bbox=dict(fc="white", ec="none", pad=0.5))
        th = np.degrees(beta)
        t1, t2 = (90 - th, 90) if sgn > 0 else (90, 90 + th)
        ax.add_patch(Arc((px, py), 1.3, 1.3, angle=0, theta1=t1, theta2=t2, color=GREEN, lw=1.5))
        ax.text(px + sgn * 0.42, py + 0.74, lab_beta, color=GREEN, fontsize=13, ha="center")
        ax.plot(px, py, "ko", ms=3)

    local_axes(xc + 1.6, 1.0, +1, r"$+\beta$")
    local_axes(xc - 1.6, 1.0, -1, r"$-\beta$")
    # plate axes (bottom-left)
    ax.annotate("", (-0.6, 0.9), (-0.6, -0.45), arrowprops=dict(arrowstyle="->", color="k", lw=1.5))
    ax.annotate("", (0.9, -0.45), (-0.6, -0.45), arrowprops=dict(arrowstyle="->", color="k", lw=1.5))
    ax.text(-0.57, 0.92, r"$y_p$", fontsize=11); ax.text(0.95, -0.45, r"$x_p$", fontsize=11, va="center")

    # adjacent plate n+1 (rotated 180 deg): inverted V
    x1, w1, h1 = W + 1.1, 2.4, 3.2
    xc1 = x1 + w1 / 2
    ax.add_patch(Polygon([[x1, 0], [x1 + w1, 0], [x1 + w1, h1], [x1, h1]], closed=True, fc="#fdf2f0", ec="k", lw=1.2))
    xs1 = np.linspace(x1, x1 + w1, 200)
    for k in np.arange(-0.3, h1 + w1 / 2 * cot + 0.3, 0.30):
        ys = k - np.abs(xs1 - xc1) * cot
        m = (ys >= 0) & (ys <= h1)
        ax.plot(np.where(m, xs1, np.nan), np.where(m, ys, np.nan), color=RED, lw=0.9, alpha=0.8)
    ax.plot([xc1, xc1], [0, h1], color=GREY, ls="-.", lw=0.9)
    ax.text(x1, h1 + 0.45, "Plate $n+1$", fontsize=11, fontweight="bold")
    ax.text(x1, h1 + 0.10, "(= plate $n$ rotated 180°)", fontsize=9.5, color=GREY)
    ax.annotate("", (x1 - 0.1, h1 / 2), (W + 0.1, h1 / 2), arrowprops=dict(arrowstyle="<->", color=GREY, lw=1.0))
    ax.text((W + x1) / 2, h1 / 2 + 0.08, "stacked", fontsize=8.5, color=GREY, ha="center")

    ax.text(0.0, -0.95, r"$\beta$ = angle between a ridge line ($y$) and the flow axis ($y_p$), measured in the plane of the plate.  "
                        r"If a drawing gives the angle from $x_p$ (horizontal), then $\beta = 90° -$ that angle.", fontsize=9, color=GREY)
    ax.text(0.0, -1.20, r"Equivalent stiffness: compute $\bar A,\bar D$ in the local $(x,y)$ axes of one half, then rotate by $+\beta$ (right half) or $-\beta$ (left half) into plate axes $(x_p,y_p)$.", fontsize=9, color=GREY)
    ax.text(0.0, -1.45, r"Stacking: plate $n+1$ is plate $n$ rotated 180°, so the ridges of the two plates cross at an angle $2\beta$ and touch at crest-to-crest points (cross-corrugated channel).", fontsize=9, color=GREY)
    ax.set_xlim(-0.9, x1 + w1 + 0.3); ax.set_ylim(-1.6, Hh + 0.8); ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(r"(c) Axis convention and chevron angle $\beta$", loc="left", fontsize=12, fontweight="bold")


# ----------------------------------------------------------------------------- (d) stack
def panel_stack(ax):
    p, H, R, t = 2.0, 1.0, 0.36, 0.12
    x1, z1, *_ = arc_tangent_profile(p, H, R)
    X = np.concatenate([x1 - p, x1[1:], x1[1:] + p]); Z = np.concatenate([z1, z1[1:], z1[1:]])
    sp = H + t  # plate pitch (crest-to-crest contact)
    for k, (dz, flip) in enumerate(((2 * sp, 1), (sp, -1), (0, 1))):
        zz = flip * Z + dz
        xo, zo = offset_curve(X, zz, t / 2); xi, zi = offset_curve(X, zz, -t / 2)
        ax.fill(np.concatenate([xo, xi[::-1]]), np.concatenate([zo, zi[::-1]]), color="#d6e4f0" if k != 1 else "#f9e0d9", ec="k", lw=1.0)
        ax.plot(X, zz, color=BLUE if k != 1 else RED, lw=0.8, ls="--")
    # contact points (crest of lower vs valley of upper in this schematic)
    for xc in (-p, 0, p, 2 * p):
        ax.plot(xc, sp + H / 2 + t / 2, "o", color=GREEN, ms=5, zorder=5)
        ax.plot(xc, sp - H / 2 - t / 2, "o", color=GREEN, ms=5, zorder=5)
    ax.text(2 * p + 0.15, sp + H / 2 + t / 2, "contact point", fontsize=9, color=GREEN, va="center")
    xd = -p - 0.5
    dim_v(ax, xd, sp, 2 * sp, r"$s_p$ (plate pitch)", off=-0.06, ha="right", fs=11)
    dim_v(ax, 2 * p + 0.5, sp - H / 2 + t / 2, sp + H / 2 - t / 2, r"$b_{ch}$ (channel gap)", off=0.06, fs=11)
    ax.text(-p - 1.2, -H / 2 - 0.50, "Adjacent plates are rotated 180° so their ridges cross at $2\\beta$; this section is schematic (cut is not perpendicular to both ridge sets).", fontsize=8.5, color=GREY)
    ax.text(-p - 1.2, -H / 2 - 0.68, "Plate count $N_p$ and channel pattern (1-pass / 2-pass) → fill in the form.", fontsize=8.5, color=GREY)
    ax.set_xlim(-p - 1.3, 2 * p + 1.5); ax.set_ylim(-H / 2 - 0.85, 2 * sp + H / 2 + 0.5); ax.set_aspect("equal"); ax.axis("off")
    ax.set_title("(d) Stack section (schematic)  —  plate pitch and channel gap", loc="left", fontsize=12, fontweight="bold")


# ----------------------------------------------------------------------------- (e) edge / apex details
def panel_edge(ax):
    # ---- (e1) rim section: paired plates n (top) / n+1 (bottom) form one channel; rims meet and are welded at the edge
    H, t, lam = 1.0, 0.12, 1.5
    xr1, xr2, xr3 = 3.0, 4.2, 5.8                     # end of full corrugation | end of run-out | plate edge
    zmid = (H + t) / 2 + 0.15                          # corrugation mid-plane height of plate n above channel centre
    zrim = t / 2 + 0.02                                # rim plane (rims touch at the weld)
    def profile(x):
        amp = np.where(x < xr1, 1.0, np.where(x < xr2, 0.5 * (1 + np.cos(np.pi * (x - xr1) / (xr2 - xr1))), 0.0))
        drop = np.where(x < xr1, 0.0, np.where(x < xr2, 0.5 * (1 - np.cos(np.pi * (x - xr1) / (xr2 - xr1))), 1.0))
        return zmid + H / 2 * np.sin(2 * np.pi * x / lam) * amp - (zmid - zrim) * drop
    x = np.linspace(0.0, xr3, 700)
    for sgn, col in ((+1, "#d6e4f0"), (-1, "#f9e0d9")):
        z = sgn * profile(x)
        xo, zo = offset_curve(x, z, t / 2); xi, zi = offset_curve(x, z, -t / 2)
        ax.fill(np.concatenate([xo, xi[::-1]]), np.concatenate([zo, zi[::-1]]), color=col, ec="k", lw=1.0)
    ax.text(0.1, zmid + H / 2 + 0.28, "plate $n$", fontsize=9.5, color=BLUE)
    ax.text(0.1, -zmid - H / 2 - 0.28, "plate $n+1$", fontsize=9.5, color=RED, va="top")
    ax.add_patch(Circle((xr3 + 0.05, 0.0), 0.17, fc="#555555", ec="k", lw=0.8))
    ax.text(xr3 + 0.3, -0.55, "periphery weld\n(or gasket)", fontsize=9, va="top")
    ax.plot([-0.4, xr3 + 0.3], [0, 0], color=GREY, lw=0.6, ls=":")
    ax.text(-0.45, 0.0, "channel\ncentre", fontsize=8, color=GREY, ha="right", va="center")
    ytop = zmid + H / 2 + 0.75
    dim_h(ax, xr1, xr2, ytop, r"$w_{ro}$", off=0.04, fs=11)
    dim_h(ax, xr2, xr3, ytop, r"$b_m$", off=0.04, fs=11)
    ax.text((xr1 + xr2) / 2, ytop + 0.42, "run-out", ha="center", fontsize=9, color=RED)
    ax.text((xr2 + xr3) / 2, ytop + 0.42, "flat rim", ha="center", fontsize=9, color=RED)
    for xx in (xr1, xr2, xr3):
        ax.plot([xx, xx], [zrim + t, ytop], color=RED, lw=0.6, ls=":")
    dim_v(ax, xr3 + 1.15, zrim, zmid, r"$h_{rim}$", off=0.08, fs=11)
    ax.plot([xr1, xr3 + 1.15], [zmid, zmid], color=RED, lw=0.6, ls=":"); ax.plot([xr3, xr3 + 1.15], [zrim, zrim], color=RED, lw=0.6, ls=":")
    ax.text(xr3 + 1.15, zmid + 0.25, "rim offset from\ncorrugation mid-plane", fontsize=8.5, color=RED, va="bottom")
    # ---- (e2) apex detail (plan-view zoom): continuous vs half-pitch offset ridges
    bx, w, h, y0 = 10.2, 3.6, 2.6, -1.1
    beta = np.radians(60); cot = 1 / np.tan(beta)
    for j, (title, shift) in enumerate((("continuous: crest meets crest", 0.0), ("offset: crest meets valley", 0.25))):
        x0 = bx + j * 4.6; xa = x0 + w / 2
        ax.add_patch(Polygon([[x0, y0], [x0 + w, y0], [x0 + w, y0 + h], [x0, y0 + h]], closed=True, fc="#f4f6f7", ec="k", lw=1.2))
        xs = np.linspace(x0, x0 + w, 300)
        for k in np.arange(-3.0, 3.5, 0.5):
            for phase, ls in ((0.0, "-"), (0.25, "--")):
                kk = y0 + 1.3 + k + phase
                ys_l = kk + (xa - xs) * cot
                ys_r = kk + shift + (xs - xa) * cot
                yl = np.where((xs <= xa) & (ys_l >= y0) & (ys_l <= y0 + h), ys_l, np.nan)
                yr = np.where((xs >= xa) & (ys_r >= y0) & (ys_r <= y0 + h), ys_r, np.nan)
                ax.plot(xs, yl, color=BLUE, lw=1.0, ls=ls); ax.plot(xs, yr, color=BLUE, lw=1.0, ls=ls)
        ax.plot([xa, xa], [y0, y0 + h], color=GREY, ls="-.", lw=0.9)
        ax.text(x0 + w / 2, y0 + h + 0.12, title, ha="center", fontsize=10, color=BLUE)
    ax.plot([], [], color=BLUE, ls="-", label="crest line"); ax.plot([], [], color=BLUE, ls="--", label="valley line")
    ax.legend(loc="upper right", bbox_to_anchor=(0.995, 1.0), fontsize=9, frameon=False, ncol=2)
    xa0 = bx + w / 2
    ax.annotate("", (xa0 + 0.3, y0 - 0.22), (xa0 - 0.3, y0 - 0.22), arrowprops=dict(arrowstyle="<->", color=RED, lw=1.1, shrinkA=0, shrinkB=0))
    ax.plot([xa0 - 0.3, xa0 - 0.3], [y0, y0 - 0.22], color=RED, lw=0.6, ls=":"); ax.plot([xa0 + 0.3, xa0 + 0.3], [y0, y0 - 0.22], color=RED, lw=0.6, ls=":")
    ax.text(xa0, y0 - 0.32, r"$w_{apex}$ (flat band, if any) or $r_{apex}$ (kink fillet radius)", ha="center", va="top", fontsize=9, color=RED)
    # captions
    ax.text(0.0, -2.35, "(e1) Rim: corrugation depth tapers to zero over the run-out width $w_{ro}$, then a flat rim of width $b_m$ offset by $h_{rim}$ meets the paired plate and is welded (or gasketed).  "
                        "Give the same items for the port rim ($b_{port}$, $h_{port}$).", fontsize=9, color=GREY)
    ax.text(0.0, -2.62, "(e2) Apex: state whether crest lines run continuously across the apex line or are offset by half a pitch, and the kink geometry ($r_{apex}$ or $w_{apex}$).  "
                        "For a W pattern give the number and spacing of apex lines.", fontsize=9, color=GREY)
    ax.set_xlim(-1.4, 19.2); ax.set_ylim(-2.85, 2.75); ax.set_aspect("equal"); ax.axis("off")
    ax.set_title("(e) Edge (rim) and apex details", loc="left", fontsize=12, fontweight="bold")


if __name__ == "__main__":
    fig = plt.figure(figsize=(17, 25))
    gs = fig.add_gridspec(4, 2, width_ratios=[1, 1.3], height_ratios=[1.25, 1.0, 0.85, 0.75])
    panel_plan(fig.add_subplot(gs[0, 0])); panel_section(fig.add_subplot(gs[0, 1]))
    panel_axes(fig.add_subplot(gs[1, :])); panel_stack(fig.add_subplot(gs[2, :])); panel_edge(fig.add_subplot(gs[3, :]))
    fig.suptitle("SPHX chevron plate — dimension definition sheet  (fill values in docs/04_plate_dimension_form.md)", fontsize=14, fontweight="bold", y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.975))
    fig.savefig("/home/user/SPHX_EQV/docs/fig/plate_dimensions.png", dpi=150)
    fig.savefig("/home/user/SPHX_EQV/docs/fig/plate_dimensions.svg")
    print("saved")
