# ─────────────────────────────────────────────────────────────────────────────
# STDLIB / THIRD-PARTY IMPORTS
# ─────────────────────────────────────────────────────────────────────────────
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
import numpy as np
from scipy.linalg import inv, eigh
from scipy.optimize import minimize, differential_evolution
import json
import os
import math
import warnings
warnings.filterwarnings("ignore")

import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch
import matplotlib.ticker as mticker
from matplotlib.colors import Normalize
import matplotlib.cm as cm

import openpyxl
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side

# ─────────────────────────────────────────────────────────────────────────────
# THEME  (dark engineering palette — v2 with per-module colors)
# ─────────────────────────────────────────────────────────────────────────────
BG_PRIMARY   = "#0a0e17"
BG_SECONDARY = "#111827"
BG_PANEL     = "#161d2e"
BG_CARD      = "#1c2333"
BORDER       = "#252f42"
BORDER_BRIGHT= "#3d4f6b"
ACCENT       = "#00C8FF"      # default cyan-blue
ACCENT2      = "#f59e0b"      # amber
ACCENT_DIM   = "#003d52"
SUCCESS      = "#10d98a"
WARNING      = "#fbbf24"
DANGER       = "#f87171"
TEXT_PRIMARY = "#e8edf5"
TEXT_SEC     = "#8fa3bf"
TEXT_MUTED   = "#3d4f6b"

PLT_BG       = "#0a0e17"
PLT_AX_BG    = "#111827"
PLT_GRID     = "#1e2d42"
PLT_TEXT     = "#8fa3bf"

# ── Per-module color palette ──────────────────────────────────────────────────
# Each module gets: (accent_color, dim_bg, header_gradient_stop)
MODULE_COLORS = {
    "database":     ("#00C8FF", "#002a3d", "#001f2e"),   # Cyan
    "constitutive": ("#a78bfa", "#2a1f4a", "#1e1638"),   # Violet
    "micro":        ("#34d399", "#0d3326", "#071f18"),   # Emerald
    "laminate":     ("#f59e0b", "#3a2700", "#291b00"),   # Amber
    "failure":      ("#f87171", "#3d1414", "#2a0d0d"),   # Red
    "stress_strain":("#fb923c", "#3a1900", "#291200"),   # Orange
    "optimization": ("#4ade80", "#0d3318", "#07200f"),   # Green
    "sensitivity":  ("#e879f9", "#3a0d4a", "#270635"),   # Fuchsia
    "visualization":("#38bdf8", "#002a42", "#001d2e"),   # Sky Blue
    "export":       ("#94a3b8", "#1e2535", "#151c28"),   # Slate
    "education":    ("#fcd34d", "#3a2e00", "#271f00"),   # Yellow
    "calculator":   ("#67e8f9", "#003540", "#002530"),   # Teal
}

def module_accent(mid):
    return MODULE_COLORS.get(mid, (ACCENT, ACCENT_DIM, BG_SECONDARY))[0]

def module_dim(mid):
    return MODULE_COLORS.get(mid, (ACCENT, ACCENT_DIM, BG_SECONDARY))[1]

def module_hdr(mid):
    return MODULE_COLORS.get(mid, (ACCENT, ACCENT_DIM, BG_SECONDARY))[2]

# Active module tracker (set when a page is built)
_ACTIVE_MODULE = "database"

matplotlib.rcParams.update({
    "font.family": "monospace", "font.size": 9,
    "figure.facecolor": PLT_BG, "axes.facecolor": PLT_AX_BG,
    "axes.edgecolor": PLT_GRID, "xtick.color": PLT_TEXT,
    "ytick.color": PLT_TEXT, "text.color": PLT_TEXT,
    "grid.color": PLT_GRID, "legend.facecolor": BG_CARD,
    "legend.edgecolor": BORDER,
})

# ─────────────────────────────────────────────────────────────────────────────
# MATERIAL DATABASE
# ─────────────────────────────────────────────────────────────────────────────
FIBER_DB = {
    "Carbon Fiber (IM7)": {
        "E1":230.0,"E2":15.0,"nu12":0.20,"nu23":0.25,"G12":27.0,"G23":7.0,
        "density":1.78,"Xt":3500.0,"Xc":1600.0,"Yt":50.0,"Yc":170.0,"S12":90.0,
    },
    "Carbon Fiber (AS4)": {
        "E1":228.0,"E2":14.0,"nu12":0.20,"nu23":0.25,"G12":25.0,"G23":6.5,
        "density":1.79,"Xt":3500.0,"Xc":1480.0,"Yt":48.0,"Yc":200.0,"S12":79.0,
    },
    "E-Glass": {
        "E1":72.0,"E2":72.0,"nu12":0.22,"nu23":0.22,"G12":30.0,"G23":30.0,
        "density":2.54,"Xt":3400.0,"Xc":3400.0,"Yt":3400.0,"Yc":3400.0,"S12":1400.0,
    },
    "S-Glass": {
        "E1":86.0,"E2":86.0,"nu12":0.23,"nu23":0.23,"G12":36.0,"G23":36.0,
        "density":2.49,"Xt":4600.0,"Xc":4600.0,"Yt":4600.0,"Yc":4600.0,"S12":1800.0,
    },
    "Kevlar-49": {
        "E1":125.0,"E2":8.0,"nu12":0.36,"nu23":0.40,"G12":2.9,"G23":2.5,
        "density":1.44,"Xt":3600.0,"Xc":480.0,"Yt":30.0,"Yc":138.0,"S12":49.0,
    },
    "Basalt Fiber": {
        "E1":89.0,"E2":89.0,"nu12":0.26,"nu23":0.26,"G12":35.0,"G23":35.0,
        "density":2.65,"Xt":4840.0,"Xc":4840.0,"Yt":4840.0,"Yc":4840.0,"S12":1900.0,
    },
    "Boron Fiber": {
        "E1":400.0,"E2":400.0,"nu12":0.20,"nu23":0.20,"G12":166.0,"G23":166.0,
        "density":2.63,"Xt":3600.0,"Xc":6900.0,"Yt":3600.0,"Yc":6900.0,"S12":2500.0,
    },
}

MATRIX_DB = {
    "Epoxy (3501-6)": {
        "E":4.2,"nu":0.34,"G":1.57,"density":1.265,
        "Xt":69.0,"Xc":250.0,"S":41.0,
    },
    "Epoxy (8552)": {
        "E":4.67,"nu":0.37,"G":1.70,"density":1.30,
        "Xt":121.0,"Xc":291.0,"S":90.0,
    },
    "Polyester": {
        "E":3.5,"nu":0.38,"G":1.27,"density":1.20,
        "Xt":65.0,"Xc":170.0,"S":35.0,
    },
    "Vinyl Ester": {
        "E":3.8,"nu":0.36,"G":1.40,"density":1.15,
        "Xt":80.0,"Xc":200.0,"S":50.0,
    },
    "PEEK": {
        "E":3.6,"nu":0.40,"G":1.29,"density":1.32,
        "Xt":100.0,"Xc":120.0,"S":55.0,
    },
    "PEKK": {
        "E":3.9,"nu":0.39,"G":1.40,"density":1.30,
        "Xt":105.0,"Xc":130.0,"S":60.0,
    },
    "Polyimide (PMR-15)": {
        "E":4.5,"nu":0.36,"G":1.65,"density":1.32,
        "Xt":38.0,"Xc":180.0,"S":48.0,
    },
}

# ─────────────────────────────────────────────────────────────────────────────
# CORE MECHANICS: STIFFNESS / COMPLIANCE
# ─────────────────────────────────────────────────────────────────────────────

def compliance_from_engineering(E1,E2,E3,nu12,nu13,nu23,G12,G13,G23):
    nu21 = nu12*E2/E1; nu31 = nu13*E3/E1; nu32 = nu23*E3/E2
    S = np.array([
        [ 1/E1,    -nu21/E2, -nu31/E3,  0,      0,      0     ],
        [-nu12/E1,  1/E2,    -nu32/E3,  0,      0,      0     ],
        [-nu13/E1, -nu23/E2,  1/E3,     0,      0,      0     ],
        [ 0,        0,        0,        1/G23,  0,      0     ],
        [ 0,        0,        0,        0,      1/G13,  0     ],
        [ 0,        0,        0,        0,      0,      1/G12 ],
    ])
    return S

def stiffness_isotropic(E, nu):
    G = E/(2*(1+nu)); lam = nu*E/((1+nu)*(1-2*nu)); mu = G
    C = np.zeros((6,6))
    C[0,0]=C[1,1]=C[2,2] = lam+2*mu
    C[0,1]=C[1,0]=C[0,2]=C[2,0]=C[1,2]=C[2,1] = lam
    C[3,3]=C[4,4]=C[5,5] = mu
    return C

def stiffness_from_compliance(E1,E2,E3,nu12,nu13,nu23,G12,G13,G23):
    return inv(compliance_from_engineering(E1,E2,E3,nu12,nu13,nu23,G12,G13,G23))

def stability_check(C):
    """Return eigenvalues and PASS/FAIL for positive definiteness."""
    evals,_ = eigh(C)
    passed = bool(np.all(evals > 0))
    return evals, passed

def voigt_notation_labels():
    return ["11","22","33","23","13","12"]

# ─────────────────────────────────────────────────────────────────────────────
# MICROMECHANICS MODELS
# ─────────────────────────────────────────────────────────────────────────────

def rule_of_mixtures(Ef1, Ef2, nuf12, Gf12, Em, num, Vf):
    """Voigt (longitudinal) + Reuss (transverse) + Halpin-Tsai."""
    Vm = 1.0 - Vf
    E1  = Ef1*Vf + Em*Vm                      # Voigt (longitudinal)
    nu12 = nuf12*Vf + num*Vm
    E2_inv = Vf/Ef2 + Vm/Em                   # Reuss (transverse)
    E2  = 1.0/E2_inv if E2_inv else Ef2
    G12_inv = Vf/Gf12 + Vm/(Em/(2*(1+num)))
    G12 = 1.0/G12_inv if G12_inv else Gf12
    nu21 = nu12*E2/E1
    return {"E1":E1,"E2":E2,"nu12":nu12,"nu21":nu21,"G12":G12}

def halpin_tsai(Ef1, Ef2, nuf12, Gf12, Em, num, Vf, xi_E=2.0, xi_G=1.0):
    """Halpin-Tsai semi-empirical model."""
    Vm = 1.0-Vf
    Gm = Em/(2*(1+num))

    def ht(Pf, Pm, xi):
        eta = (Pf/Pm - 1)/(Pf/Pm + xi)
        return Pm*(1 + xi*eta*Vf)/(1 - eta*Vf)

    E1  = Ef1*Vf + Em*Vm
    nu12 = nuf12*Vf + num*Vm
    E2  = ht(Ef2, Em, xi_E)
    G12 = ht(Gf12, Gm, xi_G)
    nu21 = nu12*E2/E1
    return {"E1":E1,"E2":E2,"nu12":nu12,"nu21":nu21,"G12":G12}

def mori_tanaka(Ef1, Ef2, nuf12, Gf12, Em, num, Vf):
    """Mori-Tanaka (simplified isotropic fiber in isotropic matrix)."""
    Vm = 1.0-Vf
    Km = Em/(3*(1-2*num)); Gm = Em/(2*(1+num))
    # Treating fiber as isotropic with E_f ~ average of Ef1, Ef2
    Ef_avg = (2*Ef1 + Ef2*3)/5
    nuf_avg = nuf12  # approx
    Kf = Ef_avg/(3*(1-2*nuf_avg)); Gff = Ef_avg/(2*(1+nuf_avg))
    # MT effective
    alpha_K = 3*Km+4*Gm
    alpha_G = Gm*(9*Km+8*Gm)/(6*(Km+2*Gm))
    K_eff = Km + Vf*(Kf-Km)/(1 + Vm*(Kf-Km)/(Km+alpha_K/3))
    G_eff = Gm + Vf*(Gff-Gm)/(1 + Vm*(Gff-Gm)/(Gm+alpha_G))
    E_eff = 9*K_eff*G_eff/(3*K_eff+G_eff)
    nu_eff = (3*K_eff-2*G_eff)/(2*(3*K_eff+G_eff))
    # Anisotropy correction for longitudinal
    E1_rom = Ef1*Vf + Em*Vm
    nu12_rom = nuf12*Vf + num*Vm
    return {"E1":E1_rom,"E2":E_eff,"nu12":nu12_rom,"nu21":nu12_rom*E_eff/E1_rom,"G12":G_eff}

def self_consistent(Ef1, Ef2, nuf12, Gf12, Em, num, Vf, iters=50):
    """Self-Consistent (iterative, simplified)."""
    # Start from Halpin-Tsai
    res = halpin_tsai(Ef1, Ef2, nuf12, Gf12, Em, num, Vf)
    Vm = 1.0-Vf
    Gm = Em/(2*(1+num))
    Gf = Gf12
    G = res["G12"]
    for _ in range(iters):
        G_new = (Vf*Gf + Vm*Gm + (Vf*Gm + Vm*Gf)*G/max(G,1e-9)) / \
                (Vf + Vm + (Vf*Gm + Vm*Gf)/max(G,1e-9))
        if abs(G_new - G) < 1e-6: break
        G = G_new
    res["G12"] = G
    return res

# ─────────────────────────────────────────────────────────────────────────────
# CLASSICAL LAMINATION THEORY
# ─────────────────────────────────────────────────────────────────────────────

def Q_matrix(E1, E2, nu12, G12):
    """Reduced stiffness matrix [Q] for a lamina in principal axes."""
    nu21 = nu12*E2/E1
    denom = 1.0 - nu12*nu21
    Q = np.zeros((3,3))
    Q[0,0] = E1/denom
    Q[1,1] = E2/denom
    Q[0,1] = Q[1,0] = nu12*E2/denom
    Q[2,2] = G12
    return Q

def Qbar_matrix(Q, theta_deg):
    """Transformed reduced stiffness [Q̄] at angle theta."""
    t = math.radians(theta_deg)
    c = math.cos(t); s = math.sin(t)
    c2=c*c; s2=s*s; cs=c*s
    Q11,Q22,Q12,Q66 = Q[0,0],Q[1,1],Q[0,1],Q[2,2]
    Q16_ = Q[0,2] if Q.shape[0]>2 else 0
    Q26_ = Q[1,2] if Q.shape[0]>2 else 0

    Qb = np.zeros((3,3))
    Qb[0,0] = Q11*c2**2 + 2*(Q12+2*Q66)*s2*c2 + Q22*s2**2
    Qb[1,1] = Q11*s2**2 + 2*(Q12+2*Q66)*s2*c2 + Q22*c2**2
    Qb[0,1] = Qb[1,0] = (Q11+Q22-4*Q66)*s2*c2 + Q12*(c2**2+s2**2)
    Qb[2,2] = (Q11+Q22-2*Q12-2*Q66)*s2*c2 + Q66*(c2**2+s2**2)
    Qb[0,2] = Qb[2,0] = (Q11-Q12-2*Q66)*c**3*s - (Q22-Q12-2*Q66)*c*s**3
    Qb[1,2] = Qb[2,1] = (Q11-Q12-2*Q66)*c*s**3 - (Q22-Q12-2*Q66)*c**3*s
    return Qb

def ABD_matrix(plies):
    """
    plies: list of dicts with keys E1,E2,nu12,G12,theta,t (thickness mm → m)
    Returns A (3x3), B (3x3), D (3x3) in N/m, N, N·m.
    """
    A = np.zeros((3,3)); B = np.zeros((3,3)); D = np.zeros((3,3))
    total_t = sum(p["t"]*1e-3 for p in plies)
    z0 = -total_t/2.0
    z = z0
    for p in plies:
        t_m = p["t"]*1e-3
        z1 = z + t_m
        zm = (z+z1)/2.0
        Qb = Qbar_matrix(Q_matrix(p["E1"],p["E2"],p["nu12"],p["G12"]), p["theta"])
        A += Qb*(z1-z)
        B += Qb*(z1**2-z**2)/2.0
        D += Qb*(z1**3-z**3)/3.0
        z = z1
    return A,B,D

def laminate_engineering_constants(A,B,D,t_total_mm):
    """Engineering constants from ABD (t in mm)."""
    t = t_total_mm*1e-3
    try:
        a = inv(A)
        Ex  = 1.0/(a[0,0]*t)
        Ey  = 1.0/(a[1,1]*t)
        Gxy = 1.0/(a[2,2]*t)
        nuxy = -a[0,1]/a[0,0]
        return {"Ex":Ex/1e9,"Ey":Ey/1e9,"Gxy":Gxy/1e9,"nuxy":nuxy}
    except Exception:
        return {}

# ─────────────────────────────────────────────────────────────────────────────
# FAILURE CRITERIA
# ─────────────────────────────────────────────────────────────────────────────

def failure_tsai_hill(s1,s2,s12, Xt,Xc,Yt,Yc,S12):
    X = Xt if s1>=0 else Xc
    Y = Yt if s2>=0 else Yc
    FI = (s1/X)**2 - s1*s2/X**2 + (s2/Y)**2 + (s12/S12)**2
    return FI

def failure_tsai_wu(s1,s2,s12, Xt,Xc,Yt,Yc,S12):
    F1 = 1/Xt - 1/Xc; F2 = 1/Yt - 1/Yc
    F11 = 1/(Xt*Xc); F22 = 1/(Yt*Yc); F66 = 1/S12**2
    F12 = -0.5*math.sqrt(F11*F22)
    FI = F1*s1 + F2*s2 + F11*s1**2 + F22*s2**2 + F66*s12**2 + 2*F12*s1*s2
    return FI

def failure_max_stress(s1,s2,s12, Xt,Xc,Yt,Yc,S12):
    fi1 = s1/Xt if s1>=0 else -s1/Xc
    fi2 = s2/Yt if s2>=0 else -s2/Yc
    fi12 = abs(s12)/S12
    return max(fi1, fi2, fi12)

def failure_hashin(s1,s2,s12, Xt,Xc,Yt,Yc,S12):
    if s1 >= 0:
        FI_ft = (s1/Xt)**2 + (s12/S12)**2   # fiber tension
    else:
        FI_ft = (s1/Xc)**2                   # fiber compression
    if s2 >= 0:
        FI_mt = (s2/Yt)**2 + (s12/S12)**2   # matrix tension
    else:
        FI_mt = (s2/(2*S12))**2 + ((Yc/(2*S12))**2-1)*s2/Yc + (s12/S12)**2
    return max(FI_ft, FI_mt)

# ─────────────────────────────────────────────────────────────────────────────
# COORDINATE TRANSFORMATION
# ─────────────────────────────────────────────────────────────────────────────

def stress_transform(sigma_glob, theta_deg):
    """Transform global stress to local (fiber) axes."""
    t = math.radians(theta_deg)
    c=math.cos(t); s=math.sin(t)
    T = np.array([
        [ c*c,   s*s,   2*c*s ],
        [ s*s,   c*c,  -2*c*s ],
        [-c*s,   c*s,   c*c-s*s],
    ])
    return T @ sigma_glob

# ─────────────────────────────────────────────────────────────────────────────
# HELPER: styled matplotlib axes
# ─────────────────────────────────────────────────────────────────────────────

def style_ax(ax):
    ax.set_facecolor(PLT_AX_BG)
    ax.tick_params(colors=PLT_TEXT, labelsize=8)
    ax.xaxis.label.set_color(PLT_TEXT); ax.yaxis.label.set_color(PLT_TEXT)
    for sp in ax.spines.values(): sp.set_color(PLT_GRID)
    ax.grid(True, color=PLT_GRID, lw=0.6, ls="--", alpha=0.5)

def fmt_val(v):
    av = abs(v)
    if av==0: return "0"
    if av >= 1e4 or (av < 1e-3 and av > 0): return f"{v:.3e}"
    if av >= 100: return f"{v:.2f}"
    if av >= 10:  return f"{v:.3f}"
    return f"{v:.4f}"

# ─────────────────────────────────────────────────────────────────────────────
# WIDGET HELPERS
# ─────────────────────────────────────────────────────────────────────────────

def make_label(parent, text, font=("Helvetica",10), fg=TEXT_SEC, bg=None, **kw):
    if bg is None: bg = parent.cget("bg") if hasattr(parent,"cget") else BG_PRIMARY
    return tk.Label(parent, text=text, font=font, fg=fg, bg=bg, **kw)

def make_entry(parent, default="", width=12):
    e = tk.Entry(parent, font=("Courier",10), width=width,
                 bg=BG_PANEL, fg=TEXT_PRIMARY, insertbackground=ACCENT,
                 relief="flat", highlightthickness=1,
                 highlightbackground=BORDER, highlightcolor=ACCENT)
    e.insert(0, str(default))
    return e

def section_header(parent, text, color=None):
    f = tk.Frame(parent, bg=BG_SECONDARY)
    f.pack(fill="x", padx=10, pady=(10,3))
    clr = color or TEXT_MUTED
    tk.Label(f, text=text, font=("Courier",8,"bold"), bg=BG_SECONDARY, fg=clr).pack(side="left")
    tk.Frame(f, bg=color or BORDER, height=1).pack(side="left", fill="x", expand=True, padx=(6,0))

def accent_btn(parent, text, command, bg=ACCENT, fg="#000000", **kw):
    return tk.Button(parent, text=text, command=command,
                     font=("Helvetica",10,"bold"), bg=bg, fg=fg,
                     activebackground=ACCENT2, activeforeground="#000",
                     relief="flat", cursor="hand2", pady=8, **kw)

def nav_btn(parent, text, command, active=False):
    bg = BG_PANEL if not active else ACCENT_DIM
    fg = TEXT_PRIMARY if active else TEXT_SEC
    btn = tk.Button(parent, text=text, command=command,
                    font=("Helvetica",10), bg=bg, fg=fg,
                    activebackground=ACCENT_DIM, activeforeground=ACCENT,
                    relief="flat", cursor="hand2", anchor="w",
                    padx=14, pady=8)
    btn.pack(fill="x", pady=1)
    return btn

def scrollable_frame(parent):
    canvas = tk.Canvas(parent, bg=parent.cget("bg"), highlightthickness=0)
    vsb = ttk.Scrollbar(parent, orient="vertical", command=canvas.yview)
    frame = tk.Frame(canvas, bg=parent.cget("bg"))
    frame.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
    canvas.create_window((0,0), window=frame, anchor="nw")
    canvas.configure(yscrollcommand=vsb.set)
    vsb.pack(side="right", fill="y")
    canvas.pack(side="left", fill="both", expand=True)
    return frame

def matrix_grid(parent, M, row_labels=None, col_labels=None, title="", accent=None):
    """Render a matrix as a colored grid of labels."""
    for w in parent.winfo_children(): w.destroy()
    r,c = M.shape
    diag_color = accent or ACCENT
    offdiag_color = ACCENT2
    if title:
        tk.Label(parent, text=title, font=("Helvetica",11,"bold"),
                 bg=parent.cget("bg"), fg=TEXT_PRIMARY).grid(
                     row=0, column=0, columnspan=c+1, sticky="w", pady=(0,6))
    off=1 if title else 0
    if col_labels:
        for j,lbl in enumerate(col_labels):
            tk.Label(parent, text=lbl, font=("Courier",9,"bold"),
                     bg=BG_PANEL, fg=TEXT_MUTED, padx=8, pady=4,
                     width=10, anchor="center", relief="flat"
                     ).grid(row=off, column=j+1, padx=2, pady=1)
    max_abs = np.max(np.abs(M)) if np.max(np.abs(M))>0 else 1.0
    for i in range(r):
        if row_labels:
            tk.Label(parent, text=row_labels[i], font=("Courier",9,"bold"),
                     bg=BG_PANEL, fg=TEXT_MUTED, padx=8, pady=4,
                     width=4, anchor="center").grid(row=i+off+1, column=0, padx=2, pady=1)
        for j in range(c):
            v = M[i,j]
            is_diag = (i==j)
            is_nz = abs(v) > max_abs*1e-9
            # Tinted backgrounds
            if is_diag and is_nz:
                # Slightly tint BG_CARD toward the accent color
                bg_c = "#1a2a10"; fg_c = diag_color
            elif is_nz:
                bg_c = "#0e1f2f"; fg_c = offdiag_color
            else:
                bg_c = BG_CARD; fg_c = TEXT_MUTED
            tk.Label(parent, text=fmt_val(v), font=("Courier",10,"bold"),
                     bg=bg_c, fg=fg_c, padx=10, pady=5, width=10,
                     anchor="center", relief="flat"
                     ).grid(row=i+off+1, column=j+1, padx=2, pady=1)

# ─────────────────────────────────────────────────────────────────────────────
# MODULE PAGES (each is a Frame subclass)
# ─────────────────────────────────────────────────────────────────────────────

class BasePage(tk.Frame):
    def __init__(self, parent, app, mid="database"):
        super().__init__(parent, bg=BG_PRIMARY)
        self.app = app
        self.mid = mid
        self.color = module_accent(mid)
        self.dim   = module_dim(mid)
        self.hdr_bg= module_hdr(mid)

    def page_header(self, title, subtitle=""):
        # Outer header with module-tinted background
        h = tk.Frame(self, bg=self.hdr_bg, padx=20, pady=14)
        h.pack(fill="x")

        # Colored left accent bar
        tk.Frame(h, bg=self.color, width=4).pack(side="left", fill="y", padx=(0,14))

        inner = tk.Frame(h, bg=self.hdr_bg)
        inner.pack(side="left", fill="both", expand=True)

        tk.Label(inner, text=title, font=("Helvetica",15,"bold"),
                 bg=self.hdr_bg, fg=TEXT_PRIMARY).pack(anchor="w")
        if subtitle:
            tk.Label(inner, text=subtitle, font=("Courier",9),
                     bg=self.hdr_bg, fg=self.color).pack(anchor="w", pady=(2,0))

        # Thin colored bottom border
        tk.Frame(self, bg=self.color, height=2).pack(fill="x")

    def embed_figure(self, parent, fig, toolbar=True):
        canvas = FigureCanvasTkAgg(fig, master=parent)
        canvas.get_tk_widget().pack(fill="both", expand=True)
        if toolbar:
            tb_frame = tk.Frame(parent, bg=BG_PRIMARY)
            tb_frame.pack(fill="x")
            tb = NavigationToolbar2Tk(canvas, tb_frame)
            tb.config(bg=BG_PRIMARY)
            tb.update()
        canvas.draw()
        return canvas


# ══════════════════════════════════════════════════════════════════════════════
# 1. MATERIAL DATABASE MODULE
# ══════════════════════════════════════════════════════════════════════════════

class MaterialDBPage(BasePage):
    def __init__(self, parent, app):
        super().__init__(parent, app, mid="database")
        self.page_header("Material Database", "Fibers · Matrices · Custom Materials")
        self._build()

    def _build(self):
        pw = tk.PanedWindow(self, orient="horizontal", bg=BG_PRIMARY, sashwidth=4)
        pw.pack(fill="both", expand=True, padx=10, pady=10)

        # LEFT – fiber list
        left = tk.Frame(pw, bg=BG_SECONDARY, width=360)
        pw.add(left, minsize=300)

        section_header(left, "FIBER DATABASE")
        self.fiber_list = tk.Listbox(left, bg=BG_PANEL, fg=TEXT_PRIMARY,
                                     selectbackground=self.dim, selectforeground=self.color,
                                     font=("Courier",10), relief="flat",
                                     activestyle="none", height=10, bd=0,
                                     highlightthickness=0)
        for f in FIBER_DB: self.fiber_list.insert("end", f)
        self.fiber_list.pack(fill="x", padx=10, pady=(0,6))
        self.fiber_list.bind("<<ListboxSelect>>", self._show_fiber)

        section_header(left, "MATRIX DATABASE")
        self.matrix_list = tk.Listbox(left, bg=BG_PANEL, fg=TEXT_PRIMARY,
                                      selectbackground=self.dim, selectforeground=self.color,
                                      font=("Courier",10), relief="flat",
                                      activestyle="none", height=8, bd=0,
                                      highlightthickness=0)
        for m in MATRIX_DB: self.matrix_list.insert("end", m)
        self.matrix_list.pack(fill="x", padx=10, pady=(0,6))
        self.matrix_list.bind("<<ListboxSelect>>", self._show_matrix)

        # RIGHT – property display
        right = tk.Frame(pw, bg=BG_PRIMARY)
        pw.add(right, minsize=400)

        self.prop_title = tk.Label(right, text="Select a material from the list",
                                   font=("Helvetica",13,"bold"),
                                   bg=BG_PRIMARY, fg=self.color)
        self.prop_title.pack(anchor="w", padx=20, pady=(14,6))
        self.prop_frame = tk.Frame(right, bg=BG_PRIMARY)
        self.prop_frame.pack(fill="both", expand=True, padx=20)

    def _show_fiber(self, event):
        sel = self.fiber_list.curselection()
        if not sel: return
        name = self.fiber_list.get(sel[0])
        self._render_props(name, FIBER_DB[name], "fiber")

    def _show_matrix(self, event):
        sel = self.matrix_list.curselection()
        if not sel: return
        name = self.matrix_list.get(sel[0])
        self._render_props(name, MATRIX_DB[name], "matrix")

    def _render_props(self, name, props, kind):
        self.prop_title.config(text=name)
        for w in self.prop_frame.winfo_children(): w.destroy()
        labels = {
            "E1":"E₁  Longitudinal Modulus","E2":"E₂  Transverse Modulus",
            "nu12":"ν₁₂  Poisson Ratio","G12":"G₁₂  Shear Modulus",
            "density":"ρ  Density","Xt":"Xₜ  Tensile Strength (long.)",
            "Xc":"Xc  Compressive Strength (long.)","Yt":"Yₜ  Tensile Strength (trans.)",
            "Yc":"Yc  Compressive Strength (trans.)","S12":"S₁₂  Shear Strength",
            "E":"E  Young's Modulus","nu":"ν  Poisson Ratio","G":"G  Shear Modulus",
            "Xt":"Xₜ  Tensile Strength","Xc":"Xc  Compressive Strength","S":"S  Shear Strength",
        }
        units = {
            "E1":"GPa","E2":"GPa","G12":"GPa","G23":"GPa","density":"g/cm³",
            "Xt":"MPa","Xc":"MPa","Yt":"MPa","Yc":"MPa","S12":"MPa",
            "E":"GPa","G":"GPa","S":"MPa",
        }
        grid = tk.Frame(self.prop_frame, bg=BG_PRIMARY)
        grid.pack(fill="both")
        for idx,(k,v) in enumerate(props.items()):
            row = idx//2; col = idx%2
            card = tk.Frame(grid, bg=BG_CARD, padx=12, pady=8)
            card.grid(row=row, column=col, padx=5, pady=4, sticky="nsew")
            grid.columnconfigure(col, weight=1)
            lbl = labels.get(k, k)
            unit = units.get(k,"")
            tk.Label(card, text=lbl, font=("Helvetica",9), bg=BG_CARD, fg=TEXT_MUTED,
                     anchor="w").pack(anchor="w")
            tk.Label(card, text=f"{v}  {unit}", font=("Courier",13,"bold"),
                     bg=BG_CARD, fg=self.color, anchor="w").pack(anchor="w", pady=(2,0))


# ══════════════════════════════════════════════════════════════════════════════
# 2. STIFFNESS / CONSTITUTIVE MATRIX MODULE  (enhanced from existing app)
# ══════════════════════════════════════════════════════════════════════════════

MATERIAL_PARAMS = {
    "Isotropic": {
        "Young's Modulus E [GPa]":  {"default":200.0,"key":"E"},
        "Poisson's Ratio ν":         {"default":0.30, "key":"nu"},
    },
    "Transversely Isotropic": {
        "E₁ (fiber) [GPa]":         {"default":150.0,"key":"E1"},
        "E₂ (transverse) [GPa]":    {"default":10.0, "key":"E2"},
        "ν₁₂":                       {"default":0.25, "key":"nu12"},
        "ν₂₃":                       {"default":0.40, "key":"nu23"},
        "G₁₂ [GPa]":                 {"default":7.0,  "key":"G12"},
        "G₂₃ [GPa]":                 {"default":4.0,  "key":"G23"},
    },
    "Orthotropic": {
        "E₁ [GPa]":  {"default":210.0,"key":"E1"}, "E₂ [GPa]":{"default":120.0,"key":"E2"},
        "E₃ [GPa]":  {"default":80.0, "key":"E3"}, "ν₁₂":      {"default":0.28, "key":"nu12"},
        "ν₁₃":       {"default":0.32, "key":"nu13"},"ν₂₃":      {"default":0.35, "key":"nu23"},
        "G₁₂ [GPa]": {"default":80.0, "key":"G12"}, "G₁₃ [GPa]":{"default":65.0,"key":"G13"},
        "G₂₃ [GPa]": {"default":50.0, "key":"G23"},
    },
    "Monoclinic": {
        "E₁ [GPa]":{"default":180.0,"key":"E1"},"E₂ [GPa]":{"default":130.0,"key":"E2"},
        "E₃ [GPa]":{"default":90.0,"key":"E3"},"ν₁₂":{"default":0.30,"key":"nu12"},
        "ν₁₃":{"default":0.28,"key":"nu13"},"ν₂₃":{"default":0.25,"key":"nu23"},
        "G₁₂ [GPa]":{"default":70.0,"key":"G12"},"G₁₃ [GPa]":{"default":55.0,"key":"G13"},
        "G₂₃ [GPa]":{"default":45.0,"key":"G23"},"C₁₆ [GPa]":{"default":5.0,"key":"C16"},
        "C₂₆ [GPa]":{"default":3.0,"key":"C26"},"C₃₆ [GPa]":{"default":2.0,"key":"C36"},
        "C₄₅ [GPa]":{"default":4.0,"key":"C45"},
    },
    "Cubic": {
        "C₁₁ [GPa]":{"default":165.0,"key":"C11"},
        "C₁₂ [GPa]":{"default":64.0, "key":"C12"},
        "C₄₄ [GPa]":{"default":79.5, "key":"C44"},
    },
    "Hexagonal": {
        "C₁₁ [GPa]":{"default":160.0,"key":"C11"},"C₁₂ [GPa]":{"default":90.0,"key":"C12"},
        "C₁₃ [GPa]":{"default":66.0,"key":"C13"},"C₃₃ [GPa]":{"default":181.0,"key":"C33"},
        "C₄₄ [GPa]":{"default":46.0,"key":"C44"},
    },
    "Anisotropic (21 constants)": {
        "C₁₁":{"default":200.0,"key":"C11"},"C₁₂":{"default":80.0,"key":"C12"},
        "C₁₃":{"default":60.0,"key":"C13"},"C₁₄":{"default":5.0,"key":"C14"},
        "C₁₅":{"default":3.0,"key":"C15"},"C₁₆":{"default":2.0,"key":"C16"},
        "C₂₂":{"default":180.0,"key":"C22"},"C₂₃":{"default":55.0,"key":"C23"},
        "C₂₄":{"default":4.0,"key":"C24"},"C₂₅":{"default":2.0,"key":"C25"},
        "C₂₆":{"default":1.5,"key":"C26"},"C₃₃":{"default":150.0,"key":"C33"},
        "C₃₄":{"default":3.0,"key":"C34"},"C₃₅":{"default":2.0,"key":"C35"},
        "C₃₆":{"default":1.0,"key":"C36"},"C₄₄":{"default":70.0,"key":"C44"},
        "C₄₅":{"default":4.0,"key":"C45"},"C₄₆":{"default":2.0,"key":"C46"},
        "C₅₅":{"default":65.0,"key":"C55"},"C₅₆":{"default":3.0,"key":"C56"},
        "C₆₆":{"default":60.0,"key":"C66"},
    },
}

def build_stiffness(mat_type, p):
    if mat_type == "Isotropic":
        return stiffness_isotropic(p["E"], p["nu"])
    elif mat_type == "Transversely Isotropic":
        S = compliance_from_engineering(p["E1"],p["E2"],p["E2"],p["nu12"],p["nu12"],p["nu23"],p["G12"],p["G12"],p["G23"])
        return inv(S)
    elif mat_type == "Orthotropic":
        S = compliance_from_engineering(p["E1"],p["E2"],p["E3"],p["nu12"],p["nu13"],p["nu23"],p["G12"],p["G13"],p["G23"])
        return inv(S)
    elif mat_type == "Monoclinic":
        S = compliance_from_engineering(p["E1"],p["E2"],p["E3"],p["nu12"],p["nu13"],p["nu23"],p["G12"],p["G13"],p["G23"])
        C = inv(S)
        C[0,5]=C[5,0]=p["C16"]; C[1,5]=C[5,1]=p["C26"]
        C[2,5]=C[5,2]=p["C36"]; C[3,4]=C[4,3]=p["C45"]
        return C
    elif mat_type == "Cubic":
        C = np.zeros((6,6))
        C[0,0]=C[1,1]=C[2,2]=p["C11"]
        C[0,1]=C[1,0]=C[0,2]=C[2,0]=C[1,2]=C[2,1]=p["C12"]
        C[3,3]=C[4,4]=C[5,5]=p["C44"]
        return C
    elif mat_type == "Hexagonal":
        C11,C12,C13,C33,C44 = p["C11"],p["C12"],p["C13"],p["C33"],p["C44"]
        C66 = (C11-C12)/2
        C = np.zeros((6,6))
        C[0,0]=C[1,1]=C11; C[2,2]=C33
        C[0,1]=C[1,0]=C12; C[0,2]=C[2,0]=C[1,2]=C[2,1]=C13
        C[3,3]=C[4,4]=C44; C[5,5]=C66
        return C
    elif mat_type == "Anisotropic (21 constants)":
        keys = ["C11","C12","C13","C14","C15","C16",
                "C22","C23","C24","C25","C26",
                "C33","C34","C35","C36",
                "C44","C45","C46",
                "C55","C56","C66"]
        idx = [(0,0),(0,1),(0,2),(0,3),(0,4),(0,5),
               (1,1),(1,2),(1,3),(1,4),(1,5),
               (2,2),(2,3),(2,4),(2,5),
               (3,3),(3,4),(3,5),(4,4),(4,5),(5,5)]
        C = np.zeros((6,6))
        for k,ij in zip(keys,idx):
            i,j = ij
            C[i,j] = C[j,i] = p.get(k,0.0)
        return C
    return np.zeros((6,6))


class ConstitutivePage(BasePage):
    def __init__(self, parent, app):
        super().__init__(parent, app, mid="constitutive")
        self.C = None; self.S_mat = None
        self.param_entries = {}
        self.page_header("Constitutive Matrix Laboratory",
                         "[C] Stiffness · [S] Compliance · Stability · Derived Properties")
        self._build()

    def _build(self):
        pw = tk.PanedWindow(self, orient="horizontal", bg=BG_PRIMARY, sashwidth=4)
        pw.pack(fill="both", expand=True)

        # SIDEBAR
        sidebar = tk.Frame(pw, bg=BG_SECONDARY, width=300)
        pw.add(sidebar, minsize=260)
        sidebar.pack_propagate(False)

        sf = scrollable_frame(sidebar)

        section_header(sf, "MATERIAL SYMMETRY CLASS")
        self.mat_var = tk.StringVar(value="")
        self.mat_btns = {}
        for m in MATERIAL_PARAMS:
            b = tk.Button(sf, text=m, font=("Helvetica",10),
                          bg=BG_PANEL, fg=TEXT_SEC, relief="flat",
                          cursor="hand2", anchor="w", padx=10, pady=6,
                          command=lambda mm=m: self._select_mat(mm))
            b.pack(fill="x", padx=8, pady=1)
            self.mat_btns[m] = b

        section_header(sf, "INPUT PARAMETERS")
        self.inp_frame = tk.Frame(sf, bg=BG_SECONDARY)
        self.inp_frame.pack(fill="x", padx=8, pady=4)

        tk.Frame(sf, bg=BORDER, height=1).pack(fill="x", padx=8, pady=8)
        accent_btn(sf, "▶  COMPUTE", self._compute).pack(fill="x", padx=8, pady=4)

        # CONTENT
        content = tk.Frame(pw, bg=BG_PRIMARY)
        pw.add(content, minsize=600)

        # Tab bar
        tab_bar = tk.Frame(content, bg=BG_SECONDARY, height=40)
        tab_bar.pack(fill="x"); tab_bar.pack_propagate(False)
        tk.Frame(content, bg=BORDER, height=1).pack(fill="x")

        self.tab_content = tk.Frame(content, bg=BG_PRIMARY)
        self.tab_content.pack(fill="both", expand=True)

        self.tabs = {}
        self.tab_btns = {}
        tab_defs = [("stiff","[C]  Stiffness"),("compl","[S]  Compliance"),
                    ("stab","⚡  Stability"),("derived","∿  Derived")]
        for key,label in tab_defs:
            btn = tk.Button(tab_bar, text=label, font=("Courier",9),
                            bg=BG_SECONDARY, fg=TEXT_MUTED, relief="flat",
                            cursor="hand2", padx=14, pady=10,
                            command=lambda k=key: self._switch_tab(k))
            btn.pack(side="left")
            self.tab_btns[key] = btn
            fr = tk.Frame(self.tab_content, bg=BG_PRIMARY)
            self.tabs[key] = fr

        self.stiff_grid   = tk.Frame(self.tabs["stiff"], bg=BG_PRIMARY)
        self.stiff_grid.pack(padx=20, pady=20)
        self.compl_grid   = tk.Frame(self.tabs["compl"], bg=BG_PRIMARY)
        self.compl_grid.pack(padx=20, pady=20)
        self.stab_frame   = self.tabs["stab"]
        self.derived_frame= self.tabs["derived"]

        self._switch_tab("stiff")

    def _switch_tab(self, key):
        for k,f in self.tabs.items():
            f.pack_forget()
        self.tabs[key].pack(fill="both", expand=True)
        for k,b in self.tab_btns.items():
            b.config(fg=self.color if k==key else TEXT_MUTED,
                     bg=BG_PRIMARY if k==key else BG_SECONDARY)

    def _select_mat(self, m):
        self.mat_var.set(m)
        for k,b in self.mat_btns.items():
            b.config(bg=self.dim if k==m else BG_PANEL,
                     fg=self.color if k==m else TEXT_SEC)
        for w in self.inp_frame.winfo_children(): w.destroy()
        self.param_entries.clear()
        for lbl,info in MATERIAL_PARAMS[m].items():
            row = tk.Frame(self.inp_frame, bg=BG_SECONDARY)
            row.pack(fill="x", pady=2)
            tk.Label(row, text=lbl, font=("Helvetica",9), bg=BG_SECONDARY,
                     fg=TEXT_SEC, anchor="w", width=22).pack(side="left")
            e = make_entry(row, info["default"], width=9)
            e.pack(side="right", padx=2)
            self.param_entries[info["key"]] = e

    def _compute(self):
        m = self.mat_var.get()
        if not m:
            messagebox.showwarning("CMVL","Select a material type first."); return
        try:
            p = {k: float(e.get()) for k,e in self.param_entries.items()}
            self.C = build_stiffness(m, p)
            self.S_mat = inv(self.C)
        except Exception as ex:
            messagebox.showerror("Compute Error", str(ex)); return

        lbl6 = ["σ₁","σ₂","σ₃","σ₄","σ₅","σ₆"]
        eps6 = ["ε₁","ε₂","ε₃","ε₄","ε₅","ε₆"]
        matrix_grid(self.stiff_grid, self.C, lbl6, eps6, f"Stiffness [C]  (GPa)   —  {m}", accent=self.color)
        matrix_grid(self.compl_grid, self.S_mat, lbl6, eps6, "Compliance [S] = [C]⁻¹  (GPa⁻¹)", accent=self.color)

        # Stability tab
        for w in self.stab_frame.winfo_children(): w.destroy()
        evals, passed = stability_check(self.C)
        color = SUCCESS if passed else DANGER
        status = "✔  STABLE — Positive Definite" if passed else "✗  UNSTABLE — Not Positive Definite"
        tk.Label(self.stab_frame, text=status, font=("Helvetica",13,"bold"),
                 bg=BG_PRIMARY, fg=color).pack(padx=20, pady=14, anchor="w")
        tk.Label(self.stab_frame, text="Eigenvalues of [C]:", font=("Courier",10),
                 bg=BG_PRIMARY, fg=TEXT_SEC).pack(padx=20, anchor="w")
        for i,ev in enumerate(evals):
            clr = SUCCESS if ev>0 else DANGER
            tk.Label(self.stab_frame, text=f"  λ{i+1} = {ev:.4f} GPa",
                     font=("Courier",11,"bold"), bg=BG_PRIMARY, fg=clr
                     ).pack(padx=20, anchor="w", pady=1)
        cond = np.linalg.cond(self.C)
        det  = np.linalg.det(self.C)
        tk.Label(self.stab_frame, text=f"\nCondition Number: {cond:.3e}\nDeterminant: {det:.3e}",
                 font=("Courier",10), bg=BG_PRIMARY, fg=TEXT_SEC).pack(padx=20, anchor="w")

        # Derived tab
        for w in self.derived_frame.winfo_children(): w.destroy()
        S = self.S_mat
        props = [
            ("E₁  Young's Modulus","1/S[0,0]",1/S[0,0],"GPa"),
            ("E₂  Young's Modulus","1/S[1,1]",1/S[1,1],"GPa"),
            ("E₃  Young's Modulus","1/S[2,2]",1/S[2,2],"GPa"),
            ("G₂₃  Shear Modulus","1/S[3,3]",1/S[3,3],"GPa"),
            ("G₁₃  Shear Modulus","1/S[4,4]",1/S[4,4],"GPa"),
            ("G₁₂  Shear Modulus","1/S[5,5]",1/S[5,5],"GPa"),
            ("ν₁₂  Poisson Ratio","-S[1,0]/S[0,0]",-S[1,0]/S[0,0],"—"),
            ("ν₁₃  Poisson Ratio","-S[2,0]/S[0,0]",-S[2,0]/S[0,0],"—"),
            ("ν₂₃  Poisson Ratio","-S[2,1]/S[1,1]",-S[2,1]/S[1,1],"—"),
            ("K  Bulk Modulus","Voigt avg.",(self.C[0,0]+self.C[1,1]+self.C[2,2]+
             2*(self.C[0,1]+self.C[0,2]+self.C[1,2]))/9,"GPa"),
        ]
        g = tk.Frame(self.derived_frame, bg=BG_PRIMARY)
        g.pack(fill="both", padx=20, pady=14)
        for idx,(lbl,_,val,unit) in enumerate(props):
            row=idx//3; col=idx%3
            c2 = tk.Frame(g, bg=BG_CARD, padx=14, pady=10)
            c2.grid(row=row, column=col, padx=5, pady=4, sticky="nsew")
            g.columnconfigure(col, weight=1)
            tk.Label(c2, text=lbl, font=("Helvetica",9), bg=BG_CARD, fg=TEXT_MUTED).pack(anchor="w")
            tk.Label(c2, text=f"{fmt_val(val)}  {unit}", font=("Courier",14,"bold"),
                     bg=BG_CARD, fg=self.color).pack(anchor="w")


# ══════════════════════════════════════════════════════════════════════════════
# 3. MICROMECHANICS MODULE
# ══════════════════════════════════════════════════════════════════════════════

class MicromechanicsPage(BasePage):
    def __init__(self, parent, app):
        super().__init__(parent, app, mid="micro")
        self.page_header("Micromechanics Laboratory",
                         "ROM · Halpin-Tsai · Mori-Tanaka · Self-Consistent")
        self._build()

    def _build(self):
        pw = tk.PanedWindow(self, orient="horizontal", bg=BG_PRIMARY, sashwidth=4)
        pw.pack(fill="both", expand=True)

        left = tk.Frame(pw, bg=BG_SECONDARY, width=320)
        pw.add(left, minsize=280)
        left.pack_propagate(False)
        sf = scrollable_frame(left)

        section_header(sf, "FIBER PROPERTIES")
        fiber_params = [
            ("E₁f [GPa]","Ef1",230.0),("E₂f [GPa]","Ef2",15.0),
            ("ν₁₂f","nuf12",0.20),("G₁₂f [GPa]","Gf12",27.0),
        ]
        self.f_entries = {}
        for lbl,key,defv in fiber_params:
            row = tk.Frame(sf, bg=BG_SECONDARY); row.pack(fill="x", padx=8, pady=2)
            tk.Label(row, text=lbl, font=("Helvetica",9), bg=BG_SECONDARY,
                     fg=TEXT_SEC, anchor="w", width=18).pack(side="left")
            e = make_entry(row, defv, width=9); e.pack(side="right"); self.f_entries[key]=e

        section_header(sf, "MATRIX PROPERTIES")
        matrix_params = [("Eₘ [GPa]","Em",4.2),("νₘ","num",0.34)]
        self.m_entries = {}
        for lbl,key,defv in matrix_params:
            row = tk.Frame(sf, bg=BG_SECONDARY); row.pack(fill="x", padx=8, pady=2)
            tk.Label(row, text=lbl, font=("Helvetica",9), bg=BG_SECONDARY,
                     fg=TEXT_SEC, anchor="w", width=18).pack(side="left")
            e = make_entry(row, defv, width=9); e.pack(side="right"); self.m_entries[key]=e

        section_header(sf, "VOLUME FRACTION")
        row = tk.Frame(sf, bg=BG_SECONDARY); row.pack(fill="x", padx=8, pady=2)
        tk.Label(row, text="Vf (0–1)", font=("Helvetica",9), bg=BG_SECONDARY,
                 fg=TEXT_SEC, anchor="w", width=18).pack(side="left")
        self.vf_entry = make_entry(row, 0.60, width=9); self.vf_entry.pack(side="right")

        tk.Frame(sf, bg=BORDER, height=1).pack(fill="x", padx=8, pady=8)
        accent_btn(sf, "▶  COMPUTE ALL MODELS", self._compute, bg=self.color, fg="#000").pack(fill="x", padx=8, pady=4)
        accent_btn(sf, "📈  PLOT Vf SWEEP", self._plot_sweep, bg=self.dim, fg=self.color).pack(fill="x", padx=8, pady=2)

        # RIGHT
        right = tk.Frame(pw, bg=BG_PRIMARY)
        pw.add(right, minsize=500)

        self.results_frame = tk.Frame(right, bg=BG_PRIMARY)
        self.results_frame.pack(fill="both", expand=True, padx=16, pady=10)

        self.fig_frame = tk.Frame(right, bg=BG_PRIMARY)
        self.fig_frame.pack(fill="both", expand=True, padx=16, pady=4)

    def _get_params(self):
        fp = {k:float(e.get()) for k,e in self.f_entries.items()}
        mp = {k:float(e.get()) for k,e in self.m_entries.items()}
        Vf = float(self.vf_entry.get())
        return fp,mp,Vf

    def _compute(self):
        try:
            fp,mp,Vf = self._get_params()
        except: messagebox.showerror("Error","Invalid inputs"); return

        models = {
            "Rule of Mixtures": rule_of_mixtures(fp["Ef1"],fp["Ef2"],fp["nuf12"],fp["Gf12"],mp["Em"],mp["num"],Vf),
            "Halpin-Tsai":      halpin_tsai(fp["Ef1"],fp["Ef2"],fp["nuf12"],fp["Gf12"],mp["Em"],mp["num"],Vf),
            "Mori-Tanaka":      mori_tanaka(fp["Ef1"],fp["Ef2"],fp["nuf12"],fp["Gf12"],mp["Em"],mp["num"],Vf),
            "Self-Consistent":  self_consistent(fp["Ef1"],fp["Ef2"],fp["nuf12"],fp["Gf12"],mp["Em"],mp["num"],Vf),
        }

        for w in self.results_frame.winfo_children(): w.destroy()
        tk.Label(self.results_frame, text=f"Results  (Vf = {Vf:.2f})",
                 font=("Helvetica",11,"bold"), bg=BG_PRIMARY, fg=TEXT_PRIMARY
                 ).pack(anchor="w", pady=(0,6))

        cols = list(models.keys())
        keys = ["E1","E2","nu12","G12"]
        units= ["GPa","GPa","—","GPa"]
        labels=["E₁","E₂","ν₁₂","G₁₂"]

        g = tk.Frame(self.results_frame, bg=BG_PRIMARY)
        g.pack(fill="x")
        # Header
        tk.Label(g, text="Property", font=("Courier",9,"bold"),
                 bg=BG_PANEL, fg=TEXT_MUTED, padx=10, pady=4, width=14
                 ).grid(row=0, column=0, padx=2, pady=1, sticky="nsew")
        for j,model in enumerate(cols):
            tk.Label(g, text=model, font=("Courier",8,"bold"),
                     bg=BG_PANEL, fg=ACCENT2, padx=6, pady=4, width=16, wraplength=120
                     ).grid(row=0, column=j+1, padx=2, pady=1, sticky="nsew")

        for i,(k,lbl,unit) in enumerate(zip(keys,labels,units)):
            tk.Label(g, text=f"{lbl} [{unit}]", font=("Courier",9),
                     bg=BG_CARD, fg=TEXT_SEC, padx=10, pady=4, width=14
                     ).grid(row=i+1, column=0, padx=2, pady=1, sticky="nsew")
            vals = [models[m].get(k,0) for m in cols]
            max_v = max(abs(v) for v in vals) or 1
            for j,v in enumerate(vals):
                intensity = abs(v)/max_v
                clr = ACCENT if intensity > 0.7 else (ACCENT2 if intensity > 0.3 else TEXT_SEC)
                tk.Label(g, text=fmt_val(v), font=("Courier",10,"bold"),
                         bg=BG_CARD, fg=clr, padx=6, pady=4, width=16
                         ).grid(row=i+1, column=j+1, padx=2, pady=1, sticky="nsew")

    def _plot_sweep(self):
        try:
            fp,mp,_ = self._get_params()
        except: messagebox.showerror("Error","Invalid inputs"); return

        for w in self.fig_frame.winfo_children(): w.destroy()

        Vf_arr = np.linspace(0.01,0.80,80)
        data = {m:{"E1":[],"E2":[],"G12":[]} for m in ["ROM","HT","MT","SC"]}
        funcs = [rule_of_mixtures, halpin_tsai, mori_tanaka, self_consistent]
        keys2 = ["ROM","HT","MT","SC"]
        for Vf in Vf_arr:
            for fn,key in zip(funcs,keys2):
                r = fn(fp["Ef1"],fp["Ef2"],fp["nuf12"],fp["Gf12"],mp["Em"],mp["num"],Vf)
                data[key]["E1"].append(r["E1"])
                data[key]["E2"].append(r["E2"])
                data[key]["G12"].append(r["G12"])

        fig = Figure(figsize=(10,3.5), dpi=100, facecolor=PLT_BG)
        colors = [ACCENT, ACCENT2, SUCCESS, DANGER]
        for idx,(prop,ylabel,title) in enumerate([("E1","E₁ (GPa)","Longitudinal Modulus"),
                                                   ("E2","E₂ (GPa)","Transverse Modulus"),
                                                   ("G12","G₁₂ (GPa)","Shear Modulus")]):
            ax = fig.add_subplot(1,3,idx+1, facecolor=PLT_AX_BG)
            style_ax(ax)
            for k,clr in zip(keys2,colors):
                ax.plot(Vf_arr, data[k][prop], color=clr, lw=1.8, label=k)
            ax.set_xlabel("Fiber Volume Fraction Vf", fontsize=8)
            ax.set_ylabel(ylabel, fontsize=8)
            ax.set_title(title, fontsize=9, color=TEXT_PRIMARY)
            ax.legend(fontsize=7)
        fig.tight_layout(pad=1.5)
        self.embed_figure(self.fig_frame, fig, toolbar=False)


# ══════════════════════════════════════════════════════════════════════════════
# 4. LAMINATE ANALYSIS MODULE  (Classical Lamination Theory)
# ══════════════════════════════════════════════════════════════════════════════

class LaminatePage(BasePage):
    def __init__(self, parent, app):
        super().__init__(parent, app, mid="laminate")
        self.plies = []
        self.page_header("Laminate Analysis (CLT)", "ABD Matrix · Engineering Constants · Ply-by-Ply Stresses")
        self._build()

    def _build(self):
        pw = tk.PanedWindow(self, orient="horizontal", bg=BG_PRIMARY, sashwidth=4)
        pw.pack(fill="both", expand=True)

        # LEFT
        left = tk.Frame(pw, bg=BG_SECONDARY, width=340)
        pw.add(left, minsize=300); left.pack_propagate(False)
        sf = scrollable_frame(left)

        section_header(sf, "ADD PLY")
        ply_fields = [("E₁ [GPa]","E1",130.0),("E₂ [GPa]","E2",9.0),
                      ("ν₁₂","nu12",0.28),("G₁₂ [GPa]","G12",6.9),
                      ("θ  [deg]","theta",0.0),("t  [mm]","t",0.25)]
        self.ply_entries = {}
        for lbl,key,defv in ply_fields:
            row = tk.Frame(sf, bg=BG_SECONDARY); row.pack(fill="x", padx=8, pady=2)
            tk.Label(row, text=lbl, font=("Helvetica",9), bg=BG_SECONDARY,
                     fg=TEXT_SEC, anchor="w", width=14).pack(side="left")
            e = make_entry(row, defv, width=9); e.pack(side="right"); self.ply_entries[key]=e

        accent_btn(sf,"+ Add Ply",self._add_ply, bg=self.color, fg="#000").pack(fill="x",padx=8,pady=4)

        section_header(sf, "STACKING SEQUENCE")
        self.ply_list = tk.Listbox(sf, bg=BG_PANEL, fg=TEXT_PRIMARY,
                                   selectbackground=self.dim, font=("Courier",9),
                                   relief="flat", height=10, activestyle="none",
                                   highlightthickness=0)
        self.ply_list.pack(fill="x", padx=8, pady=2)

        btn_row = tk.Frame(sf, bg=BG_SECONDARY); btn_row.pack(fill="x",padx=8,pady=2)
        tk.Button(btn_row, text="Remove Selected", font=("Courier",9),
                  bg=DANGER, fg="#fff", relief="flat", cursor="hand2",
                  command=self._remove_ply).pack(side="left", padx=2)
        tk.Button(btn_row, text="Clear All", font=("Courier",9),
                  bg=BG_CARD, fg=TEXT_SEC, relief="flat", cursor="hand2",
                  command=self._clear_plies).pack(side="left", padx=2)

        section_header(sf, "LOADING (N/m, N·m/m)")
        load_fields = [("Nx [N/m]","Nx",0.0),("Ny [N/m]","Ny",0.0),
                       ("Nxy [N/m]","Nxy",0.0),("Mx [N·m/m]","Mx",0.0),
                       ("My [N·m/m]","My",0.0),("Mxy [N·m/m]","Mxy",0.0)]
        self.load_entries = {}
        for lbl,key,defv in load_fields:
            row = tk.Frame(sf, bg=BG_SECONDARY); row.pack(fill="x", padx=8, pady=2)
            tk.Label(row, text=lbl, font=("Helvetica",9), bg=BG_SECONDARY,
                     fg=TEXT_SEC, anchor="w", width=16).pack(side="left")
            e = make_entry(row, defv, width=9); e.pack(side="right"); self.load_entries[key]=e

        tk.Frame(sf, bg=BORDER, height=1).pack(fill="x", padx=8, pady=8)
        accent_btn(sf,"▶  ANALYZE LAMINATE", self._analyze, bg=self.color, fg="#000").pack(fill="x",padx=8,pady=4)

        # RIGHT
        right = tk.Frame(pw, bg=BG_PRIMARY)
        pw.add(right, minsize=550)

        tab_bar = tk.Frame(right, bg=BG_SECONDARY, height=40)
        tab_bar.pack(fill="x"); tab_bar.pack_propagate(False)
        tk.Frame(right, bg=BORDER, height=1).pack(fill="x")
        self.tab_area = tk.Frame(right, bg=BG_PRIMARY)
        self.tab_area.pack(fill="both", expand=True)

        self.ltabs = {}; self.ltab_btns = {}
        for key,lbl in [("ABD","ABD Matrix"),("eng","Eng. Constants"),
                        ("ply","Ply Stresses"),("polar","Polar Stiffness")]:
            btn = tk.Button(tab_bar, text=lbl, font=("Courier",9),
                            bg=BG_SECONDARY, fg=TEXT_MUTED, relief="flat",
                            cursor="hand2", padx=12, pady=10,
                            command=lambda k=key: self._ltab(k))
            btn.pack(side="left")
            self.ltab_btns[key] = btn
            self.ltabs[key] = tk.Frame(self.tab_area, bg=BG_PRIMARY)

        self._ltab("ABD")

        # Pre-load a default laminate  [0/90/±45]s
        defaults = [(0,0.25),(90,0.25),(45,0.25),(-45,0.25),
                    (-45,0.25),(45,0.25),(90,0.25),(0,0.25)]
        for th,t in defaults:
            self.plies.append({"E1":130.0,"E2":9.0,"nu12":0.28,"G12":6.9,"theta":th,"t":t})
            self.ply_list.insert("end", f"  θ={th:+6.1f}°  t={t} mm")

    def _ltab(self, key):
        for k,f in self.ltabs.items(): f.pack_forget()
        self.ltabs[key].pack(fill="both", expand=True)
        for k,b in self.ltab_btns.items():
            b.config(fg=self.color if k==key else TEXT_MUTED,
                     bg=BG_PRIMARY if k==key else BG_SECONDARY)

    def _add_ply(self):
        try:
            p = {k:float(e.get()) for k,e in self.ply_entries.items()}
        except: messagebox.showerror("Error","Invalid ply properties"); return
        self.plies.append(p)
        self.ply_list.insert("end",
            f"  θ={p['theta']:+6.1f}°  t={p['t']} mm  E₁={p['E1']} GPa")

    def _remove_ply(self):
        sel = self.ply_list.curselection()
        if not sel: return
        idx = sel[0]; self.plies.pop(idx); self.ply_list.delete(idx)

    def _clear_plies(self):
        self.plies.clear(); self.ply_list.delete(0,"end")

    def _analyze(self):
        if len(self.plies) < 1:
            messagebox.showwarning("CMVL","Add at least one ply."); return
        A,B,D = ABD_matrix(self.plies)
        t_total = sum(p["t"] for p in self.plies)
        eng = laminate_engineering_constants(A,B,D,t_total)

        # ABD tab
        for w in self.ltabs["ABD"].winfo_children(): w.destroy()
        ABD = np.block([[A,B],[B,D]])
        lbl6 = ["N₁","N₂","N₆","M₁","M₂","M₆"]
        eps6 = ["ε₁⁰","ε₂⁰","ε₆⁰","κ₁","κ₂","κ₆"]
        mf = tk.Frame(self.ltabs["ABD"], bg=BG_PRIMARY)
        mf.pack(padx=18, pady=14)
        matrix_grid(mf, ABD, lbl6, eps6, "ABD Matrix  (N/m, N, N·m)", accent=self.color)

        # Engineering constants tab
        for w in self.ltabs["eng"].winfo_children(): w.destroy()
        g = tk.Frame(self.ltabs["eng"], bg=BG_PRIMARY); g.pack(padx=20,pady=16,fill="x")
        tk.Label(g, text=f"Laminate Engineering Constants  (total t = {t_total:.3f} mm)",
                 font=("Helvetica",11,"bold"), bg=BG_PRIMARY, fg=TEXT_PRIMARY
                 ).pack(anchor="w", pady=(0,10))
        cols_labels = [("Ex","Longitudinal Modulus"),("Ey","Transverse Modulus"),
                       ("Gxy","In-Plane Shear Modulus"),("nuxy","Poisson Ratio νxy")]
        units = ["GPa","GPa","GPa","—"]
        row_f = tk.Frame(g, bg=BG_PRIMARY); row_f.pack(fill="x")
        for idx,(k,lbl) in enumerate(cols_labels):
            c2 = tk.Frame(row_f, bg=BG_CARD, padx=14, pady=12)
            c2.pack(side="left", expand=True, fill="both", padx=4)
            tk.Label(c2, text=lbl, font=("Helvetica",9), bg=BG_CARD, fg=TEXT_MUTED
                     ).pack(anchor="w")
            val = eng.get(k, 0.0)
            tk.Label(c2, text=f"{fmt_val(val)}  {units[idx]}",
                     font=("Courier",14,"bold"), bg=BG_CARD, fg=self.color
                     ).pack(anchor="w", pady=(4,0))

        # A B D sub-matrices
        sub_f = tk.Frame(g, bg=BG_PRIMARY); sub_f.pack(fill="x", pady=10)
        for mat,title in [(A,"A Matrix  (extensional)"),(B,"B Matrix  (coupling)"),(D,"D Matrix  (bending)")]:
            sf2 = tk.Frame(sub_f, bg=BG_PRIMARY); sf2.pack(side="left", expand=True, fill="both", padx=6)
            lbl3 = ["N₁","N₂","N₆"]
            matrix_grid(sf2, mat, lbl3, lbl3, title)

        # Ply stresses tab
        for w in self.ltabs["ply"].winfo_children(): w.destroy()
        loads = np.array([float(self.load_entries[k].get()) for k in
                          ["Nx","Ny","Nxy","Mx","My","Mxy"]])
        non_zero = any(abs(v) > 0 for v in loads)
        if not non_zero:
            tk.Label(self.ltabs["ply"], text="Enter loading and re-analyze for ply stresses.",
                     font=("Courier",10), bg=BG_PRIMARY, fg=TEXT_MUTED
                     ).pack(pady=30)
        else:
            self._compute_ply_stresses(loads, A, B, D)

        # Polar stiffness tab
        self._plot_polar_stiffness()

    def _compute_ply_stresses(self, loads, A, B, D):
        try:
            ABD_full = np.block([[A,B],[B,D]])
            abd = inv(ABD_full)
            deform = abd @ loads
            eps0 = deform[:3]; kappa = deform[3:]

            t_total = sum(p["t"]*1e-3 for p in self.plies)
            z = -t_total/2.0
            rows = []
            for i,ply in enumerate(self.plies):
                t_m = ply["t"]*1e-3
                z_mid = z + t_m/2.0
                eps_glob = eps0 + z_mid*kappa
                # Transform to local
                th = ply["theta"]
                eps_loc = stress_transform(eps_glob, th)  # approximate
                Qb = Qbar_matrix(Q_matrix(ply["E1"],ply["E2"],ply["nu12"],ply["G12"]), th)
                sig_glob = Qb @ eps_glob
                sig_loc = stress_transform(sig_glob, th)
                rows.append((i+1, th, sig_loc[0], sig_loc[1], sig_loc[2]))
                z += t_m

            # Draw table
            headers = ["Ply","θ (°)","σ₁ (MPa)","σ₂ (MPa)","τ₁₂ (MPa)"]
            g = tk.Frame(self.ltabs["ply"], bg=BG_PRIMARY); g.pack(padx=16, pady=10, fill="x")
            for j,h in enumerate(headers):
                tk.Label(g, text=h, font=("Courier",9,"bold"), bg=BG_PANEL,
                         fg=TEXT_MUTED, padx=8, pady=4, width=12, anchor="center"
                         ).grid(row=0, column=j, padx=2, pady=1)
            for i,(ply_n, th, s1, s2, s12) in enumerate(rows):
                vals = [str(ply_n), f"{th:.1f}", f"{s1*1e-6:.2f}", f"{s2*1e-6:.2f}", f"{s12*1e-6:.2f}"]
                for j,v in enumerate(vals):
                    tk.Label(g, text=v, font=("Courier",10), bg=BG_CARD,
                             fg=TEXT_PRIMARY, padx=8, pady=4, width=12, anchor="center"
                             ).grid(row=i+1, column=j, padx=2, pady=1)
        except Exception as e:
            tk.Label(self.ltabs["ply"], text=f"Error: {e}", bg=BG_PRIMARY, fg=DANGER,
                     font=("Courier",10)).pack(pady=20)

    def _plot_polar_stiffness(self):
        for w in self.ltabs["polar"].winfo_children(): w.destroy()
        fig = Figure(figsize=(5,5), dpi=100, facecolor=PLT_BG)
        ax = fig.add_subplot(111, projection="polar", facecolor=PLT_AX_BG)
        ax.tick_params(colors=PLT_TEXT, labelsize=7)
        for sp in ax.spines.values(): sp.set_color(PLT_GRID)

        thetas = np.linspace(0, 2*np.pi, 180)
        Ex_vals = []
        for th in thetas:
            th_deg = math.degrees(th)
            E1,E2,nu12,G12 = 130.0, 9.0, 0.28, 6.9
            nu21 = nu12*E2/E1
            t = th
            c = math.cos(t); s = math.sin(t)
            Q = Q_matrix(E1,E2,nu12,G12)
            Qb = Qbar_matrix(Q, th_deg)
            # Ex = 1/a11 for single ply
            S11 = Qb[0,0]
            Ex_vals.append(S11/1e3)  # in GPa units approx
        ax.plot(thetas, Ex_vals, color=self.color, lw=2)
        ax.fill(thetas, Ex_vals, alpha=0.15, color=self.color)
        ax.set_title("Directional Stiffness E(θ) [GPa]", color=TEXT_PRIMARY, fontsize=9, pad=12)
        ax.grid(color=PLT_GRID, lw=0.5)
        self.embed_figure(self.ltabs["polar"], fig, toolbar=False)


# ══════════════════════════════════════════════════════════════════════════════
# 5. FAILURE ANALYSIS MODULE
# ══════════════════════════════════════════════════════════════════════════════

class FailurePage(BasePage):
    def __init__(self, parent, app):
        super().__init__(parent, app, mid="failure")
        self.page_header("Failure Analysis Laboratory",
                         "Max Stress · Tsai-Hill · Tsai-Wu · Hashin")
        self._build()

    def _build(self):
        pw = tk.PanedWindow(self, orient="horizontal", bg=BG_PRIMARY, sashwidth=4)
        pw.pack(fill="both", expand=True)

        left = tk.Frame(pw, bg=BG_SECONDARY, width=320)
        pw.add(left, minsize=280); left.pack_propagate(False)
        sf = scrollable_frame(left)

        section_header(sf, "APPLIED STRESSES (MPa)")
        stress_fields = [("σ₁","s1",500.0),("σ₂","s2",-50.0),("τ₁₂","s12",30.0)]
        self.stress_entries = {}
        for lbl,key,defv in stress_fields:
            row = tk.Frame(sf, bg=BG_SECONDARY); row.pack(fill="x", padx=8, pady=2)
            tk.Label(row, text=lbl, font=("Helvetica",10), bg=BG_SECONDARY,
                     fg=TEXT_SEC, anchor="w", width=8).pack(side="left")
            e = make_entry(row, defv, width=10); e.pack(side="right"); self.stress_entries[key]=e

        section_header(sf, "STRENGTH PROPERTIES (MPa)")
        strength_fields = [
            ("Xₜ (tensile long.)","Xt",1500.0),("Xc (comp. long.)","Xc",1200.0),
            ("Yₜ (tensile trans.)","Yt",50.0),("Yc (comp. trans.)","Yc",200.0),
            ("S₁₂ (in-plane shear)","S12",90.0),
        ]
        self.str_entries = {}
        for lbl,key,defv in strength_fields:
            row = tk.Frame(sf, bg=BG_SECONDARY); row.pack(fill="x", padx=8, pady=2)
            tk.Label(row, text=lbl, font=("Helvetica",9), bg=BG_SECONDARY,
                     fg=TEXT_SEC, anchor="w", width=22).pack(side="left")
            e = make_entry(row, defv, width=9); e.pack(side="right"); self.str_entries[key]=e

        tk.Frame(sf, bg=BORDER, height=1).pack(fill="x", padx=8, pady=8)
        accent_btn(sf, "▶  COMPUTE FAILURE INDICES", self._compute, bg=self.color, fg="#000").pack(fill="x",padx=8,pady=4)
        accent_btn(sf, "📈  PLOT FAILURE ENVELOPE", self._plot_envelope, bg=self.dim, fg=self.color
                   ).pack(fill="x",padx=8,pady=2)

        right = tk.Frame(pw, bg=BG_PRIMARY)
        pw.add(right, minsize=500)

        self.results_frame = tk.Frame(right, bg=BG_PRIMARY)
        self.results_frame.pack(fill="x", padx=16, pady=10)
        self.plot_frame = tk.Frame(right, bg=BG_PRIMARY)
        self.plot_frame.pack(fill="both", expand=True, padx=16)

    def _get_stresses(self):
        s1  = float(self.stress_entries["s1"].get())
        s2  = float(self.stress_entries["s2"].get())
        s12 = float(self.stress_entries["s12"].get())
        Xt  = float(self.str_entries["Xt"].get())
        Xc  = float(self.str_entries["Xc"].get())
        Yt  = float(self.str_entries["Yt"].get())
        Yc  = float(self.str_entries["Yc"].get())
        S12 = float(self.str_entries["S12"].get())
        return s1,s2,s12,Xt,Xc,Yt,Yc,S12

    def _compute(self):
        try:
            s1,s2,s12,Xt,Xc,Yt,Yc,S12 = self._get_stresses()
        except: messagebox.showerror("Error","Invalid inputs"); return

        results = {
            "Maximum Stress":  failure_max_stress(s1,s2,s12,Xt,Xc,Yt,Yc,S12),
            "Tsai-Hill":       failure_tsai_hill(s1,s2,s12,Xt,Xc,Yt,Yc,S12),
            "Tsai-Wu":         failure_tsai_wu(s1,s2,s12,Xt,Xc,Yt,Yc,S12),
            "Hashin":          failure_hashin(s1,s2,s12,Xt,Xc,Yt,Yc,S12),
        }

        for w in self.results_frame.winfo_children(): w.destroy()
        tk.Label(self.results_frame, text="Failure Indices  (FI ≥ 1.0 → FAILURE)",
                 font=("Helvetica",11,"bold"), bg=BG_PRIMARY, fg=TEXT_PRIMARY
                 ).pack(anchor="w", pady=(0,8))

        g = tk.Frame(self.results_frame, bg=BG_PRIMARY); g.pack(fill="x")
        for j,(name,fi) in enumerate(results.items()):
            col = j%2; row = j//2
            card = tk.Frame(g, bg=BG_CARD, padx=16, pady=12)
            card.grid(row=row, column=col, padx=5, pady=4, sticky="nsew")
            g.columnconfigure(col, weight=1)
            failed = fi >= 1.0
            status = "✗  FAILED" if failed else "✔  SAFE"
            color = DANGER if failed else SUCCESS
            RF = 1.0/fi if fi > 0 else float("inf")
            tk.Label(card, text=name, font=("Helvetica",10,"bold"),
                     bg=BG_CARD, fg=TEXT_PRIMARY).pack(anchor="w")
            tk.Label(card, text=f"FI = {fi:.4f}", font=("Courier",14,"bold"),
                     bg=BG_CARD, fg=color).pack(anchor="w", pady=(4,0))
            tk.Label(card, text=f"Reserve Factor = {RF:.2f}  ·  {status}",
                     font=("Courier",9), bg=BG_CARD, fg=color).pack(anchor="w")

    def _plot_envelope(self):
        try:
            _,_,_,Xt,Xc,Yt,Yc,S12 = self._get_stresses()
        except: return

        for w in self.plot_frame.winfo_children(): w.destroy()

        fig = Figure(figsize=(9,4), dpi=100, facecolor=PLT_BG)
        s12_fixed = float(self.stress_entries["s12"].get())

        # Tsai-Wu and Tsai-Hill envelopes in σ1-σ2 space
        n = 120
        s1_arr = np.linspace(-Xc*1.1, Xt*1.1, n)

        ax1 = fig.add_subplot(121, facecolor=PLT_AX_BG)
        ax2 = fig.add_subplot(122, facecolor=PLT_AX_BG)
        style_ax(ax1); style_ax(ax2)

        for ax, crit, name, color in [
            (ax1, "hill",  "Tsai-Hill", ACCENT),
            (ax2, "wu",    "Tsai-Wu",   ACCENT2)
        ]:
            envelope_s2_pos = []; envelope_s2_neg = []
            for s1 in s1_arr:
                # Solve for s2 limits
                if crit == "hill":
                    X = Xt if s1>=0 else Xc
                    # (s1/X)^2 - s1*s2/X^2 + (s2/Y)^2 + (s12/S12)^2 = 1
                    # quadratic in s2 for each s1
                    Y = Yt
                    A_c = 1/Y**2; B_c = -s1/X**2; C_c = (s1/X)**2 + (s12_fixed/S12)**2 - 1
                    disc = B_c**2 - 4*A_c*C_c
                    if disc >= 0:
                        envelope_s2_pos.append((-B_c+math.sqrt(disc))/(2*A_c))
                        envelope_s2_neg.append((-B_c-math.sqrt(disc))/(2*A_c))
                    else:
                        envelope_s2_pos.append(np.nan); envelope_s2_neg.append(np.nan)
                else:  # Tsai-Wu
                    F1=1/Xt-1/Xc; F2=1/Yt-1/Yc
                    F11=1/(Xt*Xc); F22=1/(Yt*Yc); F66=1/S12**2
                    F12=-0.5*math.sqrt(F11*F22)
                    A_c=F22; B_c=F2+2*F12*s1; C_c=F1*s1+F11*s1**2+F66*s12_fixed**2-1
                    disc=B_c**2-4*A_c*C_c
                    if disc>=0:
                        envelope_s2_pos.append((-B_c+math.sqrt(disc))/(2*A_c))
                        envelope_s2_neg.append((-B_c-math.sqrt(disc))/(2*A_c))
                    else:
                        envelope_s2_pos.append(np.nan); envelope_s2_neg.append(np.nan)

            ax.plot(s1_arr, envelope_s2_pos, color=color, lw=2, label="Failure boundary")
            ax.plot(s1_arr, envelope_s2_neg, color=color, lw=2)
            ax.fill_between(s1_arr, envelope_s2_neg, envelope_s2_pos,
                            alpha=0.12, color=color, label="Safe zone")
            s1_pt = float(self.stress_entries["s1"].get())
            s2_pt = float(self.stress_entries["s2"].get())
            ax.scatter([s1_pt],[s2_pt], color=WARNING, s=80, zorder=5, label="Stress state")
            ax.axhline(0, color=PLT_GRID, lw=0.5); ax.axvline(0, color=PLT_GRID, lw=0.5)
            ax.set_xlabel("σ₁ (MPa)", fontsize=8); ax.set_ylabel("σ₂ (MPa)", fontsize=8)
            ax.set_title(f"{name} Envelope  (τ₁₂={s12_fixed} MPa)", fontsize=9, color=TEXT_PRIMARY)
            ax.legend(fontsize=7)

        fig.tight_layout(pad=1.5)
        self.embed_figure(self.plot_frame, fig, toolbar=False)


# ══════════════════════════════════════════════════════════════════════════════
# 6. STRESS-STRAIN SIMULATOR
# ══════════════════════════════════════════════════════════════════════════════

class StressStrainPage(BasePage):
    def __init__(self, parent, app):
        super().__init__(parent, app, mid="stress_strain")
        self.page_header("Stress-Strain Simulator",
                         "Constitutive Response · Principal Values · Strain Energy")
        self._build()

    def _build(self):
        pw = tk.PanedWindow(self, orient="horizontal", bg=BG_PRIMARY, sashwidth=4)
        pw.pack(fill="both", expand=True)

        left = tk.Frame(pw, bg=BG_SECONDARY, width=300)
        pw.add(left, minsize=260); left.pack_propagate(False)
        sf = scrollable_frame(left)

        section_header(sf, "MATERIAL")
        mat_params = [("E₁ [GPa]","E1",130.0),("E₂ [GPa]","E2",9.0),
                      ("E₃ [GPa]","E3",9.0),("ν₁₂","nu12",0.28),
                      ("ν₁₃","nu13",0.28),("ν₂₃","nu23",0.40),
                      ("G₁₂ [GPa]","G12",6.9),("G₁₃ [GPa]","G13",6.9),
                      ("G₂₃ [GPa]","G23",4.0)]
        self.mat_entries = {}
        for lbl,key,defv in mat_params:
            row = tk.Frame(sf, bg=BG_SECONDARY); row.pack(fill="x",padx=8,pady=2)
            tk.Label(row,text=lbl,font=("Helvetica",9),bg=BG_SECONDARY,
                     fg=TEXT_SEC,anchor="w",width=16).pack(side="left")
            e = make_entry(row,defv,width=9); e.pack(side="right"); self.mat_entries[key]=e

        section_header(sf, "APPLIED STRESS STATE (MPa)")
        stress_fields = [("σ₁","s1",100.0),("σ₂","s2",0.0),("σ₃","s3",0.0),
                         ("τ₂₃","s23",0.0),("τ₁₃","s13",0.0),("τ₁₂","s12",20.0)]
        self.ss_entries = {}
        for lbl,key,defv in stress_fields:
            row = tk.Frame(sf, bg=BG_SECONDARY); row.pack(fill="x",padx=8,pady=2)
            tk.Label(row,text=lbl,font=("Helvetica",9),bg=BG_SECONDARY,
                     fg=TEXT_SEC,anchor="w",width=8).pack(side="left")
            e = make_entry(row,defv,width=10); e.pack(side="right"); self.ss_entries[key]=e

        tk.Frame(sf, bg=BORDER, height=1).pack(fill="x", padx=8, pady=8)
        accent_btn(sf, "▶  COMPUTE RESPONSE", self._compute, bg=self.color, fg="#000").pack(fill="x",padx=8,pady=4)
        accent_btn(sf, "📈  PARAMETRIC SWEEP", self._sweep, bg=self.dim, fg=self.color
                   ).pack(fill="x",padx=8,pady=2)

        right = tk.Frame(pw, bg=BG_PRIMARY)
        pw.add(right, minsize=550)
        self.result_frame = tk.Frame(right, bg=BG_PRIMARY)
        self.result_frame.pack(fill="x", padx=16, pady=10)
        self.fig_frame = tk.Frame(right, bg=BG_PRIMARY)
        self.fig_frame.pack(fill="both", expand=True, padx=16, pady=4)

    def _compute(self):
        try:
            p = {k:float(e.get()) for k,e in self.mat_entries.items()}
            sigma = np.array([float(self.ss_entries[k].get()) for k in
                              ["s1","s2","s3","s23","s13","s12"]])  # MPa
        except: messagebox.showerror("Error","Invalid inputs"); return

        S = compliance_from_engineering(p["E1"],p["E2"],p["E3"],p["nu12"],p["nu13"],
                                        p["nu23"],p["G12"],p["G13"],p["G23"])
        S_scaled = S / 1e3  # convert GPa-based S to MPa: S[1/GPa]*sigma[GPa]=eps
        sigma_gpa = sigma / 1e3  # to GPa
        eps = S_scaled @ sigma_gpa  # dimensionless

        # Principal stresses
        S3x3 = np.array([[sigma[0],sigma[5],sigma[4]],
                          [sigma[5],sigma[1],sigma[3]],
                          [sigma[4],sigma[3],sigma[2]]])
        princ,_ = np.linalg.eigh(S3x3)
        von_mises = math.sqrt(0.5*((sigma[0]-sigma[1])**2 +
                                    (sigma[1]-sigma[2])**2 +
                                    (sigma[2]-sigma[0])**2 +
                                    6*(sigma[3]**2+sigma[4]**2+sigma[5]**2)))
        W = 0.5 * float(sigma_gpa @ (S / 1e3) @ sigma_gpa)  # MJ/m³ approx

        for w in self.result_frame.winfo_children(): w.destroy()
        g = tk.Frame(self.result_frame, bg=BG_PRIMARY); g.pack(fill="x")

        items = [
            ("ε₁  Strain","",f"{eps[0]*1e6:.2f} με"),
            ("ε₂  Strain","",f"{eps[1]*1e6:.2f} με"),
            ("ε₃  Strain","",f"{eps[2]*1e6:.2f} με"),
            ("γ₂₃  Shear Strain","",f"{eps[3]*1e6:.2f} με"),
            ("Von Mises  σ_vm","",f"{von_mises:.2f} MPa"),
            ("Strain Energy  W","",f"{W:.6f} GPa"),
            ("σ_p1  Principal","",f"{princ[2]:.2f} MPa"),
            ("σ_p2  Principal","",f"{princ[1]:.2f} MPa"),
            ("σ_p3  Principal","",f"{princ[0]:.2f} MPa"),
        ]
        for idx,(lbl,_,val) in enumerate(items):
            col=idx%3; row=idx//3
            c2 = tk.Frame(g, bg=BG_CARD, padx=12, pady=8)
            c2.grid(row=row,column=col,padx=4,pady=3,sticky="nsew")
            g.columnconfigure(col,weight=1)
            tk.Label(c2,text=lbl,font=("Helvetica",8),bg=BG_CARD,fg=TEXT_MUTED).pack(anchor="w")
            tk.Label(c2,text=val,font=("Courier",12,"bold"),bg=BG_CARD,fg=self.color).pack(anchor="w")

    def _sweep(self):
        try:
            p = {k:float(e.get()) for k,e in self.mat_entries.items()}
        except: return
        for w in self.fig_frame.winfo_children(): w.destroy()

        S = compliance_from_engineering(p["E1"],p["E2"],p["E3"],p["nu12"],p["nu13"],
                                        p["nu23"],p["G12"],p["G13"],p["G23"])
        sigma_range = np.linspace(0, 1500, 200)  # MPa
        eps_dirs = []
        for i in range(3):
            eps_i = []
            for s in sigma_range:
                sigma_v = np.zeros(6); sigma_v[i] = s/1e3
                eps_v = S/1e3 @ sigma_v
                eps_i.append(eps_v[i]*1e6)  # με
            eps_dirs.append(eps_i)

        fig = Figure(figsize=(9,3.5), dpi=100, facecolor=PLT_BG)
        colors_dir = [ACCENT, ACCENT2, SUCCESS]
        labels_dir = ["Direction 1 (E₁)","Direction 2 (E₂)","Direction 3 (E₃)"]
        ax = fig.add_subplot(111, facecolor=PLT_AX_BG)
        style_ax(ax)
        for i in range(3):
            ax.plot(eps_dirs[i], sigma_range, color=colors_dir[i], lw=2, label=labels_dir[i])
        ax.set_xlabel("Strain (με)", fontsize=9)
        ax.set_ylabel("Stress σ (MPa)", fontsize=9)
        ax.set_title("Stress-Strain Response (Linear Elastic)", fontsize=10, color=TEXT_PRIMARY)
        ax.legend(fontsize=8)
        fig.tight_layout(pad=1.5)
        self.embed_figure(self.fig_frame, fig, toolbar=False)


# ══════════════════════════════════════════════════════════════════════════════
# 7. OPTIMIZATION LABORATORY
# ══════════════════════════════════════════════════════════════════════════════

class OptimizationPage(BasePage):
    def __init__(self, parent, app):
        super().__init__(parent, app, mid="optimization")
        self.history = []
        self.page_header("Optimization Laboratory",
                         "Maximize Stiffness · Minimize Weight · Fiber Orientation")
        self._build()

    def _build(self):
        pw = tk.PanedWindow(self, orient="horizontal", bg=BG_PRIMARY, sashwidth=4)
        pw.pack(fill="both", expand=True)

        left = tk.Frame(pw, bg=BG_SECONDARY, width=300)
        pw.add(left, minsize=260); left.pack_propagate(False)
        sf = scrollable_frame(left)

        section_header(sf, "MATERIAL PROPERTIES")
        mat_fields = [("E₁f [GPa]","Ef1",130.0),("E₂f [GPa]","Ef2",9.0),
                      ("ν₁₂f","nuf12",0.28),("G₁₂f [GPa]","Gf12",6.9),
                      ("Eₘ [GPa]","Em",4.2),("νₘ","num",0.34),
                      ("ρf [g/cm³]","rhof",1.78),("ρm [g/cm³]","rhom",1.27)]
        self.opt_entries = {}
        for lbl,key,defv in mat_fields:
            row = tk.Frame(sf, bg=BG_SECONDARY); row.pack(fill="x",padx=8,pady=2)
            tk.Label(row,text=lbl,font=("Helvetica",9),bg=BG_SECONDARY,
                     fg=TEXT_SEC,anchor="w",width=18).pack(side="left")
            e = make_entry(row,defv,width=9); e.pack(side="right"); self.opt_entries[key]=e

        section_header(sf, "OBJECTIVE")
        self.obj_var = tk.StringVar(value="max_E1")
        objectives = [("Maximize E₁","max_E1"),("Maximize E₂","max_E2"),
                      ("Maximize E₁/ρ  (Specific Stiffness)","max_spec"),
                      ("Minimize Density","min_rho")]
        for text,val in objectives:
            tk.Radiobutton(sf, text=text, variable=self.obj_var, value=val,
                           font=("Helvetica",9), bg=BG_SECONDARY, fg=TEXT_SEC,
                           selectcolor=ACCENT_DIM, activeforeground=ACCENT,
                           activebackground=BG_SECONDARY
                           ).pack(anchor="w", padx=10, pady=2)

        section_header(sf, "BOUNDS")
        bound_fields = [("Vf min","vf_min",0.30),("Vf max","vf_max",0.70)]
        self.bound_entries = {}
        for lbl,key,defv in bound_fields:
            row = tk.Frame(sf, bg=BG_SECONDARY); row.pack(fill="x",padx=8,pady=2)
            tk.Label(row,text=lbl,font=("Helvetica",9),bg=BG_SECONDARY,
                     fg=TEXT_SEC,anchor="w",width=10).pack(side="left")
            e = make_entry(row,defv,width=9); e.pack(side="right"); self.bound_entries[key]=e

        tk.Frame(sf, bg=BORDER, height=1).pack(fill="x", padx=8, pady=8)
        accent_btn(sf, "▶  OPTIMIZE (SciPy)", self._optimize, bg=self.color, fg="#000").pack(fill="x",padx=8,pady=4)
        accent_btn(sf, "📈  SWEEP & PLOT", self._sweep, bg=self.dim, fg=self.color
                   ).pack(fill="x",padx=8,pady=2)

        right = tk.Frame(pw, bg=BG_PRIMARY)
        pw.add(right, minsize=550)
        self.opt_result = tk.Frame(right, bg=BG_PRIMARY)
        self.opt_result.pack(fill="x", padx=16, pady=10)
        self.opt_fig = tk.Frame(right, bg=BG_PRIMARY)
        self.opt_fig.pack(fill="both", expand=True, padx=16)

    def _get_p(self):
        return {k:float(e.get()) for k,e in self.opt_entries.items()}

    def _objective(self, x, p, obj):
        Vf = float(x[0])
        res = halpin_tsai(p["Ef1"],p["Ef2"],p["nuf12"],p["Gf12"],p["Em"],p["num"],Vf)
        rho = p["rhof"]*Vf + p["rhom"]*(1-Vf)
        if obj == "max_E1": return -res["E1"]
        if obj == "max_E2": return -res["E2"]
        if obj == "max_spec": return -res["E1"]/rho
        if obj == "min_rho": return rho
        return -res["E1"]

    def _optimize(self):
        try: p = self._get_p()
        except: messagebox.showerror("Error","Invalid inputs"); return
        obj = self.obj_var.get()
        vf_min = float(self.bound_entries["vf_min"].get())
        vf_max = float(self.bound_entries["vf_max"].get())
        bounds = [(vf_min, vf_max)]
        result = minimize(lambda x: self._objective(x,p,obj), x0=[0.55],
                          bounds=bounds, method="L-BFGS-B")

        Vf_opt = float(result.x[0])
        res = halpin_tsai(p["Ef1"],p["Ef2"],p["nuf12"],p["Gf12"],p["Em"],p["num"],Vf_opt)
        rho = p["rhof"]*Vf_opt + p["rhom"]*(1-Vf_opt)

        for w in self.opt_result.winfo_children(): w.destroy()
        g = tk.Frame(self.opt_result, bg=BG_PRIMARY); g.pack(fill="x")
        items = [
            (f"Optimal Vf",f"{Vf_opt:.4f}"),
            (f"E₁",f"{res['E1']:.2f} GPa"),
            (f"E₂",f"{res['E2']:.2f} GPa"),
            (f"G₁₂",f"{res['G12']:.2f} GPa"),
            (f"ν₁₂",f"{res['nu12']:.4f}"),
            (f"Density",f"{rho:.3f} g/cm³"),
            (f"E₁/ρ",f"{res['E1']/rho:.2f} GPa·cm³/g"),
        ]
        for idx,(lbl,val) in enumerate(items):
            c2 = tk.Frame(g, bg=BG_CARD, padx=12, pady=8)
            c2.grid(row=0, column=idx, padx=4, pady=4, sticky="nsew")
            g.columnconfigure(idx, weight=1)
            tk.Label(c2, text=lbl, font=("Helvetica",8), bg=BG_CARD, fg=TEXT_MUTED).pack(anchor="w")
            tk.Label(c2, text=val, font=("Courier",12,"bold"), bg=BG_CARD, fg=self.color).pack(anchor="w")

    def _sweep(self):
        try: p = self._get_p()
        except: return
        vf_min = float(self.bound_entries["vf_min"].get())
        vf_max = float(self.bound_entries["vf_max"].get())
        Vf_arr = np.linspace(vf_min, vf_max, 80)

        E1_arr=[]; E2_arr=[]; G12_arr=[]; rho_arr=[]; spec_arr=[]
        for Vf in Vf_arr:
            r = halpin_tsai(p["Ef1"],p["Ef2"],p["nuf12"],p["Gf12"],p["Em"],p["num"],Vf)
            rho = p["rhof"]*Vf + p["rhom"]*(1-Vf)
            E1_arr.append(r["E1"]); E2_arr.append(r["E2"])
            G12_arr.append(r["G12"]); rho_arr.append(rho)
            spec_arr.append(r["E1"]/rho)

        for w in self.opt_fig.winfo_children(): w.destroy()
        fig = Figure(figsize=(10,3.5), dpi=100, facecolor=PLT_BG)
        plots = [
            (E1_arr,"E₁ (GPa)","Longitudinal Modulus",ACCENT),
            (spec_arr,"E₁/ρ  [GPa·cm³/g]","Specific Stiffness",SUCCESS),
            (rho_arr,"ρ  (g/cm³)","Density",WARNING),
        ]
        for i,(data,ylabel,title,color) in enumerate(plots):
            ax = fig.add_subplot(1,3,i+1, facecolor=PLT_AX_BG)
            style_ax(ax)
            ax.plot(Vf_arr, data, color=color, lw=2)
            ax.fill_between(Vf_arr, data, alpha=0.1, color=color)
            ax.set_xlabel("Fiber Volume Fraction Vf", fontsize=8)
            ax.set_ylabel(ylabel, fontsize=8)
            ax.set_title(title, fontsize=9, color=TEXT_PRIMARY)
        fig.tight_layout(pad=1.5)
        self.embed_figure(self.opt_fig, fig, toolbar=False)


# ══════════════════════════════════════════════════════════════════════════════
# 8. SENSITIVITY ANALYSIS
# ══════════════════════════════════════════════════════════════════════════════

class SensitivityPage(BasePage):
    def __init__(self, parent, app):
        super().__init__(parent, app, mid="sensitivity")
        self.page_header("Sensitivity Analysis",
                         "Tornado Charts · Spider Plots · Parameter Influence")
        self._build()

    def _build(self):
        pw = tk.PanedWindow(self, orient="horizontal", bg=BG_PRIMARY, sashwidth=4)
        pw.pack(fill="both", expand=True)

        left = tk.Frame(pw, bg=BG_SECONDARY, width=280)
        pw.add(left, minsize=240); left.pack_propagate(False)
        sf = scrollable_frame(left)

        section_header(sf, "BASE CASE")
        base_fields = [("E₁f [GPa]","Ef1",130.0),("E₂f [GPa]","Ef2",9.0),
                       ("ν₁₂f","nuf12",0.28),("G₁₂f [GPa]","Gf12",6.9),
                       ("Eₘ [GPa]","Em",4.2),("νₘ","num",0.34),("Vf","Vf",0.60)]
        self.base_entries = {}
        for lbl,key,defv in base_fields:
            row = tk.Frame(sf, bg=BG_SECONDARY); row.pack(fill="x",padx=8,pady=2)
            tk.Label(row,text=lbl,font=("Helvetica",9),bg=BG_SECONDARY,
                     fg=TEXT_SEC,anchor="w",width=16).pack(side="left")
            e = make_entry(row,defv,width=9); e.pack(side="right"); self.base_entries[key]=e

        section_header(sf, "PERTURBATION")
        row = tk.Frame(sf, bg=BG_SECONDARY); row.pack(fill="x",padx=8,pady=2)
        tk.Label(row, text="±% change", font=("Helvetica",9), bg=BG_SECONDARY,
                 fg=TEXT_SEC, anchor="w", width=12).pack(side="left")
        self.pct_entry = make_entry(row, 20.0, width=9); self.pct_entry.pack(side="right")

        tk.Frame(sf, bg=BORDER, height=1).pack(fill="x", padx=8, pady=8)
        accent_btn(sf, "▶  RUN SENSITIVITY", self._run, bg=self.color, fg="#000").pack(fill="x",padx=8,pady=4)

        right = tk.Frame(pw, bg=BG_PRIMARY)
        pw.add(right, minsize=550)
        self.sa_fig = right

    def _run(self):
        try:
            p0 = {k:float(e.get()) for k,e in self.base_entries.items()}
            pct = float(self.pct_entry.get())/100.0
        except: messagebox.showerror("Error","Invalid inputs"); return

        def get_E1(p):
            r = halpin_tsai(p["Ef1"],p["Ef2"],p["nuf12"],p["Gf12"],p["Em"],p["num"],p["Vf"])
            return r["E1"]

        E1_base = get_E1(p0)
        param_keys = ["Ef1","Ef2","nuf12","Gf12","Em","num","Vf"]
        param_names = ["E₁f","E₂f","ν₁₂f","G₁₂f","Eₘ","νₘ","Vf"]
        deltas = []
        for k in param_keys:
            p_up = dict(p0); p_up[k] = p0[k]*(1+pct)
            p_dn = dict(p0); p_dn[k] = p0[k]*(1-pct)
            delta = (get_E1(p_up) - get_E1(p_dn)) / 2.0
            deltas.append(delta)

        # Sort by absolute sensitivity
        order = sorted(range(len(deltas)), key=lambda i: abs(deltas[i]), reverse=True)

        for w in self.sa_fig.winfo_children(): w.destroy()
        fig = Figure(figsize=(9,4.5), dpi=100, facecolor=PLT_BG)
        ax1 = fig.add_subplot(121, facecolor=PLT_AX_BG)
        ax2 = fig.add_subplot(122, facecolor=PLT_AX_BG)
        style_ax(ax1); style_ax(ax2)

        # Tornado chart
        sorted_names = [param_names[i] for i in order]
        sorted_deltas = [deltas[i] for i in order]
        colors = [ACCENT if d>=0 else DANGER for d in sorted_deltas]
        y_pos = range(len(sorted_names))
        ax1.barh(list(y_pos), sorted_deltas, color=colors, alpha=0.85, edgecolor=PLT_GRID, lw=0.5)
        ax1.set_yticks(list(y_pos)); ax1.set_yticklabels(sorted_names, fontsize=9)
        ax1.set_xlabel("ΔE₁ (GPa) per ±" + f"{int(pct*100)}% change", fontsize=8)
        ax1.set_title("Tornado Chart — E₁ Sensitivity", fontsize=9, color=TEXT_PRIMARY)
        ax1.axvline(0, color=PLT_TEXT, lw=0.8)

        # Spider chart
        Vf_arr = np.linspace(0.3, 0.75, 60)
        for k,name,color in [("Ef1","E₁f",ACCENT),("Em","Eₘ",ACCENT2),("Vf","Vf",SUCCESS)]:
            E1_vals = []
            for v in Vf_arr:
                p_t = dict(p0); p_t["Vf"] = v
                E1_vals.append(get_E1(p_t))
            ax2.plot(Vf_arr, E1_vals, color=color, lw=2, label=name)
        ax2.set_xlabel("Fiber Volume Fraction Vf", fontsize=8)
        ax2.set_ylabel("E₁ (GPa)", fontsize=8)
        ax2.set_title("Parameter Sweep — E₁ vs Vf", fontsize=9, color=TEXT_PRIMARY)
        ax2.legend(fontsize=8)

        fig.tight_layout(pad=1.5)
        self.embed_figure(self.sa_fig, fig, toolbar=False)


# ══════════════════════════════════════════════════════════════════════════════
# 9. VISUALIZATION STUDIO
# ══════════════════════════════════════════════════════════════════════════════

class VisualizationPage(BasePage):
    def __init__(self, parent, app):
        super().__init__(parent, app, mid="visualization")
        self.page_header("Visualization Studio",
                         "Elastic Surfaces · Mohr's Circle · Property Maps")
        self._build()

    def _build(self):
        controls = tk.Frame(self, bg=BG_SECONDARY, padx=14, pady=10)
        controls.pack(fill="x")

        self.viz_var = tk.StringVar(value="elastic_surface")
        options = [("Directional Modulus","elastic_surface"),
                   ("Mohr's Circle (2D)","mohrs"),
                   ("Failure Envelope","fail_env"),
                   ("Stiffness Heatmap","stiff_heat")]
        tk.Label(controls, text="Visualization:", font=("Helvetica",10),
                 bg=BG_SECONDARY, fg=TEXT_SEC).pack(side="left", padx=(0,10))
        self.viz_btns = {}
        for text,val in options:
            b = tk.Button(controls, text=text, font=("Courier",9),
                          bg=BG_CARD, fg=TEXT_MUTED, relief="flat", cursor="hand2",
                          padx=10, pady=4,
                          command=lambda v=val: self._select_viz(v))
            b.pack(side="left", padx=3)
            self.viz_btns[val] = b

        accent_btn(controls, "▶  GENERATE", self._generate, bg=self.color, fg="#000", padx=14
                   ).pack(side="right", padx=10)

        tk.Frame(self, bg=BORDER, height=1).pack(fill="x")

        # Parameter inputs row
        inp_row = tk.Frame(self, bg=BG_PANEL, padx=14, pady=8)
        inp_row.pack(fill="x")
        params_viz = [("E₁ [GPa]","E1",130.0),("E₂ [GPa]","E2",9.0),
                      ("ν₁₂","nu12",0.28),("G₁₂ [GPa]","G12",6.9),
                      ("σ₁ [MPa]","s1",100.0),("σ₂ [MPa]","s2",30.0),
                      ("τ₁₂ [MPa]","s12",20.0)]
        self.viz_entries = {}
        for lbl,key,defv in params_viz:
            tk.Label(inp_row, text=lbl, font=("Helvetica",9), bg=BG_PANEL,
                     fg=TEXT_SEC).pack(side="left", padx=(6,2))
            e = make_entry(inp_row, defv, width=7); e.pack(side="left", padx=2)
            self.viz_entries[key] = e

        self.fig_frame = tk.Frame(self, bg=BG_PRIMARY)
        self.fig_frame.pack(fill="both", expand=True)
        self._select_viz("elastic_surface")

    def _select_viz(self, val):
        self.viz_var.set(val)
        for v,b in self.viz_btns.items():
            b.config(bg=ACCENT_DIM if v==val else BG_CARD,
                     fg=ACCENT if v==val else TEXT_MUTED)

    def _generate(self):
        viz = self.viz_var.get()
        for w in self.fig_frame.winfo_children(): w.destroy()
        try:
            p = {k:float(e.get()) for k,e in self.viz_entries.items()}
        except: return

        if viz == "elastic_surface":
            self._plot_elastic_surface(p)
        elif viz == "mohrs":
            self._plot_mohrs(p)
        elif viz == "fail_env":
            self._plot_fail_env(p)
        elif viz == "stiff_heat":
            self._plot_stiff_heat(p)

    def _plot_elastic_surface(self, p):
        E1,E2,nu12,G12 = p["E1"],p["E2"],p["nu12"],p["G12"]
        thetas = np.linspace(0, 2*np.pi, 360)
        Q = Q_matrix(E1,E2,nu12,G12)
        Ex_vals = []
        for th in thetas:
            Qb = Qbar_matrix(Q, math.degrees(th))
            S_bar = inv(Qb)
            Ex_vals.append(abs(1/S_bar[0,0]))
        Ex_vals = np.array(Ex_vals)

        fig = Figure(figsize=(8,4.5), dpi=100, facecolor=PLT_BG)
        ax1 = fig.add_subplot(121, projection="polar", facecolor=PLT_AX_BG)
        ax2 = fig.add_subplot(122, facecolor=PLT_AX_BG)
        style_ax(ax2)

        ax1.plot(thetas, Ex_vals, color=ACCENT, lw=2)
        ax1.fill(thetas, Ex_vals, alpha=0.15, color=ACCENT)
        ax1.set_title("E(θ) Polar Plot [GPa]", color=TEXT_PRIMARY, fontsize=9, pad=12)
        ax1.tick_params(colors=PLT_TEXT, labelsize=7)

        theta_deg = np.degrees(thetas)
        ax2.plot(theta_deg, Ex_vals, color=ACCENT, lw=2, label="E(θ)")
        ax2.axhline(E1, color=ACCENT2, lw=1, ls="--", label=f"E₁={E1} GPa")
        ax2.axhline(E2, color=SUCCESS, lw=1, ls="--", label=f"E₂={E2} GPa")
        ax2.set_xlabel("Angle θ (degrees)", fontsize=9)
        ax2.set_ylabel("E(θ) (GPa)", fontsize=9)
        ax2.set_title("Directional Modulus", fontsize=10, color=TEXT_PRIMARY)
        ax2.legend(fontsize=8)
        fig.tight_layout(pad=1.5)
        self.embed_figure(self.fig_frame, fig)

    def _plot_mohrs(self, p):
        s1,s2,s12 = p["s1"],p["s2"],p["s12"]
        sc = (s1+s2)/2.0; r = math.sqrt(((s1-s2)/2)**2 + s12**2)
        sp1 = sc+r; sp2 = sc-r

        theta_c = np.linspace(0, 2*np.pi, 200)
        x_c = sc + r*np.cos(theta_c)
        y_c = r*np.sin(theta_c)

        fig = Figure(figsize=(7,4.5), dpi=100, facecolor=PLT_BG)
        ax = fig.add_subplot(111, facecolor=PLT_AX_BG)
        style_ax(ax); ax.set_aspect("equal")

        ax.plot(x_c, y_c, color=ACCENT, lw=2, label="Mohr's Circle")
        ax.fill(x_c, y_c, alpha=0.08, color=ACCENT)
        ax.scatter([s1,s2], [s12,-s12], color=ACCENT2, s=60, zorder=5, label="Stress points")
        ax.scatter([sp1,sp2],[0,0], color=SUCCESS, s=80, zorder=5, marker="*",
                   label=f"σ_p = {sp1:.1f}, {sp2:.1f} MPa")
        ax.scatter([sc],[r], color=WARNING, s=60, zorder=5, label=f"τ_max = {r:.1f} MPa")
        ax.axhline(0, color=PLT_TEXT, lw=0.7); ax.axvline(sc, color=PLT_GRID, lw=0.5, ls=":")
        ax.set_xlabel("Normal Stress σ (MPa)", fontsize=9)
        ax.set_ylabel("Shear Stress τ (MPa)", fontsize=9)
        ax.set_title("Mohr's Circle", fontsize=11, color=TEXT_PRIMARY, pad=10)
        ax.legend(fontsize=8)
        fig.tight_layout(pad=1.5)
        self.embed_figure(self.fig_frame, fig)

    def _plot_fail_env(self, p):
        Xt,Xc,Yt,Yc,S12 = 1500,1200,50,200,90
        E1,E2,nu12,G12 = p["E1"],p["E2"],p["nu12"],p["G12"]
        s12_f = p["s12"]
        s1_arr = np.linspace(-Xc*1.1, Xt*1.1, 120)
        fig = Figure(figsize=(7,4.5), dpi=100, facecolor=PLT_BG)
        ax = fig.add_subplot(111, facecolor=PLT_AX_BG)
        style_ax(ax)
        for crit, color, name in [("hill",ACCENT,"Tsai-Hill"),("wu",ACCENT2,"Tsai-Wu")]:
            pos=[]; neg=[]
            for s1v in s1_arr:
                if crit=="hill":
                    X=Xt if s1v>=0 else Xc; Y=Yt
                    Ac=1/Y**2; Bc=-s1v/X**2; Cc=(s1v/X)**2+(s12_f/S12)**2-1
                else:
                    F1=1/Xt-1/Xc; F2=1/Yt-1/Yc
                    F11=1/(Xt*Xc); F22=1/(Yt*Yc); F12=-0.5*math.sqrt(F11*F22)
                    Ac=F22; Bc=F2+2*F12*s1v; Cc=F1*s1v+F11*s1v**2+(s12_f/S12)**2*1/(S12**2/(S12**2))-1
                    Cc=F1*s1v+F11*s1v**2-1
                disc=Bc**2-4*Ac*Cc
                if disc>=0:
                    pos.append((-Bc+math.sqrt(disc))/(2*Ac))
                    neg.append((-Bc-math.sqrt(disc))/(2*Ac))
                else: pos.append(np.nan); neg.append(np.nan)
            ax.plot(s1_arr, pos, color=color, lw=2, label=name)
            ax.plot(s1_arr, neg, color=color, lw=2)
            ax.fill_between(s1_arr, neg, pos, alpha=0.07, color=color)
        ax.scatter([p["s1"]],[p["s2"]], color=WARNING, s=80, zorder=5, label="Stress state")
        ax.axhline(0,color=PLT_GRID,lw=0.5); ax.axvline(0,color=PLT_GRID,lw=0.5)
        ax.set_xlabel("σ₁ (MPa)",fontsize=9); ax.set_ylabel("σ₂ (MPa)",fontsize=9)
        ax.set_title("Failure Envelope Comparison",fontsize=10,color=TEXT_PRIMARY)
        ax.legend(fontsize=8)
        fig.tight_layout(pad=1.5)
        self.embed_figure(self.fig_frame, fig)

    def _plot_stiff_heat(self, p):
        E1,E2,nu12,G12 = p["E1"],p["E2"],p["nu12"],p["G12"]
        Vf_arr = np.linspace(0.20,0.75,40)
        th_arr  = np.linspace(0,90,40)
        E1_grid = np.zeros((40,40))
        for i,Vf in enumerate(Vf_arr):
            r = halpin_tsai(E1*Vf+4.2*(1-Vf), E2, nu12, G12, 4.2, 0.34, Vf)
            for j,th in enumerate(th_arr):
                Q = Q_matrix(r["E1"],r["E2"],r["nu12"],r["G12"])
                Qb = Qbar_matrix(Q, th)
                S_b = inv(Qb)
                E1_grid[i,j] = abs(1/S_b[0,0])
        fig = Figure(figsize=(7,4.5), dpi=100, facecolor=PLT_BG)
        ax = fig.add_subplot(111, facecolor=PLT_AX_BG)
        im = ax.contourf(th_arr, Vf_arr, E1_grid, levels=20, cmap="plasma")
        fig.colorbar(im, ax=ax, label="E₁ (GPa)")
        ax.set_xlabel("Fiber Orientation θ (°)",fontsize=9)
        ax.set_ylabel("Fiber Volume Fraction Vf",fontsize=9)
        ax.set_title("Stiffness E₁ Heatmap  (θ vs Vf)",fontsize=10,color=TEXT_PRIMARY)
        ax.tick_params(colors=PLT_TEXT,labelsize=8)
        fig.tight_layout(pad=1.5)
        self.embed_figure(self.fig_frame, fig)


# ══════════════════════════════════════════════════════════════════════════════
# 10. EXPORT CENTER
# ══════════════════════════════════════════════════════════════════════════════

class ExportPage(BasePage):
    def __init__(self, parent, app):
        super().__init__(parent, app, mid="export")
        self.page_header("Export Center", "Excel · CSV · JSON Report")
        self._build()

    def _build(self):
        inner = tk.Frame(self, bg=BG_PRIMARY)
        inner.pack(fill="both", expand=True, padx=30, pady=20)

        tk.Label(inner, text="Export computed results to various formats.",
                 font=("Helvetica",11), bg=BG_PRIMARY, fg=TEXT_SEC
                 ).pack(anchor="w", pady=(0,16))

        btns = [
            ("📊  Export Stiffness Matrix to Excel",   self._export_excel),  # colorized below
            ("📄  Export Material Database to JSON",   self._export_json),
            ("📋  Export Summary CSV",                 self._export_csv),
        ]
        for text, cmd in btns:
            accent_btn(inner, text, cmd, padx=20).pack(anchor="w", pady=5)

        tk.Frame(inner, bg=BORDER, height=1).pack(fill="x", pady=16)
        self.log = scrolledtext.ScrolledText(inner, height=12, bg=BG_PANEL,
                                             fg=TEXT_SEC, font=("Courier",9),
                                             relief="flat", highlightthickness=0)
        self.log.pack(fill="both", expand=True)
        self.log.insert("end","Export log will appear here.\n")

    def _log(self, msg):
        self.log.insert("end", msg+"\n")
        self.log.see("end")

    def _export_excel(self):
        C = getattr(self.app, "last_C", None)
        if C is None:
            # Use a default
            C = stiffness_isotropic(200.0, 0.30)
            self._log("No matrix computed — using default isotropic steel.")
        path = filedialog.asksaveasfilename(defaultextension=".xlsx",
                                            filetypes=[("Excel","*.xlsx")])
        if not path: return
        wb = openpyxl.Workbook(); ws = wb.active; ws.title = "Stiffness Matrix"
        thin = Side(style="thin", color="2D3748")
        brd  = Border(left=thin, right=thin, top=thin, bottom=thin)
        ws.cell(1,1,"CMVL — Stiffness Tensor [C] (GPa)").font = Font(bold=True, color="00C8FF", size=12)
        ws.merge_cells("A1:H1")
        lbl = ["σ₁","σ₂","σ₃","σ₄","σ₅","σ₆"]
        for j,l in enumerate([""] + lbl):
            c = ws.cell(3, j+1, l); c.font = Font(bold=True, color="94A3B8"); c.border = brd
        for i,row in enumerate(C):
            ws.cell(4+i, 1, lbl[i]).font = Font(bold=True, color="94A3B8")
            for j,v in enumerate(row):
                c = ws.cell(4+i, j+2, round(float(v),6))
                c.border = brd
                c.alignment = Alignment(horizontal="center")
                is_nz = abs(v) > np.max(np.abs(C))*1e-9
                is_d  = (i==j)
                if is_d and is_nz:
                    c.font = Font(bold=True, color="00C8FF")
                    c.fill = PatternFill("solid", fgColor="1A2A10")
                elif is_nz:
                    c.font = Font(bold=True, color="F59E0B")
                    c.fill = PatternFill("solid", fgColor="0E1F2F")
        wb.save(path)
        self._log(f"✔ Excel exported: {path}")

    def _export_json(self):
        path = filedialog.asksaveasfilename(defaultextension=".json",
                                            filetypes=[("JSON","*.json")])
        if not path: return
        data = {"fibers": FIBER_DB, "matrices": MATRIX_DB}
        with open(path,"w") as f:
            json.dump(data, f, indent=2)
        self._log(f"✔ JSON exported: {path}")

    def _export_csv(self):
        path = filedialog.asksaveasfilename(defaultextension=".csv",
                                            filetypes=[("CSV","*.csv")])
        if not path: return
        lines = ["CMVL Export — Material Database","",
                 "FIBER,E1_GPa,E2_GPa,nu12,G12_GPa,density_gcm3,Xt_MPa,Xc_MPa"]
        for name,props in FIBER_DB.items():
            lines.append(f"{name},{props.get('E1','')},{props.get('E2','')},{props.get('nu12','')},{props.get('G12','')},{props.get('density','')},{props.get('Xt','')},{props.get('Xc','')}")
        lines += ["","MATRIX,E_GPa,nu,G_GPa,density_gcm3,Xt_MPa"]
        for name,props in MATRIX_DB.items():
            lines.append(f"{name},{props.get('E','')},{props.get('nu','')},{props.get('G','')},{props.get('density','')},{props.get('Xt','')}")
        with open(path,"w") as f:
            f.write("\n".join(lines))
        self._log(f"✔ CSV exported: {path}")


# ══════════════════════════════════════════════════════════════════════════════
# 11. EDUCATIONAL MODE
# ══════════════════════════════════════════════════════════════════════════════

THEORY_CONTENT = {
    "Classical Lamination Theory": """
CLASSICAL LAMINATION THEORY (CLT)
══════════════════════════════════════

The ABD matrix relates midplane strains/curvatures to in-plane forces/moments:

  ⌈ N ⌉   ⌈ A  B ⌉ ⌈ ε⁰ ⌉
  ⌊ M ⌋ = ⌊ B  D ⌋ ⌊  κ ⌋

A (extensional stiffness):  Aᵢⱼ = Σ Q̄ᵢⱼ · (zₖ - zₖ₋₁)
B (coupling stiffness):      Bᵢⱼ = ½ Σ Q̄ᵢⱼ · (zₖ² - zₖ₋₁²)
D (bending stiffness):       Dᵢⱼ = ⅓ Σ Q̄ᵢⱼ · (zₖ³ - zₖ₋₁³)

Q̄ = transformed reduced stiffness matrix (depends on ply angle θ)

For a symmetric laminate: B = 0 (no bending-extension coupling)

Engineering constants from A matrix:
  Ex  = 1 / (a₁₁ · h)
  Ey  = 1 / (a₂₂ · h)
  Gxy = 1 / (a₃₃ · h)
  νxy = -a₁₂ / a₁₁
where [a] = [A]⁻¹ and h = total laminate thickness
""",

    "Micromechanics Models": """
MICROMECHANICS MODELS
═══════════════════════════════════

RULE OF MIXTURES (Voigt / Reuss):
  E₁ = Ef·Vf + Em·Vm                  (parallel, exact)
  E₂ = 1/(Vf/Ef + Vm/Em)             (series, lower bound)
  ν₁₂ = νf·Vf + νm·Vm

HALPIN-TSAI (semi-empirical):
  P = Pm · (1 + ξ·η·Vf) / (1 - η·Vf)
  η = (Pf/Pm - 1) / (Pf/Pm + ξ)
  ξ = reinforcement geometry parameter
    (ξ=2 for E₂, ξ=1 for G₁₂)

MORI-TANAKA:
  Uses Eshelby inclusion tensor
  More accurate for higher Vf
  Accounts for fiber interactions

SELF-CONSISTENT:
  Fiber embedded in "effective medium"
  Iterative solution
  Good for polycrystalline materials
""",

    "Failure Criteria": """
FAILURE CRITERIA FOR COMPOSITES
═══════════════════════════════════════

TSAI-HILL:
  (σ₁/X)² - σ₁σ₂/X² + (σ₂/Y)² + (τ₁₂/S)² = 1
  Single interactive parameter, no sign distinction

TSAI-WU (tensor polynomial):
  F₁σ₁ + F₂σ₂ + F₁₁σ₁² + F₂₂σ₂² + F₆₆τ₁₂² + 2F₁₂σ₁σ₂ = 1
  F₁₂ = -½√(F₁₁·F₂₂)  (interaction term)
  Distinguishes tension/compression

MAXIMUM STRESS:
  σ₁/X < 1,  σ₂/Y < 1,  τ₁₂/S < 1
  No interaction — conservative

HASHIN (mode-based):
  Fiber tension:    (σ₁/Xₜ)² + (τ₁₂/S)² = 1
  Matrix tension:   (σ₂/Yₜ)² + (τ₁₂/S)² = 1
  Different failure modes tracked

Reserve Factor (RF) = 1/FI
Safety: RF > 1.0  →  FI < 1.0
""",

    "Tensor Notation": """
ELASTICITY TENSOR NOTATION
══════════════════════════════════

FOURTH-ORDER TENSOR (full form):
  σᵢⱼ = Cᵢⱼₖₗ · εₖₗ   (Einstein summation)

VOIGT NOTATION (6×6 matrix):
  Map: 11→1, 22→2, 33→3, 23→4, 13→5, 12→6
  
  ⌈ σ₁ ⌉   ⌈ C₁₁ C₁₂ C₁₃ C₁₄ C₁₅ C₁₆ ⌉ ⌈ ε₁ ⌉
  ⌊ σ₆ ⌋ = ⌊ C₁₆ C₂₆ C₃₆ C₄₆ C₅₆ C₆₆ ⌋ ⌊ ε₆ ⌋

SYMMETRY CLASSES (independent constants):
  Triclinic:             21
  Monoclinic:            13
  Orthotropic:            9
  Transversely Isotropic: 5
  Cubic:                  3
  Isotropic:              2

STABILITY CONDITIONS (Born):
  All eigenvalues of C must be positive
  Necessary: C₁₁>0, C₁₁·C₂₂>C₁₂², etc.
""",
}

class EducationalPage(BasePage):
    def __init__(self, parent, app):
        super().__init__(parent, app, mid="education")
        self.page_header("Educational Mode", "Theory · Equations · Derivations")
        self._build()

    def _build(self):
        pw = tk.PanedWindow(self, orient="horizontal", bg=BG_PRIMARY, sashwidth=4)
        pw.pack(fill="both", expand=True)

        left = tk.Frame(pw, bg=BG_SECONDARY, width=240)
        pw.add(left, minsize=200); left.pack_propagate(False)

        section_header(left, "TOPICS")
        self.topic_var = tk.StringVar()
        for topic in THEORY_CONTENT:
            b = tk.Button(left, text=topic, font=("Helvetica",10),
                          bg=BG_PANEL, fg=TEXT_SEC, relief="flat",
                          cursor="hand2", anchor="w", padx=10, pady=8,
                          command=lambda t=topic: self._show_topic(t),
                          wraplength=200)
            b.pack(fill="x", padx=8, pady=2)

        right = tk.Frame(pw, bg=BG_PRIMARY)
        pw.add(right, minsize=500)
        self.text_area = scrolledtext.ScrolledText(
            right, bg=BG_PANEL, fg=TEXT_PRIMARY, font=("Courier",11),
            relief="flat", highlightthickness=0, padx=20, pady=16,
            state="disabled", wrap="word")
        self.text_area.pack(fill="both", expand=True, padx=10, pady=10)
        self._show_topic("Classical Lamination Theory")

    def _show_topic(self, topic):
        content = THEORY_CONTENT.get(topic, "No content available.")
        self.text_area.config(state="normal")
        self.text_area.delete("1.0","end")
        self.text_area.insert("end", content)
        self.text_area.config(state="disabled")


# ══════════════════════════════════════════════════════════════════════════════
# 12. ENGINEERING CALCULATOR
# ══════════════════════════════════════════════════════════════════════════════

class CalculatorPage(BasePage):
    def __init__(self, parent, app):
        super().__init__(parent, app, mid="calculator")
        self.page_header("Engineering Calculator Suite",
                         "Specific Modulus · Safety Factor · Pressure Vessel · Unit Converter")
        self._build()

    def _build(self):
        pw = tk.PanedWindow(self, orient="horizontal", bg=BG_PRIMARY, sashwidth=4)
        pw.pack(fill="both", expand=True)

        left = tk.Frame(pw, bg=BG_SECONDARY, width=260)
        pw.add(left, minsize=220); left.pack_propagate(False)

        section_header(left, "CALCULATOR")
        calcs = [
            ("Specific Modulus","spec_mod"),
            ("Specific Strength","spec_str"),
            ("Safety Factor","safety"),
            ("Pressure Vessel (thin-wall)","pv"),
            ("Buckling (plate)","buckle"),
            ("Unit Converter","units"),
        ]
        self.calc_var = tk.StringVar(value="spec_mod")
        for label,val in calcs:
            b = tk.Button(left, text=label, font=("Helvetica",10),
                          bg=BG_PANEL, fg=TEXT_SEC, relief="flat",
                          cursor="hand2", anchor="w", padx=12, pady=8,
                          command=lambda v=val: self._select_calc(v))
            b.pack(fill="x", padx=8, pady=2)

        right = tk.Frame(pw, bg=BG_PRIMARY)
        pw.add(right, minsize=500)
        self.calc_area = right
        self._select_calc("spec_mod")

    def _select_calc(self, calc_id):
        self.calc_var.set(calc_id)
        for w in self.calc_area.winfo_children(): w.destroy()
        builders = {
            "spec_mod":  self._build_spec_mod,
            "spec_str":  self._build_spec_str,
            "safety":    self._build_safety,
            "pv":        self._build_pv,
            "buckle":    self._build_buckle,
            "units":     self._build_units,
        }
        builders.get(calc_id, lambda: None)()

    def _calc_frame(self, title):
        tk.Label(self.calc_area, text=title, font=("Helvetica",13,"bold"),
                 bg=BG_PRIMARY, fg=ACCENT).pack(anchor="w", padx=20, pady=(14,4))
        tk.Frame(self.calc_area, bg=BORDER, height=1).pack(fill="x", padx=20, pady=4)
        f = tk.Frame(self.calc_area, bg=BG_PRIMARY)
        f.pack(fill="x", padx=20, pady=10)
        return f

    def _row(self, parent, label, default, key, result_var=None):
        row = tk.Frame(parent, bg=BG_PRIMARY); row.pack(fill="x", pady=3)
        tk.Label(row, text=label, font=("Helvetica",10), bg=BG_PRIMARY,
                 fg=TEXT_SEC, anchor="w", width=28).pack(side="left")
        e = make_entry(row, default, width=12); e.pack(side="left", padx=6)
        return e

    def _build_spec_mod(self):
        f = self._calc_frame("Specific Modulus  E/ρ")
        entries = {}
        entries["E"] = self._row(f,"Young's Modulus E (GPa)", 130.0,"E")
        entries["rho"] = self._row(f,"Density ρ (g/cm³)", 1.58,"rho")
        res_lbl = tk.Label(f, text="", font=("Courier",15,"bold"), bg=BG_PRIMARY, fg=SUCCESS)
        res_lbl.pack(anchor="w", pady=8)
        def compute():
            try:
                E = float(entries["E"].get()); rho = float(entries["rho"].get())
                spec = E/rho
                res_lbl.config(text=f"E/ρ  =  {spec:.2f}  GPa·cm³/g  =  {spec*1e6:.2f}  m²/s²")
            except: res_lbl.config(text="Invalid input", fg=DANGER)
        accent_btn(f,"Compute", compute).pack(anchor="w", pady=4)

    def _build_spec_str(self):
        f = self._calc_frame("Specific Strength  σ/ρ")
        entries = {}
        entries["s"] = self._row(f,"Tensile Strength σ (MPa)", 1500.0,"s")
        entries["rho"] = self._row(f,"Density ρ (g/cm³)", 1.58,"rho")
        res_lbl = tk.Label(f, text="", font=("Courier",15,"bold"), bg=BG_PRIMARY, fg=SUCCESS)
        res_lbl.pack(anchor="w", pady=8)
        def compute():
            try:
                s = float(entries["s"].get()); rho = float(entries["rho"].get())
                spec = s/rho
                res_lbl.config(text=f"σ/ρ  =  {spec:.2f}  MPa·cm³/g  =  {spec*1e3:.0f}  m²/s²")
            except: res_lbl.config(text="Invalid input", fg=DANGER)
        accent_btn(f,"Compute", compute).pack(anchor="w", pady=4)

    def _build_safety(self):
        f = self._calc_frame("Safety Factor")
        entries = {}
        entries["fail"] = self._row(f,"Failure Load / Strength", 1500.0,"fail")
        entries["appl"] = self._row(f,"Applied Load / Stress", 600.0,"appl")
        res_lbl = tk.Label(f, text="", font=("Courier",15,"bold"), bg=BG_PRIMARY, fg=SUCCESS)
        res_lbl.pack(anchor="w", pady=8)
        def compute():
            try:
                fail = float(entries["fail"].get()); appl = float(entries["appl"].get())
                sf = fail/appl
                clr = SUCCESS if sf>=2.0 else (WARNING if sf>=1.5 else DANGER)
                status = "✔ SAFE" if sf>=1.5 else "⚠ MARGINAL" if sf>=1.0 else "✗ FAILED"
                res_lbl.config(text=f"SF  =  {sf:.3f}   {status}", fg=clr)
            except: res_lbl.config(text="Invalid input", fg=DANGER)
        accent_btn(f,"Compute", compute).pack(anchor="w", pady=4)

    def _build_pv(self):
        f = self._calc_frame("Thin-Wall Pressure Vessel")
        entries = {}
        entries["p"] = self._row(f,"Internal Pressure p (MPa)", 1.5,"p")
        entries["r"] = self._row(f,"Inner Radius r (mm)", 150.0,"r")
        entries["t"] = self._row(f,"Wall Thickness t (mm)", 3.0,"t")
        res_lbl = tk.Label(f, text="", font=("Courier",13,"bold"), bg=BG_PRIMARY, fg=SUCCESS,
                           justify="left")
        res_lbl.pack(anchor="w", pady=8)
        def compute():
            try:
                p = float(entries["p"].get())
                r = float(entries["r"].get())
                t = float(entries["t"].get())
                s_hoop = p*r/t
                s_long = p*r/(2*t)
                res_lbl.config(text=(
                    f"σ_hoop (circumferential)  =  {s_hoop:.2f} MPa\n"
                    f"σ_long (longitudinal)     =  {s_long:.2f} MPa\n"
                    f"R/t ratio                 =  {r/t:.2f}"
                ), fg=SUCCESS)
            except: res_lbl.config(text="Invalid input", fg=DANGER)
        accent_btn(f,"Compute", compute).pack(anchor="w", pady=4)

    def _build_buckle(self):
        f = self._calc_frame("Plate Buckling (Simply Supported)")
        entries = {}
        entries["D"] = self._row(f,"Bending Stiffness D (N·m)", 100.0,"D")
        entries["a"] = self._row(f,"Plate Length a (m)", 0.5,"a")
        entries["b"] = self._row(f,"Plate Width b (m)", 0.25,"b")
        entries["t"] = self._row(f,"Thickness h (mm)", 3.0,"t")
        res_lbl = tk.Label(f, text="", font=("Courier",13,"bold"), bg=BG_PRIMARY, fg=SUCCESS)
        res_lbl.pack(anchor="w", pady=8)
        def compute():
            try:
                D = float(entries["D"].get())
                a = float(entries["a"].get())
                b = float(entries["b"].get())
                # Ncr = k*pi²*D/b²  where k=4 for uniaxial
                Ncr = 4 * math.pi**2 * D / b**2
                res_lbl.config(text=f"Ncr  =  {Ncr:.2f}  N/m  (uniaxial, k=4)", fg=SUCCESS)
            except: res_lbl.config(text="Invalid input", fg=DANGER)
        accent_btn(f,"Compute", compute).pack(anchor="w", pady=4)

    def _build_units(self):
        f = self._calc_frame("Unit Converter")
        entries = {}
        entries["val"] = self._row(f,"Value", 1.0,"val")
        conv_options = ["GPa → MPa","MPa → GPa","GPa → psi","psi → MPa",
                        "mm → inch","inch → mm","kg/m³ → g/cm³","g/cm³ → kg/m³"]
        self.conv_var = tk.StringVar(value=conv_options[0])
        row = tk.Frame(f, bg=BG_PRIMARY); row.pack(fill="x", pady=3)
        tk.Label(row, text="Conversion", font=("Helvetica",10), bg=BG_PRIMARY,
                 fg=TEXT_SEC, anchor="w", width=28).pack(side="left")
        om = ttk.OptionMenu(row, self.conv_var, conv_options[0], *conv_options)
        om.pack(side="left")
        res_lbl = tk.Label(f, text="", font=("Courier",15,"bold"), bg=BG_PRIMARY, fg=SUCCESS)
        res_lbl.pack(anchor="w", pady=8)
        factors = {"GPa → MPa":1e3,"MPa → GPa":1e-3,"GPa → psi":145037.7,
                   "psi → MPa":6.895e-3,"mm → inch":1/25.4,"inch → mm":25.4,
                   "kg/m³ → g/cm³":1e-3,"g/cm³ → kg/m³":1e3}
        def compute():
            try:
                val = float(entries["val"].get())
                conv = self.conv_var.get()
                result = val * factors[conv]
                res_lbl.config(text=f"{val}  {conv.split('→')[0].strip()}  →  {result:.6g}  {conv.split('→')[1].strip()}")
            except: res_lbl.config(text="Invalid input", fg=DANGER)
        accent_btn(f,"Convert", compute).pack(anchor="w", pady=4)


# ══════════════════════════════════════════════════════════════════════════════
# MAIN APPLICATION WINDOW
# ══════════════════════════════════════════════════════════════════════════════

class CMVLApp:
    def __init__(self, root):
        self.root = root
        self.root.title("CMVL — Composite Materials Virtual Laboratory  v1.0")
        self.root.configure(bg=BG_PRIMARY)
        self.root.minsize(1200, 750)
        self._maximise()

        self.last_C = None  # shared state for export

        self._build_header()
        self._build_body()
        self._select_module("database")

    def _maximise(self):
        """Maximise the window cross-platform (Windows, macOS, Linux)."""
        import sys
        try:
            # Windows & most Linux window managers
            self.root.state("zoomed")
        except tk.TclError:
            try:
                # macOS (Aqua) and some Linux WMs
                self.root.attributes("-zoomed", True)
            except tk.TclError:
                # Fallback: fill the screen manually
                w = self.root.winfo_screenwidth()
                h = self.root.winfo_screenheight()
                self.root.geometry(f"{w}x{h}+0+0")

    # ── HEADER ─────────────────────────────────────────────────────────────

    def _build_header(self):
        hdr = tk.Frame(self.root, bg=BG_SECONDARY, height=60)
        hdr.pack(fill="x", side="top"); hdr.pack_propagate(False)

        left = tk.Frame(hdr, bg=BG_SECONDARY)
        left.pack(side="left", padx=18, fill="y")

        # Logo block — rainbow gradient approximated with colored segments
        logo_outer = tk.Frame(left, bg=BG_SECONDARY)
        logo_outer.pack(side="left", anchor="center", pady=12)
        for color in ["#00C8FF","#a78bfa","#34d399","#f59e0b"]:
            tk.Frame(logo_outer, bg=color, width=6, height=32).pack(side="left")
        tk.Label(logo_outer, text=" C", font=("Helvetica",17,"bold"),
                 bg=BG_SECONDARY, fg=TEXT_PRIMARY).pack(side="left")

        tk.Label(left, text="omposite Materials", font=("Helvetica",15,"bold"),
                 bg=BG_SECONDARY, fg=TEXT_PRIMARY).pack(side="left", anchor="center")
        tk.Label(left, text=" Virtual Lab",
                 font=("Helvetica",15,"bold"), bg=BG_SECONDARY, fg=ACCENT
                 ).pack(side="left", anchor="center")

        # Version badge
        badge = tk.Frame(left, bg="#003d52", padx=8, pady=2)
        badge.pack(side="left", padx=(14,0), anchor="center")
        tk.Label(badge, text="v1.0", font=("Courier",8,"bold"),
                 bg="#003d52", fg=ACCENT).pack()

        right = tk.Frame(hdr, bg=BG_SECONDARY)
        right.pack(side="right", padx=18, fill="y")
        tk.Label(right,
                 text="CLT  ·  Micromechanics  ·  Failure  ·  Optimization  ·  Python",
                 font=("Courier",8), bg=BG_SECONDARY, fg=TEXT_MUTED
                 ).pack(anchor="center", pady=(8,0))
        # Colored module dots row
        dot_row = tk.Frame(right, bg=BG_SECONDARY)
        dot_row.pack(anchor="e", pady=(2,8))
        for color in [mc[0] for mc in MODULE_COLORS.values()]:
            tk.Frame(dot_row, bg=color, width=8, height=8).pack(side="left", padx=1)

        # Bottom rainbow accent line
        rainbow = tk.Frame(self.root, height=3, bg=BG_SECONDARY)
        rainbow.pack(fill="x")
        colors_bar = ["#00C8FF","#a78bfa","#34d399","#f59e0b","#f87171",
                      "#fb923c","#4ade80","#e879f9","#38bdf8","#94a3b8","#fcd34d","#67e8f9"]
        for c in colors_bar:
            tk.Frame(rainbow, bg=c, height=3).pack(side="left", fill="both", expand=True)

    # ── BODY ───────────────────────────────────────────────────────────────

    def _build_body(self):
        paned = tk.PanedWindow(self.root, orient="horizontal",
                               bg=BG_PRIMARY, sashwidth=4)
        paned.pack(fill="both", expand=True)

        # SIDEBAR NAV
        nav = tk.Frame(paned, bg=BG_SECONDARY, width=220)
        paned.add(nav, minsize=180); nav.pack_propagate(False)
        self._build_nav(nav)

        # MAIN CONTENT
        self.content = tk.Frame(paned, bg=BG_PRIMARY)
        paned.add(self.content, minsize=900)

        # status bar
        sb = tk.Frame(self.root, bg="#0d1420", height=26)
        sb.pack(fill="x", side="bottom"); sb.pack_propagate(False)
        # Left colored dot + status text
        self.status_dot = tk.Frame(sb, bg=ACCENT, width=8, height=8)
        self.status_dot.pack(side="left", padx=(12,6), pady=9)
        self.status_var = tk.StringVar(value="Ready  ·  CMVL v1.0  ·  University of Colorado Boulder")
        tk.Label(sb, textvariable=self.status_var, font=("Courier",8),
                 bg="#0d1420", fg=TEXT_SEC, anchor="w"
                 ).pack(side="left", fill="y")
        tk.Label(sb, text="NumPy · SciPy · Matplotlib · Python",
                 font=("Courier",8), bg="#0d1420", fg=TEXT_MUTED, padx=12
                 ).pack(side="right", fill="y")

    def _build_nav(self, nav):
        canvas = tk.Canvas(nav, bg=BG_SECONDARY, highlightthickness=0)
        vsb = ttk.Scrollbar(nav, orient="vertical", command=canvas.yview)
        nav_inner = tk.Frame(canvas, bg=BG_SECONDARY)
        nav_inner.bind("<Configure>", lambda e: canvas.configure(
            scrollregion=canvas.bbox("all")))
        canvas.create_window((0,0), window=nav_inner, anchor="nw")
        canvas.configure(yscrollcommand=vsb.set)
        vsb.pack(side="right", fill="y")
        canvas.pack(side="left", fill="both", expand=True)

        modules = [
            ("🗄  Material Database",      "database"),
            ("⚙  Constitutive Matrix",     "constitutive"),
            ("🔬  Micromechanics",          "micro"),
            ("📐  Laminate Analysis (CLT)", "laminate"),
            ("💥  Failure Analysis",        "failure"),
            ("📊  Stress-Strain Simulator", "stress_strain"),
            ("🎯  Optimization",            "optimization"),
            ("📈  Sensitivity Analysis",    "sensitivity"),
            ("🎨  Visualization Studio",    "visualization"),
            ("💾  Export Center",           "export"),
            ("📚  Educational Mode",        "education"),
            ("🔧  Engineering Calculator",  "calculator"),
        ]
        self.nav_btns = {}
        self.nav_indicators = {}
        self.pages = {}

        section_label = tk.Frame(nav_inner, bg=BG_SECONDARY)
        section_label.pack(fill="x", padx=10, pady=(12,4))
        tk.Label(section_label, text="MODULES", font=("Courier",8),
                 bg=BG_SECONDARY, fg=TEXT_MUTED).pack(side="left")
        tk.Frame(section_label, bg=BORDER, height=1
                 ).pack(side="left", fill="x", expand=True, padx=(6,0))

        for label, mid in modules:
            color = module_accent(mid)
            dim   = module_dim(mid)

            row_f = tk.Frame(nav_inner, bg=BG_SECONDARY)
            row_f.pack(fill="x", padx=6, pady=1)

            # Left color indicator bar (hidden when inactive)
            ind = tk.Frame(row_f, bg=BG_SECONDARY, width=3)
            ind.pack(side="left", fill="y")
            self.nav_indicators[mid] = ind

            btn = tk.Button(row_f, text=label, font=("Helvetica",10),
                            bg=BG_PANEL, fg=TEXT_SEC, relief="flat",
                            cursor="hand2", anchor="w", padx=10, pady=8,
                            command=lambda m=mid: self._select_module(m))
            btn.pack(side="left", fill="x", expand=True)
            self.nav_btns[mid] = (btn, color, dim, row_f)

    def _select_module(self, mid):
        # Highlight active button with its module color
        for k, (btn, color, dim, row_f) in self.nav_btns.items():
            if k == mid:
                btn.config(bg=dim, fg=color)
                row_f.config(bg=dim)
                self.nav_indicators[k].config(bg=color)
            else:
                btn.config(bg=BG_PANEL, fg=TEXT_SEC)
                row_f.config(bg=BG_SECONDARY)
                self.nav_indicators[k].config(bg=BG_SECONDARY)

        # Clear current content
        for w in self.content.winfo_children(): w.destroy()

        # Lazy-build pages
        page_classes = {
            "database":    MaterialDBPage,
            "constitutive":ConstitutivePage,
            "micro":       MicromechanicsPage,
            "laminate":    LaminatePage,
            "failure":     FailurePage,
            "stress_strain":StressStrainPage,
            "optimization":OptimizationPage,
            "sensitivity": SensitivityPage,
            "visualization":VisualizationPage,
            "export":      ExportPage,
            "education":   EducationalPage,
            "calculator":  CalculatorPage,
        }
        cls = page_classes.get(mid)
        if cls:
            page = cls(self.content, self)
            page.pack(fill="both", expand=True)

        self.status_var.set(f"Module: {mid.replace('_',' ').title()}  ·  CMVL v1.0  ·  University of Colorado Boulder")
        self.status_dot.config(bg=module_accent(mid))


# ─────────────────────────────────────────────────────────────────────────────
# ENTRY POINT
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    root = tk.Tk()
    try:
        root.tk.call("wm", "iconphoto", root._w)
    except Exception:
        pass
    app = CMVLApp(root)
    # Re-assert maximised state after the event loop initialises
    # (needed on some Linux desktops that override geometry on first draw)
    root.after(50, app._maximise)
    root.mainloop()
