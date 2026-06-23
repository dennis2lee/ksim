"""
Generate all manuscript figures from simulation data.
Each figure uses the same seed/parameters as reproduce_manuscript_numbers.py.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import matplotlib.patches as mpatches
import os

OUT = os.path.dirname(__file__)
# Colorblind-safe palette (Wong 2011)
C_BLUE = '#0072B2'
C_ORANGE = '#E69F00'
C_GREEN = '#009E73'
C_RED = '#D55E00'
C_PURPLE = '#CC79A7'
C_GRAY = '#999999'

plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Times New Roman', 'DejaVu Serif'],
    'font.size': 9,
    'axes.labelsize': 10,
    'axes.titlesize': 11,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'legend.fontsize': 8,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'axes.spines.top': False,
    'axes.spines.right': False,
})

def save(fig, name):
    for ext in ['png', 'pdf', 'svg']:
        fig.savefig(os.path.join(OUT, f'{name}.{ext}'))
    plt.close(fig)
    print(f'  Saved {name}.png/pdf/svg')

# =====================================================================
# FIGURE 1: Protocol flowchart — clean top-to-bottom layout
# =====================================================================
print("Figure 1: Protocol flowchart")
fig, ax = plt.subplots(figsize=(7.5, 9.0))
ax.set_xlim(0, 15)
ax.set_ylim(0, 13)
ax.axis('off')

BW = 3.8   # box width
BH = 1.0   # box height
BFS = 8    # box font size
VGAP = 0.4 # vertical gap between box edge and next arrow start

def box(x, y, w, h, text, color='white', ec='black', lw=1.2, fs=BFS):
    rect = FancyBboxPatch((x-w/2, y-h/2), w, h,
        boxstyle="round,pad=0.15", facecolor=color, edgecolor=ec, linewidth=lw)
    ax.add_patch(rect)
    ax.text(x, y, text, ha='center', va='center', fontsize=fs, fontfamily='serif',
            linespacing=1.35)

def arrow(x1, y1, x2, y2, label='', color='black', ha='left', dx=0.2):
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
        arrowprops=dict(arrowstyle='->', color=color, lw=1.2))
    if label:
        mx, my = (x1+x2)/2 + dx, (y1+y2)/2
        ax.text(mx, my, label, fontsize=7, color=color, va='center', ha=ha)

cx = 4.5   # center x for main column
rx = 11.5  # center x for outcome column
OW = 3.4   # outcome box width
OH = 0.85  # outcome box height

# Derived edges
main_R = cx + BW/2    # right edge of main boxes
main_B = lambda y: y - BH/2   # bottom edge
main_T = lambda y: y + BH/2   # top edge
out_L = rx - OW/2     # left edge of outcome boxes

# Step 0
y0 = 12.0
box(cx, y0, BW, BH, 'Step 0: Eligibility\nCKD 3b–4, stable eGFR\nIS assay available', '#E8F5E9', C_GREEN)

# Step 1
y1 = y0 - BH - 2*VGAP
arrow(cx, main_B(y0), cx, main_T(y1))
box(cx, y1, BW, BH, 'Step 1: Measurement standardization\nFasting AM, duplicate draws\n→ target CV ≈ 0.15', '#E3F2FD', C_BLUE)

# Step 2
y2 = y1 - BH - 2*VGAP
arrow(cx, main_B(y1), cx, main_T(y2))
box(cx, y2, BW, BH, 'Step 2: Stage 1 (wk 0–24)\n2 × 3 AB crossover\n12 visits, 24 duplicate draws', '#FFF3E0', C_ORANGE)

# Step 3: Decision
y3 = y2 - BH - 2*VGAP
arrow(cx, main_B(y2), cx, main_T(y3))
box(cx, y3, BW, BH, 'Step 3: Classification\nobs_red = (mean_A − mean_B) / mean_A\nvs MDE threshold (14%)', '#FCE4EC', C_RED)

# Branch right: RESPONDER — arrow from main right edge to outcome left edge
yr = y3 + 0.15
arrow(main_R, yr, out_L, yr, '', C_GREEN)
ax.text((main_R+out_L)/2, yr+0.25, 'obs > 14%', fontsize=7, color=C_GREEN, ha='center')
box(rx, yr, OW, OH, 'RESPONDER\nContinue regimen', '#E8F5E9', C_GREEN, fs=8)

# Branch right-down: NON-RESPONDER
ynr = y3 - 1.8
arrow(main_R, y3-0.15, out_L, ynr, '', C_RED)
ax.text(8.5, y3-1.0, 'obs < 0%', fontsize=7, color=C_RED, ha='center')
box(rx, ynr, OW, OH, 'NON-RESPONDER\nStop, deprescribe', '#FFEBEE', C_RED, fs=8)

# Branch down: BORDERLINE
y4 = y3 - BH - 2*VGAP
arrow(cx, main_B(y3), cx, main_T(y4))
ax.text(cx+0.25, (main_B(y3)+main_T(y4))/2, '0–14%', fontsize=7, color=C_ORANGE, ha='left')

box(cx, y4, BW, BH, 'Step 4: Stage 2 (wk 24–36)\n+1 AB cycle for borderline only\nCombine 9+9 measures, MDE = 12%', '#FFF3E0', C_ORANGE)

# Stage 2 branches to same outcome boxes
arrow(main_R, y4+0.15, out_L, yr-0.3, '', C_GREEN)
ax.text(8.8, (y4+yr)/2+0.1, '> 12%', fontsize=7, color=C_GREEN, ha='center')
arrow(main_R, y4-0.15, out_L, ynr+0.3, '', C_RED)
ax.text(8.8, (y4+ynr)/2-0.1, '≤ 12%', fontsize=7, color=C_RED, ha='center')

# Step 5: Follow-up
y5 = y4 - BH - 2*VGAP
arrow(cx, main_B(y4), cx, main_T(y5))
box(cx, y5, BW, BH, 'Step 5: Follow-up\nResponders: IS every 3 months\nNon-resp: deprescribe review', '#F3E5F5', C_PURPLE)

# Phase labels on left margin
for label, y in [('Screening', y0), ('Preparation', y1), ('Trial', y2),
                 ('Decision', y3), ('Adaptive', y4), ('Monitoring', y5)]:
    ax.text(0.3, y, label, fontsize=7.5, color=C_GRAY, style='italic', va='center')

ax.set_title('Figure 1. N-of-1 protocol overview and decision rules',
             fontsize=10, fontweight='bold', loc='left', pad=12)
save(fig, 'figure1_protocol')

# =====================================================================
# FIGURE 2: MDE by design + operating characteristics
# =====================================================================
print("Figure 2: MDE and operating characteristics")
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6.5, 3.0))

# Panel A: MDE by CV and design
cvs_sweep = np.arange(0.10, 0.31, 0.01)
for nc, km, ls, label in [(2,3,'-','2×3 (24 wk)'), (3,3,'--','3×3 (36 wk)'), (4,3,':','4×3 (48 wk)')]:
    n = nc * km
    mdes = [1.645 * cv * np.sqrt(2/n) * 100 for cv in cvs_sweep]
    ax1.plot(cvs_sweep*100, mdes, ls, color=C_BLUE, linewidth=1.5, label=label)

ax1.axhline(y=14, color=C_GREEN, linewidth=0.8, linestyle='-', alpha=0.5)
ax1.text(28, 14.5, 'MDE = 14%', fontsize=7, color=C_GREEN)
ax1.axvline(x=15, color=C_ORANGE, linewidth=0.8, linestyle='--', alpha=0.5)
ax1.axvline(x=25, color=C_RED, linewidth=0.8, linestyle='--', alpha=0.5)
ax1.text(15.5, 30, 'CV = 0.15\n(target)', fontsize=6, color=C_ORANGE)
ax1.text(25.5, 30, 'CV = 0.25\n(native)', fontsize=6, color=C_RED)
ax1.set_xlabel('Within-person CV (%)')
ax1.set_ylabel('Minimum detectable effect (%)')
ax1.legend(frameon=False, loc='upper left')
ax1.set_title('A', fontweight='bold', loc='left')

# Panel B: Subgroup detection rates (single-run from reproduce script)
subgroups = ['Weak\n(10–20%)', 'Moderate\n(20–35%)', 'Strong\n(≥35%)']
rates = [76, 98, 100]  # from reproduce_manuscript_numbers.py
colors = [C_ORANGE, C_BLUE, C_GREEN]
bars = ax2.bar(subgroups, rates, color=colors, width=0.5, edgecolor='black', linewidth=0.5)
for bar, rate in zip(bars, rates):
    ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1.5,
             f'{rate}%', ha='center', fontsize=8, fontweight='bold')
ax2.set_ylabel('Single-run detection rate (%)')
ax2.set_ylim(0, 115)
ax2.axhline(y=80, color=C_GRAY, linewidth=0.8, linestyle='--', alpha=0.6)
ax2.text(2.4, 81, '80%', fontsize=7, color=C_GRAY)
ax2.set_title('B', fontweight='bold', loc='left')

fig.suptitle('Figure 2. MDE across designs (A) and subgroup detection rates (B)',
             fontsize=9, fontweight='bold', y=1.02)
plt.tight_layout()
save(fig, 'figure2_performance')

# =====================================================================
# FIGURE 3: EVSI / classification-efficiency frontier
# =====================================================================
print("Figure 3: Classification-efficiency frontier")

# Data from evsi_analysis.py (seed=800)
cv_pts = [0.30, 0.28, 0.25, 0.22, 0.20, 0.18, 0.15, 0.13, 0.12, 0.10]
ncc_pts = [801, 815, 853, 886, 891, 897, 924, 926, 929, 928]
sens_pts = [75, 77, 84, 88, 90, 91, 95, 97, 98, 99]
spec_pts = [93, 91, 89, 90, 88, 87, 87, 82, 82, 78]

fig, ax1 = plt.subplots(figsize=(5.5, 4.0))
ax2 = ax1.twinx()

cv_pct = [c*100 for c in cv_pts]

ln1 = ax1.plot(cv_pct, ncc_pts, 'o-', color=C_BLUE, linewidth=1.8,
               markersize=6, label='NCC (TP + TN)', zorder=3)
ln2 = ax2.plot(cv_pct, sens_pts, 's--', color=C_GREEN, linewidth=1.2,
               markersize=4, label='Sensitivity (%)', zorder=2)
ln3 = ax2.plot(cv_pct, spec_pts, '^--', color=C_RED, linewidth=1.2,
               markersize=4, label='Specificity (%)', zorder=2)

# Optimal zone shading (behind data)
ax1.axvspan(12, 15, alpha=0.10, color=C_ORANGE, zorder=0)

# Label in the middle of the chart where lines are spread apart
ax1.text(13.5, 868, 'Optimal\nCV zone', fontsize=8, ha='center',
    color=C_ORANGE, fontweight='bold', va='center', alpha=0.9)

peak_idx = ncc_pts.index(max(ncc_pts))
ax1.annotate(f'Peak = {max(ncc_pts)}',
    xy=(cv_pts[peak_idx]*100, max(ncc_pts)),
    xytext=(17, 935), fontsize=7, ha='left',
    arrowprops=dict(arrowstyle='->', color=C_BLUE, lw=0.8))

ax1.annotate('$180 per patient\n(12 extra draws)',
    xy=(15, 924), xytext=(23, 815), fontsize=7, color=C_GRAY, ha='center',
    arrowprops=dict(arrowstyle='->', color=C_GRAY, lw=0.8))

ax1.set_xlabel('Within-person CV (%)')
ax1.set_ylabel('Net correct classifications (NCC)', color=C_BLUE)
ax2.set_ylabel('Sensitivity / Specificity (%)')
ax1.set_xlim(8, 32)
ax1.set_ylim(785, 945)
ax2.set_ylim(70, 105)
ax1.invert_xaxis()

lns = ln1 + ln2 + ln3
labs = [l.get_label() for l in lns]
ax1.legend(lns, labs, frameon=True, facecolor='white', edgecolor='none',
           loc='lower right', fontsize=7)

ax1.set_title('Figure 3. Classification-efficiency frontier across CV levels',
              fontsize=9, fontweight='bold', loc='left', pad=8)
plt.tight_layout()
save(fig, 'figure3_evsi')

# =====================================================================
# FIGURE 4: Robustness forest plot
# =====================================================================
print("Figure 4: Robustness across distributional violations")

scenarios = [
    '(0) Reference',
    '(A) Log-normal IS',
    '(B) Bimodal τ',
    '(C) Heteroscedastic CV',
    '(D) AR(1) ρ=0.3',
    "(D') AR(1) ρ=0.5",
    '(E) 20% carryover',
    '(F) Combined worst',
]
# Data from robustness_experiments.py
sens_r = [95.1, 96.1, 94.9, 96.1, 94.8, 94.2, 93.6, 92.0]
spec_r = [85.5, 87.9, 86.9, 86.2, 81.9, 80.1, 89.4, 81.1]

fig, ax = plt.subplots(figsize=(5.5, 3.5))
y = np.arange(len(scenarios))

ax.barh(y + 0.15, sens_r, height=0.28, color=C_BLUE, label='Sensitivity (%)', zorder=3)
ax.barh(y - 0.15, spec_r, height=0.28, color=C_ORANGE, label='Specificity (%)', zorder=3)

ax.axvline(x=80, color=C_RED, linewidth=1, linestyle='--', alpha=0.6, zorder=2)
ax.text(80.5, -0.7, 'Threshold\n(80%)', fontsize=6, color=C_RED)

for i in range(len(scenarios)):
    ax.text(sens_r[i] + 0.3, y[i] + 0.15, f'{sens_r[i]:.0f}', fontsize=6, va='center')
    ax.text(spec_r[i] + 0.3, y[i] - 0.15, f'{spec_r[i]:.0f}', fontsize=6, va='center')

ax.set_yticks(y)
ax.set_yticklabels(scenarios)
ax.set_xlabel('Operating characteristic (%)')
ax.set_xlim(75, 102)
ax.legend(frameon=False, loc='lower right', fontsize=7)
ax.invert_yaxis()
ax.set_title('Figure 4. Protocol robustness under distributional stress tests',
             fontsize=9, fontweight='bold', loc='left')
plt.tight_layout()
save(fig, 'figure4_robustness')

print("\nAll figures generated.")
