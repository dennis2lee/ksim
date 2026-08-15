# ARCHIVED. NOT PART OF THE PAPER.
#
# validation/reproduce_all.py does not run this file, and no result in
# the manuscript or its supplement depends on it. It is kept as a record
# of earlier work and is not maintained.
#
# This file is exploratory work on toxin models and intervention
# engineering. It informed no number in the paper.
#
# For what the paper actually does, read validation/nof1_core.py and the
# thirteen scripts the README lists under "What the paper reproduces".

"""
Gut-clearance ("enteric dialysis") simulation for protein-bound uremic toxins (PBUTs)
in a stage-4 CKD patient (eGFR ~25, hypertensive, non-diabetic, 74F).

PURPOSE: hypothesis-exploration model, NOT a clinical dosing tool.
We model serum steady-state and trajectory of three wastes:
  - Indoxyl sulfate (IS)   : tryptophan -> gut bacterial tryptophanase -> indole -> liver sulfation
  - p-Cresyl sulfate (PCS) : tyrosine/phe -> gut bacteria -> p-cresol -> sulfation
  - Urea (BUN)             : protein nitrogen -> urea

Core equation (per toxin):
    serum_relative = G_factor / CL_factor
where G_factor = product of intervention generation-reduction multipliers,
and CL_factor = clearance relative to baseline (residual renal + small gut sink).

Intervention efficacies are mid-range estimates anchored to published trial signals
(AST-120 ~30-40% IS drop; prebiotic fiber ~20-30% IS/PCS; synbiotic modest;
low-protein diet reduces nitrogenous substrate). Treat as ORDER-OF-MAGNITUDE.
"""

import numpy as np

# ---------------------------------------------------------------------------
# 1. Intervention library: per-toxin GENERATION multipliers (1.0 = no effect)
#    Lower = more reduction. Mechanism noted for overlap handling later.
# ---------------------------------------------------------------------------
INTERVENTIONS = {
    # name: (IS_mult, PCS_mult, urea_mult, mechanism_tag, note)
    "protein_moderation_0.6": (0.80, 0.80, 0.70, "substrate",
        "Protein 0.8->0.6 g/kg/d: cuts Trp/Tyr/N substrate. Frailty caution in elderly."),
    "fiber_prebiotic":        (0.75, 0.68, 0.90, "microbiota_shift",
        "+25-40 g/d fiber/inulin/resistant starch: saccharolytic shift; PCS very fiber-sensitive."),
    "AST120_sorbent":         (0.62, 0.75, 0.98, "lumen_sorbent",
        "AST-120 6-9 g/d: adsorbs indole/p-cresol in lumen before absorption; urea not adsorbed."),
    "synbiotic_renadyl":      (0.88, 0.84, 0.90, "microbiota_shift",
        "Probiotic+prebiotic: consumes N precursors, makes SCFA. Modest, overlaps with fiber."),
    # --- NOVEL drug targets (hypothetical R&D candidates) ---
    "NOVEL_tryptophanase_inh":(0.40, 1.00, 1.00, "enzyme_block",
        "Gut tryptophanase (TnaA) inhibitor: blocks indole AT SOURCE. IS-specific, strong."),
    "NOVEL_urease_strain":    (1.00, 1.00, 0.60, "engineered_strain",
        "Engineered urea-consuming strain: diverts urea-N to gut. Urea-specific."),
}

# ---------------------------------------------------------------------------
# 2. Combine interventions. Independent mechanisms multiply; interventions
#    sharing a mechanism_tag get a diminishing-returns dampening (no naive stacking).
# ---------------------------------------------------------------------------
def combine(selected, toxin_idx):
    by_mech = {}
    for name in selected:
        is_m, pcs_m, urea_m, mech, _ = INTERVENTIONS[name]
        mult = (is_m, pcs_m, urea_m)[toxin_idx]
        by_mech.setdefault(mech, []).append(mult)

    total = 1.0
    for mech, mults in by_mech.items():
        # within same mechanism: take the best, then add only damped extra from others
        mults_sorted = sorted(mults)  # smallest (strongest) first
        eff = mults_sorted[0]
        for m in mults_sorted[1:]:
            extra_reduction = (1 - m) * 0.5   # overlapping -> only half-credit
            eff *= (1 - extra_reduction)
        total *= eff
    return total

# ---------------------------------------------------------------------------
# 3. Residual renal clearance scaling (stage 4 baseline eGFR=25).
#    PBUTs are >90% albumin-bound -> cleared mainly by tubular secretion,
#    which falls FASTER than GFR. We model CL ~ (eGFR/25)^1.2 for PBUTs.
# ---------------------------------------------------------------------------
def cl_factor(eGFR, baseline=25.0, exp=1.2):
    return (eGFR / baseline) ** exp

def serum_relative(selected, eGFR=25.0):
    out = {}
    for i, tox in enumerate(["IS", "PCS", "urea"]):
        G = combine(selected, i)
        CL = cl_factor(eGFR) if tox != "urea" else (eGFR / 25.0)  # urea ~ linear w/ GFR
        out[tox] = 100.0 * G / CL   # relative to untreated baseline at eGFR 25
    return out

# ---------------------------------------------------------------------------
# 4. Scenarios
# ---------------------------------------------------------------------------
scenarios = {
    "Baseline (no Tx)":                 [],
    "Fiber only":                       ["fiber_prebiotic"],
    "AST-120 only":                     ["AST120_sorbent"],
    "Synbiotic only":                   ["synbiotic_renadyl"],
    "Protein moderation only":          ["protein_moderation_0.6"],
    "STACK: fiber+AST120":              ["fiber_prebiotic", "AST120_sorbent"],
    "STACK: fiber+AST120+synbiotic":    ["fiber_prebiotic", "AST120_sorbent", "synbiotic_renadyl"],
    "STACK + protein moderation":       ["fiber_prebiotic", "AST120_sorbent", "synbiotic_renadyl", "protein_moderation_0.6"],
    "ADD novel tryptophanase inh":      ["fiber_prebiotic", "AST120_sorbent", "NOVEL_tryptophanase_inh"],
    "FULL incl. novel targets":         ["fiber_prebiotic", "AST120_sorbent", "synbiotic_renadyl",
                                         "protein_moderation_0.6", "NOVEL_tryptophanase_inh", "NOVEL_urease_strain"],
}

print("="*78)
print("SERUM UREMIC TOXIN LEVEL (relative %, baseline untreated stage-4 = 100)")
print("eGFR = 25 (stage 4). Lower = better.")
print("="*78)
print(f"{'Scenario':<34}{'IS':>8}{'PCS':>8}{'urea':>8}   {'mean':>6}")
print("-"*78)
for name, sel in scenarios.items():
    r = serum_relative(sel)
    mean = np.mean([r['IS'], r['PCS'], r['urea']])
    print(f"{name:<34}{r['IS']:>8.0f}{r['PCS']:>8.0f}{r['urea']:>8.0f}   {mean:>6.0f}")

# ---------------------------------------------------------------------------
# 5. Single-factor sensitivity (which lever moves IS the most?) at eGFR 25
# ---------------------------------------------------------------------------
print("\n" + "="*78)
print("SENSITIVITY: single-intervention reduction in INDOXYL SULFATE (IS)")
print("="*78)
base_IS = serum_relative([])['IS']
rows = []
for name in INTERVENTIONS:
    is_level = serum_relative([name])['IS']
    rows.append((name, base_IS - is_level))
for name, drop in sorted(rows, key=lambda x: -x[1]):
    bar = "#" * int(drop)
    print(f"{name:<28}{drop:>5.0f}%  {bar}")

# ---------------------------------------------------------------------------
# 6. eGFR decline: how serum IS rises as kidney fails, treated vs untreated.
#    Shows the "buffer" gut-clearance buys as residual function drops.
# ---------------------------------------------------------------------------
print("\n" + "="*78)
print("IS level vs declining eGFR  (untreated  vs  realistic stack)")
print("="*78)
stack = ["fiber_prebiotic", "AST120_sorbent", "synbiotic_renadyl"]
print(f"{'eGFR':>6}{'untreated IS':>16}{'stacked IS':>14}")
for egfr in [30, 25, 20, 15, 12]:
    u = serum_relative([], egfr)['IS']
    s = serum_relative(stack, egfr)['IS']
    print(f"{egfr:>6}{u:>16.0f}{s:>14.0f}")

print("\nNote: 'stacked' at eGFR 20 reaches a similar IS level to 'untreated' at a")
print("higher eGFR -> interpretable as buying back equivalent toxin-clearance.")
