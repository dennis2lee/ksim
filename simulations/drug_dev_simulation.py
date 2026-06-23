"""
Stage 2 of the gut-clearance drug-development simulation.
Builds on gut_clearance_model.py with three additions the single-point model lacked:

  (1) Monte Carlo uncertainty  -> turn point estimates into 90% credible intervals
  (2) Adherence sensitivity    -> real-world effect when pills are missed
  (3) Novel drug dose-response -> how hard must a gut tryptophanase (TnaA) inhibitor
                                  hit to beat / add to the existing fiber+AST120 stack

Outputs feed the go/no-go question: is a NEW molecule worth it, or does the
already-buyable stack capture most of the achievable benefit?
"""

import numpy as np
rng = np.random.default_rng(7)

# ---------------------------------------------------------------------------
# Point estimates: per-toxin GENERATION reduction r0 = 1 - multiplier
# and a relative uncertainty (CV). (IS, PCS, urea), mechanism tag.
# ---------------------------------------------------------------------------
LIB = {
    "fiber_prebiotic":  dict(r=(0.25, 0.32, 0.10), cv=0.40, mech="microbiota"),
    "AST120":           dict(r=(0.38, 0.25, 0.02), cv=0.35, mech="sorbent"),
    "synbiotic":        dict(r=(0.12, 0.16, 0.10), cv=0.50, mech="microbiota"),
    "protein_mod":      dict(r=(0.20, 0.20, 0.30), cv=0.30, mech="substrate"),
    # novel agents
    "TnaA_inhibitor":   dict(r=(0.66, 0.00, 0.00), cv=0.45, mech="enzyme"),   # IS-specific, at full dose
    "urease_strain":    dict(r=(0.00, 0.00, 0.40), cv=0.50, mech="strain"),
}
TOX = ["IS", "PCS", "urea"]

def draw_mult(name, toxin_idx, adherence=1.0):
    """Sample one generation multiplier (=1-reduction) for an intervention/toxin."""
    spec = LIB[name]
    r0 = spec["r"][toxin_idx]
    if r0 <= 0:
        return 1.0
    sd = spec["cv"] * r0
    r = np.clip(rng.normal(r0, sd), 0.0, 0.95)
    r *= adherence                      # missed doses scale the effect down
    return 1.0 - r

def combine(names, toxin_idx, adherence=1.0):
    """Combine multipliers; same-mechanism interventions get half-credit overlap."""
    by_mech = {}
    for n in names:
        m = draw_mult(n, toxin_idx, adherence)
        by_mech.setdefault(LIB[n]["mech"], []).append(m)
    total = 1.0
    for mults in by_mech.values():
        mults.sort()
        eff = mults[0]
        for m in mults[1:]:
            eff *= (1 - (1 - m) * 0.5)
        total *= eff
    return total

def cl_factor(eGFR, base=25.0, exp=1.2):
    return (eGFR / base) ** exp

def serum(names, toxin_idx, eGFR=25.0, adherence=1.0):
    G = combine(names, toxin_idx, adherence)
    CL = cl_factor(eGFR) if toxin_idx != 2 else (eGFR / 25.0)
    return 100.0 * G / CL

N = 20000
def mc(names, toxin="IS", eGFR=25.0, adherence=1.0):
    ti = TOX.index(toxin)
    s = np.array([serum(names, ti, eGFR, adherence) for _ in range(N)])
    return s

# ===========================================================================
# (1) MONTE CARLO: realistic stack vs full(+novel), with 90% CI
# ===========================================================================
stacks = {
    "Realistic stack (fiber+AST120+synbiotic)": ["fiber_prebiotic","AST120","synbiotic"],
    "+ protein moderation":                     ["fiber_prebiotic","AST120","synbiotic","protein_mod"],
    "+ NOVEL TnaA inhibitor":                   ["fiber_prebiotic","AST120","synbiotic","TnaA_inhibitor"],
    "FULL (+TnaA +urease strain +protein)":     ["fiber_prebiotic","AST120","synbiotic","protein_mod","TnaA_inhibitor","urease_strain"],
}
print("="*82)
print("MONTE CARLO  (N=20,000)  serum level vs untreated=100, median [90% CI], lower=better")
print("="*82)
print(f"{'Stack':<44}{'IS median[90% CI]':>20}{'PCS':>9}{'urea':>9}")
print("-"*82)
for name, sel in stacks.items():
    isd = mc(sel,"IS"); pc = mc(sel,"PCS"); ur = mc(sel,"urea")
    ci = np.percentile(isd,[5,95])
    print(f"{name:<44}{np.median(isd):>7.0f} [{ci[0]:>3.0f}-{ci[1]:>3.0f}]{np.median(pc):>9.0f}{np.median(ur):>9.0f}")

# ===========================================================================
# (2) ADHERENCE: IS level for realistic stack as adherence drops 100->50%
# ===========================================================================
print("\n" + "="*82)
print("ADHERENCE SENSITIVITY  (realistic stack, IS, median)")
print("="*82)
base = ["fiber_prebiotic","AST120","synbiotic"]
print(f"{'adherence':>10}{'IS median':>12}{'reduction':>12}")
for a in [1.0, 0.9, 0.8, 0.7, 0.6, 0.5]:
    isd = mc(base,"IS", adherence=a)
    print(f"{a*100:>9.0f}%{np.median(isd):>12.0f}{100-np.median(isd):>11.0f}%")

# ===========================================================================
# (3) NOVEL DRUG DOSE-RESPONSE: TnaA inhibition theta vs IS, on fiber+AST120 base
#     E_max=0.70 (not all indole is TnaA/dietary-independent). Compare to the
#     "buffer target": untreated IS at eGFR 30 = 80 (i.e. claw back ~ one stage).
# ===========================================================================
print("\n" + "="*82)
print("NOVEL TnaA-INHIBITOR DOSE-RESPONSE  (on fiber+AST120 background, IS, eGFR=25)")
print("="*82)
E_MAX = 0.70
bg = ["fiber_prebiotic","AST120"]
target = 80.0  # untreated IS at eGFR 30
print(f"{'enzyme inhib θ':>14}{'IS median':>12}{'vs target 80':>14}")
for theta in [0.0, 0.2, 0.4, 0.6, 0.8, 0.95]:
    # express the inhibitor as an intervention with reduction E_MAX*theta
    LIB["_dose"] = dict(r=(E_MAX*theta, 0.0, 0.0), cv=0.45, mech="enzyme")
    isd = mc(bg + (["_dose"] if theta>0 else []), "IS")
    med = np.median(isd)
    flag = "reached" if med <= target else ""
    print(f"{theta:>14.2f}{med:>12.0f}{('  '+flag):>14}")

print("\nInterpretation hooks printed above; see response for the go/no-go read.")
