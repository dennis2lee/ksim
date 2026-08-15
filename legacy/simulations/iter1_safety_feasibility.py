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
ITERATION 1 of autonomous loop.
Adds three layers the efficacy model lacked:
  (A) SIDE EFFECTS -> encoded as adherence penalties + hard safety gates
  (B) FEASIBILITY  -> P(deployable to mother within ~1 yr) from regulatory status x safety
  (C) DIALYSIS DELAY -> translate toxin/slope effects into months of delay (what matters)

Patient: eGFR 25 (stage 3b->4), hypertensive, non-diabetic, 74F.
Dialysis threshold modeled at eGFR 8. Untreated slope ~ -4 mL/min/yr.
All numbers illustrative for hypothesis-ranking, NOT clinical promises.
"""
import numpy as np
rng = np.random.default_rng(11)
N = 30000

# ---------------------------------------------------------------------------
# Intervention catalog. Each:
#   slope_gain  = mL/min/yr decline AVERTED (mean, cv)   [bigger = slows CKD more]
#   adh_penalty = baseline non-adherence drag (0..1)     [GI burden, pill load]
#   feas        = P(deployable within 1 yr to mother)    [regulatory+safety]
#   tox_se      = qualitative safety flag for gating
# ---------------------------------------------------------------------------
CAT = {
 "base_medical":   dict(slope=(2.0,0.25), adh=0.10, feas=0.97, note="BP+RAS+SGLT2i: established, dominant lever"),
 "acidosis_bicarb":dict(slope=(0.4,0.5),  adh=0.15, feas=0.90, note="NaHCO3: slope benefit BUT sodium load vs HTN"),
 "fiber_prebiotic":dict(slope=(0.25,0.7), adh=0.25, feas=0.95, note="OTC; GI bloating limits dose"),
 "AST120":         dict(slope=(0.30,0.8), adh=0.45, feas=0.70, note="Rx in KR; high pill burden -> poor adherence; binds other drugs"),
 "synbiotic":      dict(slope=(0.10,0.9), adh=0.20, feas=0.92, note="OTC; weak evidence"),
 "protein_mod":    dict(slope=(0.30,0.6), adh=0.30, feas=0.85, note="frailty/sarcopenia risk in 74F -> capped"),
 "TnaA_inhibitor": dict(slope=(0.20,0.9), adh=0.20, feas=0.05, note="NOVEL molecule: years of trials; not deployable in 1yr"),
 "urease_strain":  dict(slope=(0.35,0.9), adh=0.20, feas=0.03, note="engineered live strain: regulatory+infection wall"),
}

# Side-effect HARD GATES (a sampled adverse event removes that arm's benefit this draw)
SE_EVENT = {
 "acidosis_bicarb": ("fluid/BP worsening (Na load)", 0.20),
 "AST120":          ("GI intolerance / drug malabsorption", 0.35),
 "protein_mod":     ("protein-energy wasting / frailty", 0.30),
 "fiber_prebiotic": ("bloating -> self-discontinue", 0.20),
}
# Potassium gate: RAS blocker + any fiber load -> small hyperkalemia risk that caps add-ons
K_RISK = 0.12

def eff_slope(name):
    """Sampled averted-decline for one intervention, after adherence + side-effect gate."""
    s = CAT[name]
    g = max(rng.normal(*[s["slope"][0], s["slope"][0]*s["slope"][1]]), 0)
    # adherence drag
    a = np.clip(rng.normal(1 - s["adh"], 0.08), 0.3, 1.0)
    g *= a
    # side-effect hard gate
    if name in SE_EVENT and rng.random() < SE_EVENT[name][1]:
        g *= 0.2   # event -> mostly lose benefit (dose cut / stop)
    return g

def regimen_slope(names, apply_feas=True):
    """Total averted decline (mL/min/yr) for a regimen, with diminishing returns on add-ons."""
    base = -4.0
    gains = []
    for n in names:
        if apply_feas and rng.random() > CAT[n]["feas"]:
            continue   # not deployable this draw
        gains.append(eff_slope(n))
    gains.sort(reverse=True)
    total = 0.0
    for i, g in enumerate(gains):
        total += g * (0.8 ** i)   # each additional add-on ~80% credit (overlap/ceiling)
    # potassium ceiling: if many K-raising arms stacked, occasionally cap
    if len([n for n in names if n in ("fiber_prebiotic","synbiotic","AST120")]) >= 2 and rng.random() < K_RISK:
        total *= 0.7
    return base + total   # less negative = slower decline

def dialysis_delay_months(names):
    eGFR0, thr = 25.0, 8.0
    out = []
    for _ in range(N):
        sl = regimen_slope(names)
        sl = min(sl, -0.2)  # never let it go non-progressive
        yrs = (eGFR0 - thr) / (-sl)
        out.append(yrs*12)
    return np.array(out)

regimens = {
 "Base medical only":              ["base_medical"],
 "Base + gut stack":               ["base_medical","fiber_prebiotic","AST120","synbiotic"],
 "Base + gut + acidosis":          ["base_medical","fiber_prebiotic","AST120","synbiotic","acidosis_bicarb"],
 "Base + gut + protein mod":       ["base_medical","fiber_prebiotic","AST120","synbiotic","protein_mod"],
 "Base + ALL (incl novel)":        ["base_medical","fiber_prebiotic","AST120","synbiotic","acidosis_bicarb","protein_mod","TnaA_inhibitor","urease_strain"],
}

print("="*84)
print("TIME TO DIALYSIS (months from now)  median [90% CI]   |  delay vs base  |  feasibility")
print("="*84)
ref = np.median(dialysis_delay_months(["base_medical"]))
for name, sel in regimens.items():
    d = dialysis_delay_months(sel)
    md = np.median(d); ci = np.percentile(d,[5,95])
    delay = md - ref
    feas = np.prod([CAT[n]["feas"] for n in sel])
    print(f"{name:<30}{md:>6.0f} [{ci[0]:>3.0f}-{ci[1]:>3.0f}] mo   {('+%.0f mo'%delay):>9}     P(deploy)~{feas:.2f}")

print("\n" + "="*84)
print("WHERE FEASIBILITY COLLAPSES (the wall)")
print("="*84)
for n in CAT:
    bar = "#"*int(CAT[n]["feas"]*30)
    print(f"{n:<18}{CAT[n]['feas']:.2f} {bar:<30} {CAT[n]['note']}")
