"""
ITERATION 2 of autonomous loop.
The infeasible high-benefit arms (novel molecule, engineered strain, Na-alkali,
aggressive protein cut) are each REPLACED by a creative, deployable analogue that
dodges that specific feasibility wall. Question: can we recover most of the
theoretical +16-20 month dialysis delay at a REAL P(deploy)?
"""
import numpy as np
rng = np.random.default_rng(23)
N = 30000

CAT = {
 # --- established base ---
 "base_medical":   dict(slope=(2.0,0.25), adh=0.10, feas=0.97, se=None,  note="BP+RAS+SGLT2i"),
 # --- creative bypass arms ---
 "combined_sachet":dict(slope=(0.50,0.55),adh=0.12, feas=0.80, se=("GI",0.15),
                        note="BYPASS AST120 burden: fiber+selective sorbent, 1/day sachet -> high adherence"),
 "natural_urea_pro":dict(slope=(0.28,0.8),adh=0.18, feas=0.85, se=None,
                        note="BYPASS engineered strain: NATURAL urea-consuming synbiotic (Renadyl-type), OTC"),
 "Ntrap_fiber":    dict(slope=(0.22,0.7), adh=0.25, feas=0.90, se=("bloat",0.18),
                        note="BYPASS strain: high-dose fermentable fiber -> bacterial N trapped in stool"),
 "keto_lowprotein":dict(slope=(0.42,0.5), adh=0.18, feas=0.82, se=None,
                        note="BYPASS frailty: ketoacid-analogue + 0.6g/kg protein -> urea down, muscle kept"),
 "veverimer_alkali":dict(slope=(0.45,0.5),adh=0.15, feas=0.45, se=None,
                        note="BYPASS Na load: non-absorbed acid-binding polymer (investigational, KR access iffy)"),
 # --- the OLD infeasible originals, for comparison ---
 "novel_molecule": dict(slope=(0.20,0.9), adh=0.20, feas=0.05, se=None, note="orig novel drug (infeasible 1yr)"),
}

def eff_slope(name):
    s = CAT[name]
    g = max(rng.normal(s["slope"][0], s["slope"][0]*s["slope"][1]), 0)
    g *= np.clip(rng.normal(1 - s["adh"], 0.08), 0.3, 1.0)
    if s["se"] and rng.random() < s["se"][1]:
        g *= 0.3
    return g

def regimen_slope(names):
    gains = []
    for n in names:
        if rng.random() > CAT[n]["feas"]:
            continue
        gains.append(eff_slope(n))
    gains.sort(reverse=True)
    total = sum(g * (0.8**i) for i, g in enumerate(gains))
    return -4.0 + total

def delay_months(names):
    out = []
    for _ in range(N):
        sl = min(regimen_slope(names), -0.2)
        out.append((25.0-8.0)/(-sl)*12)
    return np.array(out)

regimens = {
 "Base only":                       ["base_medical"],
 "OLD infeasible ALL":              ["base_medical","combined_sachet","natural_urea_pro","novel_molecule"],
 "CREATIVE deployable (core)":      ["base_medical","combined_sachet","natural_urea_pro","keto_lowprotein"],
 "CREATIVE + N-trap fiber":         ["base_medical","combined_sachet","natural_urea_pro","keto_lowprotein","Ntrap_fiber"],
 "CREATIVE + veverimer (stretch)":  ["base_medical","combined_sachet","natural_urea_pro","keto_lowprotein","Ntrap_fiber","veverimer_alkali"],
}

print("="*88)
print("ITER 2  TIME TO DIALYSIS  median [90% CI] | delay vs base | P(deploy) | risk-adj delay")
print("="*88)
ref = np.median(delay_months(["base_medical"]))
rows=[]
for name, sel in regimens.items():
    d = delay_months(sel); md=np.median(d); ci=np.percentile(d,[5,95])
    feas = float(np.prod([CAT[n]["feas"] for n in sel]))
    delay = md-ref
    radj = delay*feas       # risk-adjusted: benefit you actually expect to realize
    rows.append((name,md,ci,delay,feas,radj))
    print(f"{name:<32}{md:>5.0f} [{ci[0]:>3.0f}-{ci[1]:>3.0f}]  {('+%.0f'%delay):>5} mo  P~{feas:.2f}   {radj:>5.0f} mo")

best = max(rows, key=lambda r: r[5])
print("\n--> Highest RISK-ADJUSTED delay:", best[0],
      f"(+{best[3]:.0f} mo nominal, P~{best[4]:.2f}, risk-adj +{best[5]:.0f} mo)")
