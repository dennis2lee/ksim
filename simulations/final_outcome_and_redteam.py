"""
STEP 1 (final outcome integration) + RED-TEAM fixes.
Part 1: dialysis-free outcome of the engineered system vs base, across f_tox
        (the unresolved 'do toxins cause death' parameter).
Part 2 (expert critique #4 - pharmacology): AST-120 ADSORBS her antihypertensives.
        Her BP control is the #1 slope lever. Does AST-120's toxin benefit survive the
        penalty of weakening her BP meds? -> keep / separate / drop decision.
Part 3 (expert critique #8 - microbiology): switching Lacto->spore for SURVIVAL may
        sacrifice the documented urea FUNCTION. Effective = survival x function.
"""
import numpy as np
rng = np.random.default_rng(91)
N = 40000

# ---------------- PART 1: outcome across f_tox ----------------
def sim_outcome(engineered, f_tox, tox_red=0.58, horizon=25):
    reach=0; dfree=[]
    for _ in range(N):
        sl = -4.0 + max(rng.normal(2.0,0.5),0)*np.clip(rng.normal(0.9,0.08),0.3,1)   # base medical
        if engineered:
            sl += max(rng.normal(0.5,0.4),0)   # small slope gain (toxin+urea routes)
        sl = min(sl,-0.2)
        eGFR=25.0; t=0.0; out=None
        for yr in range(horizon):
            h = 0.09+0.012*yr
            if engineered:
                h = h*(1-f_tox) + h*f_tox*(1-tox_red)   # toxins->death coupling
            if rng.random()<h: t=yr+rng.random(); out="comp"; break
            eGFR+=sl; t=yr+1
            if eGFR<=8: out="dial"; break
        if out=="dial": reach+=1
        dfree.append(t*12)
    return reach/N, np.median(dfree)

print("="*76)
print("PART 1  Dialysis-free survival across f_tox (toxin->death coupling)")
print("="*76)
prb,dfb = sim_outcome(False,0)
print(f"{'f_tox':>8}{'base dfree':>13}{'engineered dfree':>18}{'gain':>8}{'P(dial) b->e':>16}")
for f in [0.0,0.1,0.2,0.3]:
    pre,dfe = sim_outcome(True,f)
    print(f"{f:>8.1f}{dfb:>12.0f}mo{dfe:>16.0f}mo{('+%.0f'%(dfe-dfb)):>8}{('%.0f%%->%.0f%%'%(prb*100,pre*100)):>16}")
print("  (f_tox=0 is the trial-supported base case: gain ~0; benefit grows only IF toxins cause death)")

# ---------------- PART 2: AST-120 vs her BP meds ----------------
print("\n" + "="*76)
print("PART 2  AST-120 NET slope effect after it binds her antihypertensives")
print("="*76)
AST_TOXIN_GAIN = 0.15       # mL/min/yr from extra IS adsorption
BASE_MED = 2.0              # mL/min/yr from full BP control (the #1 lever)
print(f"{'separation quality':>22}{'BP-med binding':>16}{'BP penalty':>12}{'AST net':>10}")
for label,bind in [("co-administered",0.20),("2 h apart",0.10),("strict >3 h apart",0.04)]:
    penalty = BASE_MED*bind
    net = AST_TOXIN_GAIN - penalty
    verdict = "HARMFUL" if net<0 else "marginal+"
    print(f"{label:>22}{bind*100:>14.0f}%{('-%.2f'%penalty):>12}{('%+.2f'%net):>10}  {verdict}")
print("  -> AST-120 is NET NEGATIVE unless strictly separated from BP meds. Alternative below.")

# ---------------- PART 3: spore vs Lacto effective urea function ----------------
print("\n" + "="*76)
print("PART 3  Probiotic: survival x urea-FUNCTION (switching to spores costs function)")
print("="*76)
options = {
 "Lacto + co-administered AST120": dict(surv=0.41, func=1.00),
 "Spore + co-administered AST120": dict(surv=0.98, func=0.50),
 "Spore + bedtime":                dict(surv=0.98, func=0.50),
 "Lacto + bedtime (8h sep)":       dict(surv=0.91, func=1.00),
 "Lacto + NO AST120 (drop sorbent)":dict(surv=1.00, func=1.00),
}
print(f"{'option':>36}{'survival':>10}{'function':>10}{'effective':>11}")
best=None
for name,o in options.items():
    eff=o["surv"]*o["func"]
    if best is None or eff>best[1]: best=(name,eff)
    print(f"{name:>36}{o['surv']*100:>9.0f}%{o['func']*100:>9.0f}%{eff*100:>10.0f}%")
print(f"\n  Best effective urea function: '{best[0]}' ({best[1]*100:.0f}%)")
print("  -> Earlier 'switch to spores' was WRONG once function is counted: Lacto+bedtime")
print("     (or dropping AST-120) delivers more urea metabolism than spores.")
