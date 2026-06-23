"""
DEV component (3): optimize protein target with/without ketoacid analogues.
Tradeoff: lower protein -> less urea (good) BUT more frailty/sarcopenia in a 74F (bad).
Ketoacid analogues supply essential amino-acid skeletons WITHOUT nitrogen, so they let
protein drop further while protecting muscle -> they shift the frailty curve down.
Find the protein target that minimizes a combined urea+frailty objective.
"""
import numpy as np

p = np.linspace(0.40, 1.00, 61)   # protein g/kg/day

# urea generation (ref at 0.8 g/kg = 1.0); endogenous floor 0.2
urea = (0.2 + p) / 1.0

def frailty(p, keto):
    eff = p + (0.25 if keto else 0.0)   # ketoacid adds ~0.25 g/kg-equivalent EAA, no N
    return 1/(1+np.exp(12*(eff-0.72)))  # rises sharply as effective intake < ~0.72

W_UREA, W_FRAIL = 1.0, 1.3   # frailty weighted a bit higher (74yo, irreversible)

def optimize(keto):
    obj = W_UREA*urea + W_FRAIL*frailty(p, keto)
    i = int(np.argmin(obj))
    return p[i], urea[i], frailty(p[i], keto), obj[i]

print("="*72)
print("DEV-(3)  OPTIMAL PROTEIN TARGET  (minimize urea + frailty)")
print("="*72)
for keto in [False, True]:
    popt,u,f,o = optimize(keto)
    tag = "WITH ketoacid" if keto else "NO ketoacid  "
    print(f"  {tag}: protein* = {popt:.2f} g/kg | urea {u*100:>3.0f}% of ref | frailty risk {f*100:>3.0f}%")

pn,un,fn,_ = optimize(False)
pk,uk,fk,_ = optimize(True)
print(f"\n  Ketoacid lets protein drop {pn:.2f} -> {pk:.2f} g/kg safely")
print(f"  -> extra urea reduction: {(un-uk)*100:.0f} percentage points, frailty kept ~{fk*100:.0f}%")

# dosing for a ~58 kg patient (typical 74F); Ketosteril ~1 tab / 5 kg / day
WT = 58
print(f"\n  For {WT} kg: protein target ~{pk*WT:.0f} g/day; ketoacid ~{WT/5:.0f} tabs/day with meals")
print("  Monitor: serum calcium (Ca-salt tabs), albumin, body weight, grip strength.")

# show the tradeoff curve numerically
print("\n  TRADEOFF TABLE (with ketoacid):")
print(f"  {'protein':>8}{'urea %ref':>11}{'frailty %':>11}")
for pv in [0.4,0.5,0.6,0.7,0.8]:
    fv = 1/(1+np.exp(12*((pv+0.25)-0.72)))
    print(f"  {pv:>8.2f}{((0.2+pv))*100:>11.0f}{fv*100:>11.0f}")
