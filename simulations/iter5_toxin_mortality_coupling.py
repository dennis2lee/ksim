"""
ITERATION 5. Tests the user's challenge: if lowering uremic toxins ALSO lowers the
competing (death) hazard -- i.e. toxins and death are NOT independent -- does the
gut/toxin route start to extend dialysis-free survival after all?

Key parameter f_tox = fraction of the competing hazard that is TOXIN-DRIVEN.
  f_tox = 0.0  -> my original assumption (toxins independent of death; AST-120 trial-like)
  f_tox > 0    -> toxins drive some CV death; gut stack (which cuts toxins ~50%) then
                  also cuts that share of the death hazard => a SHARED factor.
We sweep f_tox and watch the gut-stack benefit appear.
"""
import numpy as np
rng = np.random.default_rng(53)
N = 40000

# gut stack lumped: slows CKD slope AND cuts toxins by ~50%
GUT_SLOPE = (1.0, 0.5)     # mL/min/yr averted (mean, cv)
GUT_TOXIN_CUT = 0.50       # fraction reduction in uremic toxins
BASE_SLOPE = -4.0
BASE_MED_GAIN = (2.0, 0.25)

def draw(mean_cv):
    return max(rng.normal(mean_cv[0], mean_cv[0]*mean_cv[1]), 0)

def simulate(use_gut, f_tox, horizon=25):
    reach=0; dfree=[]
    for _ in range(N):
        sl = BASE_SLOPE + draw(BASE_MED_GAIN)*np.clip(rng.normal(0.9,0.08),0.3,1)
        if use_gut:
            sl += draw(GUT_SLOPE)*np.clip(rng.normal(0.8,0.08),0.3,1)
        sl = min(sl, -0.2)
        eGFR=25.0; t=0.0; outcome=None
        for yr in range(horizon):
            h = 0.09 + 0.012*yr
            # split hazard: toxin-driven part can be reduced by the gut stack
            if use_gut:
                h = h*(1 - f_tox) + h*f_tox*(1 - GUT_TOXIN_CUT)
            if rng.random() < h:
                t=yr+rng.random(); outcome="comp"; break
            eGFR += sl; t=yr+1
            if eGFR<=8: outcome="dial"; break
        if outcome=="dial": reach+=1
        dfree.append(t*12)
    return reach/N, np.median(dfree)

print("="*78)
print("ITER 5  Does the gut/toxin route extend dialysis-free survival IF toxins")
print("        also drive death?  (sweep f_tox = toxin-driven share of mortality)")
print("="*78)
print(f"{'f_tox':>8}{'base dfree':>13}{'+gut dfree':>13}{'gut gain':>11}{'P(dial) base->gut':>22}")
for f_tox in [0.0, 0.1, 0.2, 0.3, 0.5]:
    pr_b, df_b = simulate(False, f_tox)
    pr_g, df_g = simulate(True,  f_tox)
    print(f"{f_tox:>8.1f}{df_b:>12.0f}{df_g:>13.0f}{('+%.0f mo'%(df_g-df_b)):>11}{('%.0f%%->%.0f%%'%(pr_b*100,pr_g*100)):>22}")

print("\nf_tox=0  -> my original model (gut route adds ~0-2 mo).")
print("f_tox>0  -> if toxins truly cause CV death, the gut route DOES extend the clock.")
print("This is the load-bearing assumption. Current trial evidence ~ f_tox near 0.")
