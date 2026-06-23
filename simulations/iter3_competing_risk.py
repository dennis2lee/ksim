"""
ITERATION 3 of autonomous loop. Two honesty fixes + a decision output.
  (1) COMPETING RISK: a 74yo with stage-4 CKD faces real annual hazard of death /
      CV event / AKI-forced dialysis the gut approach cannot touch. This caps the
      absurd upper tails and reframes the goal as 'dialysis-FREE survival'.
  (2) DROP-ONE value decomposition: which creative arm actually carries the benefit?
Year-by-year simulation: eGFR declines; each year a competing event may intervene;
dialysis when eGFR<8.
"""
import numpy as np
rng = np.random.default_rng(31)
N = 40000

CAT = {
 "base_medical":    dict(slope=(2.0,0.25), adh=0.10, feas=0.97),
 "combined_sachet": dict(slope=(0.50,0.55),adh=0.12, feas=0.80, se=("GI",0.15)),
 "natural_urea_pro":dict(slope=(0.28,0.8), adh=0.18, feas=0.85),
 "keto_lowprotein": dict(slope=(0.42,0.5), adh=0.18, feas=0.82),
 "Ntrap_fiber":     dict(slope=(0.22,0.7), adh=0.25, feas=0.90, se=("bloat",0.18)),
}
def eff_slope(name):
    s=CAT[name]; g=max(rng.normal(s["slope"][0], s["slope"][0]*s["slope"][1]),0)
    g*=np.clip(rng.normal(1-s["adh"],0.08),0.3,1.0)
    if s.get("se") and rng.random()<s["se"][1]: g*=0.3
    return g
def regimen_slope(names):
    gains=[]
    for n in names:
        if rng.random()>CAT[n]["feas"]: continue
        gains.append(eff_slope(n))
    gains.sort(reverse=True)
    return -4.0 + sum(g*(0.8**i) for i,g in enumerate(gains))

# competing annual hazard for 74yo CKD4 (death / CV / AKI-forced dialysis), rising slowly
def competing_event(age_offset):
    h = 0.09 + 0.012*age_offset            # ~9%/yr now, climbing with age
    return rng.random() < h

def simulate(names, horizon_yr=25):
    reach_dialysis=0; dfree_times=[]; dialysis_times=[]
    for _ in range(N):
        eGFR=25.0; sl=min(regimen_slope(names),-0.2); t=0.0; outcome=None
        for yr in range(horizon_yr):
            if competing_event(yr):
                outcome="competing"; t=yr+rng.random(); break
            eGFR += sl
            t = yr+1
            if eGFR<=8:
                outcome="dialysis"; break
        if outcome=="dialysis":
            reach_dialysis+=1; dialysis_times.append(t*12)
        dfree_times.append(t*12)   # months lived dialysis-free (until dialysis OR competing event)
    return (reach_dialysis/N, np.median(dialysis_times) if dialysis_times else np.nan,
            np.median(dfree_times), np.percentile(dfree_times,[5,95]))

regimens = {
 "Base only":                 ["base_medical"],
 "CREATIVE deployable (best)":["base_medical","combined_sachet","natural_urea_pro","keto_lowprotein","Ntrap_fiber"],
}
print("="*86)
print("ITER 3  WITH COMPETING RISK (74yo)   P(ever reach dialysis) | dialysis-free survival")
print("="*86)
for name,sel in regimens.items():
    pr,dt,dfree,ci = simulate(sel)
    dtxt = f"{dt:.0f} mo" if not np.isnan(dt) else "n/a"
    print(f"{name:<28} P(dialysis)={pr*100:>4.0f}%   dialysis-free median={dfree:>4.0f} mo [{ci[0]:.0f}-{ci[1]:.0f}]   (if reached: {dtxt})")

# ---- DROP-ONE value decomposition on the best regimen ----
print("\n" + "="*86)
print("DROP-ONE: contribution of each creative arm to dialysis-free survival (months)")
print("="*86)
best = ["base_medical","combined_sachet","natural_urea_pro","keto_lowprotein","Ntrap_fiber"]
_,_,full_df,_ = simulate(best)
for n in best[1:]:
    reduced=[x for x in best if x!=n]
    _,_,df,_ = simulate(reduced)
    print(f"removing {n:<20} -> dialysis-free {df:>4.0f} mo   (loses {full_df-df:>4.0f} mo)")
print(f"\nfull best regimen dialysis-free median = {full_df:.0f} mo")
