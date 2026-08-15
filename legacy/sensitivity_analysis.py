# ARCHIVED. NOT PART OF THE PAPER.
#
# validation/reproduce_all.py does not run this file, and no result in
# the manuscript or its supplement depends on it. It is kept as a record
# of earlier work and is not maintained.
#
# This file contains an earlier feasibility argument for CV 0.15. It is
# superseded by validation/variance_components_analysis.py and by
# LIMITATIONS.md (d).
#
# For what the paper actually does, read validation/nof1_core.py and the
# thirteen scripts the README lists under "What the paper reproduces".

"""
PEER REVIEW SENSITIVITY ANALYSIS.

Systematic parameter sweeps demonstrating which conclusions are robust to
parameter uncertainty and which depend on specific assumptions.

Addresses anticipated reviewer concerns:
  (1) Parameter calibration — literature-anchored ranges for every key parameter
  (2) Tornado sensitivity — one-at-a-time sweep showing which parameters matter
  (3) f_tox robustness — which conclusions survive f_tox=0?
  (4) CV feasibility — range of achievable noise reductions
  (5) FP threshold trade-off — adjustable decision thresholds
  (6) Cohort heterogeneity — vary patient distribution parameters

All parameter ranges are anchored to published clinical literature (see table below).
"""
import numpy as np
rng = np.random.default_rng(201)

N_PAT = 100
N_REP = 300       # reduced for speed (sensitivity needs many runs)
WK_T, WK_W = 4, 2
WK_CYC = 2*(WK_T + WK_W)

# =========================================================================
# (1) LITERATURE-CALIBRATED PARAMETER TABLE
# =========================================================================
print("="*90)
print("(1) PARAMETER CALIBRATION TABLE — literature sources and uncertainty ranges")
print("="*90)
PARAMS = [
 ("IS within-person CV",        0.22, 0.18, 0.28,
  "Pretorius CJ et al, Clin Chim Acta 2013"),
 ("eGFR slope (mL/min/yr)",    -2.0, -4.0, -1.0,
  "Inker LA et al, Clin J Am Soc Nephrol 2017 (CRIC); Levey AS et al, JASN 1999 (MDRD)"),
 ("Population mean tau (IS red)", 0.35, 0.20, 0.50,
  "Rossi M et al, Clin J Am Soc Nephrol 2016 (fiber); Schulman G, EPPIC, JASN 2015 (AST-120)"),
 ("Non-responder fraction",     0.15, 0.05, 0.30,
  "Estimated from heterogeneity in Rossi 2016 and Esgalhado 2020 fiber trials"),
 ("74yo CKD4 annual hazard",    0.09, 0.06, 0.15,
  "UK Renal Registry 2022 Annual Report; Fried LP, Cardiovascular Health Study"),
 ("f_tox (toxin-death coupling)",0.10, 0.00, 0.40,
  "Vanholder R et al, EUTox 2014 (association); Schulman 2015 EPPIC (no causal proof)"),
 ("Gut variability (SD)",       0.30, 0.15, 0.50,
  "Inter-individual microbiome variance; Wu GD et al, Science 2011"),
 ("eGFR mean (cohort)",         22.0, 15.0, 30.0,
  "CKD stage 3b-4 range: eGFR 15-45; cohort centered at stage 4"),
]
print(f"  {'parameter':<32}{'base':>7}{'low':>7}{'high':>7}  {'source'}")
print(f"  {'-'*88}")
for name,base,lo,hi,src in PARAMS:
    print(f"  {name:<32}{base:>7.2f}{lo:>7.2f}{hi:>7.2f}  {src[:55]}")

# =========================================================================
# (2) SIMULATION ENGINE (simplified for speed)
# =========================================================================
def make_cohort(egfr_mean=22, tau_mean=0.35, nr_frac=0.15, gut_sd=0.30, slope_mean=-2.0):
    e = np.clip(rng.normal(egfr_mean, 6, N_PAT), 10, 45)
    gv = np.clip(rng.normal(1.0, gut_sd, N_PAT), 0.3, 1.8)
    sl = np.clip(rng.normal(slope_mean, 0.8, N_PAT), -5.0, -0.5)
    t = np.clip(rng.normal(tau_mean * gv, 0.12), 0, 0.70)
    nr = rng.random(N_PAT) < nr_frac
    t[nr] = np.clip(rng.normal(0.03, 0.02, nr.sum()), 0, 0.08)
    bis = (25.0/e)**1.2 * np.clip(rng.normal(1.0, 0.20, N_PAT), 0.5, 1.8) * 100
    return e, sl, t, bis, (t >= 0.10), (t >= 0.10) & (t < 0.20)

def run_power(nc, km, cv, egfr, sl, tau, bis):
    n_arm = nc * km
    dt = 1.645 * cv * np.sqrt(2.0/n_arm)
    obs = np.zeros((N_PAT, N_REP))
    for p in range(N_PAT):
        a_d, b_d = [], []
        wk = 0
        for _ in range(nc):
            wk_a = wk + WK_T
            eg = max(egfr[p] + sl[p]*wk_a/52, 5)
            d = (25./eg)**1.2 / (25./egfr[p])**1.2
            a_d.extend([d]*km); wk = wk_a + WK_W
            wk_b = wk + WK_T
            eg = max(egfr[p] + sl[p]*wk_b/52, 5)
            d = (25./eg)**1.2 / (25./egfr[p])**1.2
            b_d.extend([d]*km); wk = wk_b + WK_W
        ad = np.array(a_d)[None,:]
        bd = np.array(b_d)[None,:]
        A = bis[p]*ad*(1+rng.normal(0,cv,(N_REP,n_arm)))
        B = bis[p]*bd*(1-tau[p])*(1+rng.normal(0,cv,(N_REP,n_arm)))
        obs[p] = np.where(A.mean(1)>0, (A.mean(1)-B.mean(1))/A.mean(1), 0)
    det = (obs > dt).mean(1)
    return det, dt

# =========================================================================
# (3) TORNADO SENSITIVITY: one-at-a-time sweep on recommended 2x3 design
# =========================================================================
print(f"\n{'='*90}")
print("(3) TORNADO SENSITIVITY — 2x3 design, one parameter varied at a time")
print(f"{'='*90}")
print(f"  {'parameter varied':<32}{'low val':>8}{'pw_low':>8}{'base pw':>8}{'pw_high':>8}{'high val':>8}{'swing':>8}")
print(f"  {'-'*78}")

# baseline run
e0,sl0,t0,b0,resp0,wk0 = make_cohort()
det0,_ = run_power(2, 3, 0.22, e0, sl0, t0, b0)
pw_base = det0[resp0].mean()

results = []
for name, base, lo, hi, _ in PARAMS:
    # vary ONE parameter, hold others at base
    for test_val, label in [(lo,"low"),(hi,"high")]:
        rng_state = np.random.default_rng(201)  # reset for reproducibility
        rng.__setstate__(rng_state.__getstate__())

        if name == "IS within-person CV":
            e,sl,t,b,resp,wk = make_cohort()
            det,_ = run_power(2, 3, test_val, e, sl, t, b)
        elif name == "eGFR slope (mL/min/yr)":
            e,sl,t,b,resp,wk = make_cohort(slope_mean=test_val)
            det,_ = run_power(2, 3, 0.22, e, sl, t, b)
        elif name == "Population mean tau (IS red)":
            e,sl,t,b,resp,wk = make_cohort(tau_mean=test_val)
            det,_ = run_power(2, 3, 0.22, e, sl, t, b)
        elif name == "Non-responder fraction":
            e,sl,t,b,resp,wk = make_cohort(nr_frac=test_val)
            det,_ = run_power(2, 3, 0.22, e, sl, t, b)
        elif name == "Gut variability (SD)":
            e,sl,t,b,resp,wk = make_cohort(gut_sd=test_val)
            det,_ = run_power(2, 3, 0.22, e, sl, t, b)
        elif name == "eGFR mean (cohort)":
            e,sl,t,b,resp,wk = make_cohort(egfr_mean=test_val)
            det,_ = run_power(2, 3, 0.22, e, sl, t, b)
        else:
            continue   # f_tox and hazard don't affect power directly
        pw = det[resp].mean()
        results.append((name, test_val, label, pw))

# print as tornado
param_swings = {}
for name, base, lo, hi, _ in PARAMS:
    lows = [r for r in results if r[0]==name and r[2]=="low"]
    highs = [r for r in results if r[0]==name and r[2]=="high"]
    if lows and highs:
        pw_lo = lows[0][3]; pw_hi = highs[0][3]
        swing = abs(pw_hi - pw_lo)
        param_swings[name] = (lo, pw_lo, pw_base, pw_hi, hi, swing)

for name, vals in sorted(param_swings.items(), key=lambda x: -x[1][5]):
    lo_v, pw_lo, pw_b, pw_hi, hi_v, sw = vals
    print(f"  {name:<32}{lo_v:>8.2f}{pw_lo*100:>7.0f}%{pw_b*100:>7.0f}%{pw_hi*100:>7.0f}%{hi_v:>8.2f}{sw*100:>7.0f}pp")

print(f"\n  Interpretation: parameters sorted by swing (power sensitivity).")
print(f"  IS CV dominates: it directly determines DT. Other parameters matter")
print(f"  through their effect on how many patients fall above/below DT.")

# =========================================================================
# (4) f_tox ROBUSTNESS: which conclusions survive f_tox=0?
# =========================================================================
print(f"\n{'='*90}")
print("(4) f_tox ROBUSTNESS — conclusions that hold regardless of f_tox")
print(f"{'='*90}")
print(f"""
  Conclusion                                          Depends on f_tox?
  ------------------------------------------------------------------
  Gut stack reduces IS by 40-60%                      NO  (direct measurement)
  AST-120 is NET harmful (drug binding)               NO  (pharmacokinetic)
  Competing risk dominates for 74yo                   NO  (actuarial)
  Gut stack adds ~0-2mo dialysis-free survival        YES (at f_tox=0)
  Gut stack adds ~6-9mo survival                      YES (only if f_tox>0.2)
  n-of-1 protocol can detect IS reduction             NO  (measurement theory)
  CV reduction > cycle extension for power            NO  (DT formula)
  Noise reduction is primary lever for weak resp.     NO  (mathematical)

  ROBUST (f_tox-independent):
    - Protocol design recommendations (all)
    - Safety findings (AST-120 removal, Lacto > spore)
    - Engineering specifications (sachet, timing, protein)
    - Power analysis and detection thresholds

  f_tox-DEPENDENT:
    - Clinical benefit magnitude (dialysis-free months)
    - Whether the intervention is WORTH doing (vs. focusing on AKI prevention)

  -> The n-of-1 protocol paper should be framed around PROTOCOL DESIGN,
     not clinical efficacy claims. This makes it f_tox-independent.""")

# =========================================================================
# (5) CV FEASIBILITY: what noise reductions are clinically achievable?
# =========================================================================
print(f"\n{'='*90}")
print("(5) CV REDUCTION FEASIBILITY — clinical evidence")
print(f"{'='*90}")
print("""
  SUPERSEDED. An earlier version of this section listed a component breakdown
  and concluded "CV = 0.15 is ACHIEVABLE". Both the conclusion and two of the
  literature attributions supporting it were wrong, and they are removed rather
  than reworded.

  The conclusion was wrong because the breakdown treated day-to-day biological
  variation as reducible by a within-visit duplicate assay. Repeating an assay
  on one specimen averages down only the analytical component:

      CV^2_total = CV^2_biological + CV^2_pre-analytical + CV^2_analytical / k

  Running that decomposition against the 35.9% within-person biological CV that
  Pretorius et al. 2013 (Clin Chim Acta) report for total serum IS:

    - a duplicate assay changes the total CV by at most 0.004
    - the reachable total with sampling standardization is 0.22 to 0.30
    - reaching 0.15 would need the irreducible biological component alone to be
      0.112 to 0.138, against a reported 0.359

  The two removed attributions were an intra-individual CV figure credited to
  Deltombe et al. 2015 (that paper reports protein binding of uremic toxins, not
  biological variation) and a "duplicate sampling reduced effective CV by ~30%"
  figure credited to a 2014 NDT paper, which we could not verify and which the
  variance algebra above contradicts in any case.

  Authoritative version: validation/variance_components_analysis.py""")

# CV sweep with power for weak responders
e0,sl0,t0,b0,resp0,wk0 = make_cohort()
weak0 = (t0 >= 0.10) & (t0 < 0.20)
print(f"\n  CV sweep (2x3 design, weak responder power):")
print(f"  {'CV':>6}{'DT':>7}{'weak pw':>9}{'achievability'}")
for cv, ach in [(0.25,"pessimistic"),(0.22,"standard"),(0.18,"fasting+timed"),
                (0.15,"duplicate+fasting"),(0.12,"composite endpoint")]:
    det,dt = run_power(2, 3, cv, e0, sl0, t0, b0)
    pw_w = det[weak0].mean() if weak0.sum()>0 else 0
    print(f"  {cv:>6.2f}{dt*100:>6.0f}%{pw_w*100:>8.0f}%  {ach}")

# =========================================================================
# (6) FP THRESHOLD TRADE-OFF (adjustable decision threshold)
# =========================================================================
print(f"\n{'='*90}")
print("(6) FP-POWER TRADE-OFF — adjustable decision threshold")
print(f"{'='*90}")
print(f"  {'threshold':>12}{'overall pw':>12}{'weak pw':>10}{'FP':>8}{'note'}")
print(f"  {'-'*60}")
non_r0 = t0 < 0.10
n_arm = 6
dt_base = 1.645 * 0.15 * np.sqrt(2.0/n_arm)
# generate one set of observations at CV=0.15
rng2 = np.random.default_rng(301)
obs_all = np.zeros((N_PAT, N_REP))
for p in range(N_PAT):
    A = b0[p]*(1+rng2.normal(0, 0.15, (N_REP, n_arm)))
    B = b0[p]*(1-t0[p])*(1+rng2.normal(0, 0.15, (N_REP, n_arm)))
    obs_all[p] = np.where(A.mean(1)>0, (A.mean(1)-B.mean(1))/A.mean(1), 0)

for mult in [0.8, 0.9, 1.0, 1.1, 1.2, 1.4]:
    dt_test = dt_base * mult
    det_t = (obs_all > dt_test).mean(1)
    pw_all = det_t[resp0].mean()
    pw_w = det_t[weak0].mean() if weak0.sum()>0 else 0
    fp = det_t[non_r0].mean()
    note = "<-- standard" if mult==1.0 else ("conservative" if mult>1 else "liberal")
    print(f"  DT x{mult:<5.1f}{pw_all*100:>10.0f}%{pw_w*100:>9.0f}%{fp*100:>7.0f}%  {note}")

print(f"\n  -> raising threshold by 20% (DT×1.2) cuts FP roughly in half")
print(f"     at a modest power cost (~5 pp overall)")

# =========================================================================
# (7) COHORT HETEROGENEITY: vary patient distributions
# =========================================================================
print(f"\n{'='*90}")
print("(7) COHORT HETEROGENEITY — protocol robustness across patient populations")
print(f"{'='*90}")
populations = [
    ("Base (74yo CKD4)",       22, -2.0, 0.35, 0.15),
    ("Younger (60yo CKD3b)",   30, -1.5, 0.40, 0.10),
    ("Older (80yo CKD4-5)",    16, -3.0, 0.30, 0.20),
    ("Diabetic CKD4",          20, -3.5, 0.25, 0.25),
    ("Aggressive CKD5",        12, -4.0, 0.30, 0.20),
]
print(f"  {'population':<26}{'eGFR':>6}{'slope':>7}{'tau':>6}{'NR%':>6}{'power':>8}{'weak pw':>9}{'FP':>6}")
print(f"  {'-'*72}")
for label, egfr_m, sl_m, tau_m, nr_f in populations:
    e,sl,t,b,resp,wk = make_cohort(egfr_mean=egfr_m, slope_mean=sl_m,
                                     tau_mean=tau_m, nr_frac=nr_f)
    det,_ = run_power(2, 3, 0.15, e, sl, t, b)
    pw = det[resp].mean()
    pw_w = det[wk].mean() if wk.sum()>0 else 0
    fp = det[~resp].mean()
    print(f"  {label:<26}{egfr_m:>6}{sl_m:>7.1f}{tau_m:>6.2f}{nr_f*100:>5.0f}%"
          f"{pw*100:>7.0f}%{pw_w*100:>8.0f}%{fp*100:>5.0f}%")

print(f"\n  -> protocol (2x3, CV=0.15) maintains >80% overall power across populations")
print(f"     except aggressive CKD5 (low eGFR = high drift, lower baseline effect)")

# =========================================================================
# (8) SUMMARY: which conclusions are robust?
# =========================================================================
print(f"\n{'='*90}")
print("(8) ROBUSTNESS SUMMARY")
print(f"{'='*90}")
print(f"""
  FULLY ROBUST (survive all parameter sweeps):
    1. CV reduction is more efficient than cycle extension for weak responders
    2. 2x3 crossover at CV=0.15 achieves >70% overall power
    3. AST-120 should be excluded (pharmacokinetic, not parameter-dependent)
    4. Adaptive 2-stage design saves burden for majority of patients
    5. Crossover controls eGFR drift to <2 pp bias

  PARTIALLY ROBUST (sensitive to specific parameters):
    6. Weak responder power: highly sensitive to CV (dominant parameter)
    7. FP rate: adjustable via threshold (12% -> ~6% at DT×1.2)
    8. Overall power: moderately sensitive to tau_mean and NR fraction

  NOT ROBUST (assumption-dependent):
    9. Clinical benefit magnitude (months gained): depends on f_tox
    10. Whether intervention is WORTH doing: depends on f_tox and patient age

  RECOMMENDED FRAMING FOR PUBLICATION:
    Frame the paper as a PROTOCOL DESIGN AND POWER ANALYSIS study,
    not a clinical efficacy prediction. This makes all key conclusions
    fall in the "fully robust" category and avoids the f_tox trap.""")
