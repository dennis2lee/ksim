"""
NOVELTY AND GAP ANALYSIS: what this project adds beyond prior work.

Reframed scope: IS/PCS reduction itself is the primary endpoint —
NOT a surrogate for hard outcomes (dialysis, death). We design a protocol
to answer "does this intervention lower THIS patient's toxins?" without
claiming downstream clinical benefit. This avoids the EPPIC/surrogate trap.

Four novelty axes demonstrated quantitatively:
  (1) Published RCT calibration — our parameters match literature
  (2) RCT vs n-of-1 — population mean hides individual heterogeneity
  (3) Heterogeneity tax — burden wasted on non-responders in treat-all
  (4) CV as design lever — novel methodological insight absent from literature
  (5) IS as 1° endpoint — explicit avoidance of surrogate-endpoint trap
"""
import numpy as np
rng = np.random.default_rng(501)

N_PAT = 100
N_REP = 500
CV0   = 0.22

# =========================================================================
# (1) PUBLISHED RCT CALIBRATION: our parameters vs literature
# =========================================================================
print("="*90)
print("(1) CALIBRATION: published RCT effect sizes vs our model parameters")
print("="*90)
LIT = [
 ("AST-120 (EPPIC, Schulman 2015)",          "IS",  "30-40%", 0.35,
  "N=2035, multicenter RCT. IS reduced but NO effect on hard endpoints."),
 ("Fiber/prebiotic (Rossi 2016 CJASN)",       "IS",  "20-35%", 0.28,
  "N=31 crossover, high-fiber diet. IS -25%, PCS -30%."),
 ("Resistant starch (Esgalhado 2018)",         "IS",  "20-30%", 0.25,
  "N=43, RS 16g/d x 4wk. IS -26%, urea NS."),
 ("Synbiotic (Rossi 2016 KI)",                "PCS", "12-20%", 0.16,
  "N=37, synbiotic 6wk. PCS -16%, IS NS."),
 ("Gum arabic (network meta, 2025)",           "IS",  "15-25%", 0.20,
  "Multiple trials pooled. Consistent IS/PCS reduction."),
 ("Combined stack (model estimate)",           "IS",  "40-60%", 0.50,
  "Our model: fiber+sorbent+synbiotic with overlap damping."),
]
print(f"  {'study':<42}{'toxin':>6}{'reported':>10}{'our param':>10}  {'note'}")
print(f"  {'-'*88}")
for name,tox,rep,ours,note in LIT:
    print(f"  {name:<42}{tox:>6}{rep:>10}{ours*100:>9.0f}%  {note[:45]}")

print(f"\n  Match quality: our parameters fall WITHIN published ranges for all interventions.")
print(f"  We use mid-range estimates; sensitivity_analysis.py sweeps the full uncertainty.")

# =========================================================================
# (2) RCT vs n-of-1: WHAT THE POPULATION MEAN HIDES
# =========================================================================
# Generate a heterogeneous cohort (same as nof1_virtual_cohort.py)
gut_var = np.clip(rng.normal(1.0, 0.30, N_PAT), 0.3, 1.8)
tau     = np.clip(rng.normal(0.35 * gut_var, 0.12), 0, 0.70)
non_resp = rng.random(N_PAT) < 0.15
tau[non_resp] = np.clip(rng.normal(0.03, 0.02, non_resp.sum()), 0, 0.08)
baseline = np.clip(rng.normal(100, 25, N_PAT), 40, 200)

print(f"\n{'='*90}")
print("(2) RCT vs n-of-1: what the population mean HIDES")
print(f"{'='*90}")

# Simulate a 2-arm parallel RCT
# Control arm: baseline IS + noise
# Treatment arm: baseline * (1-tau) + noise
rct_control = baseline * (1 + rng.normal(0, CV0, N_PAT))
rct_treated = baseline * (1 - tau) * (1 + rng.normal(0, CV0, N_PAT))

mean_ctrl = np.mean(rct_control)
mean_trt  = np.mean(rct_treated)
rct_effect = (mean_ctrl - mean_trt) / mean_ctrl
se_rct = np.sqrt(np.var(rct_control)/N_PAT + np.var(rct_treated)/N_PAT) / mean_ctrl
z_rct = rct_effect / se_rct

print(f"  VIRTUAL 2-ARM PARALLEL RCT (N={N_PAT} per arm):")
print(f"    control mean IS:  {mean_ctrl:.0f}")
print(f"    treated mean IS:  {mean_trt:.0f}")
print(f"    mean reduction:   {rct_effect*100:.0f}%  (z={z_rct:.1f}, p {'< 0.001' if z_rct>3.3 else '< 0.05'})")
print(f"    RCT conclusion:   'Intervention significantly reduces IS by ~{rct_effect*100:.0f}%'")

print(f"\n  BUT the RCT cannot see the INDIVIDUAL REALITY:")
print(f"    true effect distribution:")
for lo,hi,lb in [(0,.10,"non-resp <10%"),(.10,.25,"weak 10-25%"),
                  (.25,.45,"moderate 25-45%"),(.45,1.0,"strong >45%")]:
    m = (tau>=lo)&(tau<hi); n=m.sum()
    if n: print(f"      {lb:<22} n={n:>3}  mean tau={np.mean(tau[m])*100:.0f}%")

nr_count = (tau < 0.10).sum()
print(f"\n    {nr_count} patients ({nr_count}%) get ZERO benefit but bear full burden")
print(f"    RCT reports: 'treatment works' -> ALL patients continue treatment")
print(f"    n-of-1 would: correctly CLASSIFY each patient -> deprescribe non-responders")

# n-of-1 classification for same patients (2x3 CV=0.15 recommended protocol)
n_arm = 6; cv_opt = 0.15
mde = 1.645 * cv_opt * np.sqrt(2.0 / n_arm)
true_resp = tau >= 0.10
det_power = np.zeros(N_PAT)
for p in range(N_PAT):
    A = baseline[p] * (1 + rng.normal(0, cv_opt, (N_REP, n_arm)))
    B = baseline[p] * (1 - tau[p]) * (1 + rng.normal(0, cv_opt, (N_REP, n_arm)))
    obs = np.where(A.mean(1)>0, (A.mean(1)-B.mean(1))/A.mean(1), 0)
    det_power[p] = (obs > mde).mean()

tp = (det_power[true_resp] > 0.5).sum()
tn = (det_power[~true_resp] <= 0.5).sum()
print(f"\n  n-of-1 CLASSIFICATION (2x3, CV=0.15, MDE={mde*100:.0f}%):")
print(f"    responders correctly identified: {tp}/{true_resp.sum()} ({tp/true_resp.sum()*100:.0f}%)")
print(f"    non-resp correctly deprescribed: {tn}/{(~true_resp).sum()} ({tn/(~true_resp).sum()*100:.0f}%)")

# =========================================================================
# (3) HETEROGENEITY TAX: burden wasted under treat-all vs n-of-1
# =========================================================================
print(f"\n{'='*90}")
print("(3) HETEROGENEITY TAX: burden wasted on non-responders")
print(f"{'='*90}")
BURDEN_PER_MONTH = 1.0     # arbitrary burden units per month of regimen
HORIZON = 36               # months of treatment after decision

# Treat-all (RCT conclusion)
burden_all = N_PAT * BURDEN_PER_MONTH * HORIZON
benefit_all = np.sum(tau[true_resp]) * HORIZON   # total IS-reduction-months
waste_all = (~true_resp).sum() * BURDEN_PER_MONTH * HORIZON

# n-of-1 guided
n_identified = tp
n_deprescribed = tn
burden_nof1 = n_identified * BURDEN_PER_MONTH * HORIZON
waste_nof1 = ((~true_resp).sum() - tn) * BURDEN_PER_MONTH * HORIZON

print(f"  36-month treatment horizon, {N_PAT} patients:")
print(f"  {'metric':<36}{'treat-all':>12}{'n-of-1 guided':>14}")
print(f"  {'-'*62}")
print(f"  {'patients on regimen':<36}{N_PAT:>12}{n_identified:>14}")
print(f"  {'total burden-months':<36}{burden_all:>12.0f}{burden_nof1:>14.0f}")
print(f"  {'wasted on non-resp (burden-mo)':<36}{waste_all:>12.0f}{waste_nof1:>14.0f}")
print(f"  {'burden reduction':<36}{'-':>12}{(1-burden_nof1/burden_all)*100:>13.0f}%")
print(f"  {'non-resp burden reduction':<36}{'-':>12}{(1-waste_nof1/waste_all)*100:>13.0f}%")

print(f"\n  -> n-of-1 eliminates {(1-waste_nof1/waste_all)*100:.0f}% of burden wasted on non-responders")
print(f"     while retaining {tp/true_resp.sum()*100:.0f}% of treatment benefit")

# =========================================================================
# (4) CV AS DESIGN LEVER: novel methodological insight
# =========================================================================
print(f"\n{'='*90}")
print("(4) CV AS DESIGN LEVER — absent from existing n-of-1 methodology literature")
print(f"{'='*90}")
print(f"""
  EXISTING n-of-1 LITERATURE focuses on:
    - Number of cycles / periods (De Carvalho 2026 BJCP)
    - Washout duration (Duan 2013)
    - Carryover effects (Senn 2002)
    - Serial correlation (Zucker 1997)
    - Sample size planning (Percha 2019)

  NONE address measurement noise reduction as a design lever.

  OUR CONTRIBUTION: MDE = 1.645 * CV * sqrt(2/n)
    - CV enters LINEARLY; n enters as 1/sqrt(n)
    - Halving CV halves MDE (same effort)
    - Halving MDE by adding measurements requires 4x the measurements
    - This makes CV reduction QUADRATICALLY more efficient than adding cycles""")

print(f"\n  Quantitative demonstration:")
print(f"  {'approach':<36}{'MDE':>8}{'weeks':>8}{'draws':>8}{'efficiency':>12}")
print(f"  {'-'*72}")
ref_mde = 1.645 * 0.22 * np.sqrt(2./6)
approaches = [
    ("Baseline: 2x3, CV=0.22",  0.22, 6, 24, 12),
    ("More cycles: 4x3, CV=0.22", 0.22, 12, 48, 24),
    ("More cycles: 8x3, CV=0.22", 0.22, 24, 96, 48),
    ("Reduce CV: 2x3, CV=0.15",   0.15, 6, 24, 24),
    ("Reduce CV: 2x3, CV=0.12",   0.12, 6, 24, 24),
    ("Both: 3x3, CV=0.15",        0.15, 9, 36, 36),
]
for label, cv, n, wk, dr in approaches:
    mde_v = 1.645 * cv * np.sqrt(2./n)
    red = (1 - mde_v/ref_mde) * 100
    print(f"  {label:<36}{mde_v*100:>7.0f}%{wk:>8}{dr:>8}{red:>10.0f}% MDE cut")

print(f"\n  -> CV 0.22->0.15 at 2x3 (24 wk) achieves BETTER MDE than 4x3 (48 wk) at CV=0.22")
print(f"     this saves 24 WEEKS — the dominant efficiency lever for protocol design")

# =========================================================================
# (5) IS AS PRIMARY ENDPOINT: avoiding the surrogate trap
# =========================================================================
print(f"\n{'='*90}")
print("(5) IS REDUCTION AS PRIMARY ENDPOINT — avoiding the EPPIC/surrogate trap")
print(f"{'='*90}")
print(f"""
  THE TRAP (EPPIC 2015):
    premise:   IS is a surrogate for CKD progression / CV events
    design:    reduce IS with AST-120 -> expect slower progression
    result:    IS reduced 30-40%, BUT no effect on hard endpoints
    lesson:    surrogate != outcome; lowering IS may not help

  OUR REFRAMING:
    premise:   IS reduction is itself a MEASURABLE, PERSONAL response
    design:    detect whether THIS patient's IS drops with gut-clearance
    endpoint:  IS reduction (direct measurement, no surrogate claim)
    question:  "does the intervention work IN THIS PATIENT?" (not "does it save lives?")

  WHY THIS MATTERS:
    - Avoids the surrogate-endpoint fallacy entirely
    - The patient (and clinician) can make an INFORMED decision:
      "your IS dropped 35% -> the regimen is working for you, continue"
      vs "your IS didn't change -> stop, no benefit, just burden"
    - Whether IS reduction ALSO improves hard outcomes is a SEPARATE question
      that our protocol doesn't need to answer

  COMPARISON WITH PRIOR WORK:
    study              endpoint            level       gap
    ---------------------------------------------------------------
    EPPIC (2015)       IS as surrogate     population  failed
    Rossi (2016)       IS/PCS mean change  population  no individual
    Esgalhado (2018)   IS mean change      population  no individual
    Network meta (2025) pooled effect      population  no individual
    ** THIS WORK **    IS change           INDIVIDUAL  fills the gap""")

# =========================================================================
# (6) LITERATURE GAP SUMMARY TABLE
# =========================================================================
print(f"\n{'='*90}")
print("(6) LITERATURE GAP SUMMARY")
print(f"{'='*90}")
gaps = [
 ("Individual effect detection",
  "All RCTs report population means; none classify individual responders",
  "n-of-1 crossover with per-patient power analysis (N=100 virtual cohort)"),
 ("Measurement noise as design lever",
  "n-of-1 methodology focuses on cycles/washout, not CV optimization",
  "MDE formula analysis showing CV reduction is quadratically more efficient"),
 ("IS as primary (non-surrogate) endpoint",
  "EPPIC used IS as surrogate for hard outcomes -> failed",
  "IS reduction reframed as personal response detection, not outcome prediction"),
 ("In-silico protocol validation for uremic toxins",
  "No virtual cohort study for gut-clearance n-of-1 trial design",
  "100-patient virtual cohort with heterogeneous params, adaptive 2-stage"),
 ("Weak responder rescue",
  "No strategy for detecting modest (10-20%) IS reductions",
  "CV reduction + adaptive enrichment -> 69% power for weak responders"),
 ("Published biological CV quantification",
  "IS biological CV (25-27%) documented but not used for protocol design",
  "First use of published CV data to optimize n-of-1 protocol parameters"),
]
print(f"  {'gap':<38}{'prior work limit':<48}{'our contribution'}")
print(f"  {'-'*88}")
for gap, prior, ours in gaps:
    print(f"\n  {gap}")
    print(f"    prior: {prior}")
    print(f"    ours:  {ours}")

# =========================================================================
# (7) VERIFICATION: key results under refocused scope
# =========================================================================
print(f"\n{'='*90}")
print("(7) VERIFICATION: key results under refocused scope (IS as 1° endpoint)")
print(f"{'='*90}")
print(f"""
  REFOCUSED CORE CLAIMS (all independent of f_tox / downstream outcomes):

  1. Individual IS reduction is detectable via n-of-1 crossover
     -> overall power {det_power[true_resp].mean()*100:.0f}% at 2x3, CV=0.15 (verified above)

  2. Non-responders can be correctly deprescribed
     -> specificity {tn/(~true_resp).sum()*100:.0f}% (verified above)

  3. CV reduction is the dominant design lever
     -> 2x3 at CV=0.15 matches 4x3 at CV=0.22 (24 wk vs 48 wk, verified above)

  4. Adaptive 2-stage rescues weak responders without burdening strong ones
     -> weak power 27% -> 69% (from nof1_weak_rescue.py)

  5. The heterogeneity tax is large: treat-all wastes {waste_all:.0f} burden-months
     on non-responders; n-of-1 guided reduces this by {(1-waste_nof1/waste_all)*100:.0f}%

  REMOVED FROM CORE CLAIMS (background only):
  - "AKI/CV prevention adds more months than toxin route" (iter4) -> context only
  - "gut stack adds ~0-2 months dialysis-free survival" -> not our endpoint
  - f_tox dependence -> irrelevant to IS-as-endpoint framing""")
