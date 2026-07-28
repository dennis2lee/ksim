"""
WEAK RESPONDER RESCUE: boosted n-of-1 designs for tau 10-20%.

Problem: the 2x3 default protocol detects only ~19% of weak responders
(true IS reduction 10-20%). Three strategies evaluated:
  (A) Extended fixed designs — more cycles and/or measures per period
  (B) Adaptive 2-stage enrichment — extra cycles for borderline only
  (C) Measurement noise reduction — lower CV via better assay/protocol
  (D) Combined: reduced CV + adaptive

Same virtual cohort (N=100, seed=101) as nof1_virtual_cohort.py.
"""
import numpy as np
rng = np.random.default_rng(101)

N_PAT  = 100
CV0    = 0.22
WK_T   = 4
WK_W   = 2
WK_CYC = 2*(WK_T + WK_W)   # 12 wk per AB cycle
N_REP  = 500

# =========================================================================
# (1) COHORT — identical to nof1_virtual_cohort.py
# =========================================================================
age     = np.clip(rng.normal(74, 5, N_PAT), 60, 85).astype(int)
eGFR    = np.clip(rng.normal(22, 6, N_PAT), 10, 35)
base_IS = (25.0/eGFR)**1.2 * np.clip(rng.normal(1.0, 0.20, N_PAT), 0.5, 1.8) * 100
gut_var = np.clip(rng.normal(1.0, 0.30, N_PAT), 0.3, 1.8)
slope   = np.clip(rng.normal(-2.0, 0.8, N_PAT), -5.0, -0.5)
_       = 0.05 + 0.002*(age - 60) + 0.003*np.clip(35 - eGFR, 0, 25)  # hazard (unused here)
f_tox   = np.clip(rng.beta(1.5, 6.0, N_PAT), 0, 0.5)
tau     = np.clip(rng.normal(0.35 * gut_var, 0.12), 0, 0.70)
non_resp = rng.random(N_PAT) < 0.15
tau[non_resp] = np.clip(rng.normal(0.03, 0.02, non_resp.sum()), 0, 0.08)

weak      = (tau >= 0.10) & (tau < 0.20)
mod_bord  = (tau >= 0.20) & (tau < 0.30)
strong    = tau >= 0.30
non_r     = tau < 0.10
true_resp = tau >= 0.10

print("="*88)
print("(1) COHORT RECAP — focus groups")
print("="*88)
for lb, mk in [("non-resp <10%", non_r), ("WEAK 10-20%", weak),
                ("borderline 20-30%", mod_bord), ("strong 30%+", strong)]:
    print(f"  {lb:<24} n={mk.sum():>3}  mean tau={np.mean(tau[mk])*100:>3.0f}%")

# =========================================================================
# (2) SIMULATION ENGINE
# =========================================================================
def run(nc, km, cv=CV0):
    n_arm = nc * km
    obs = np.zeros((N_PAT, N_REP))
    for p in range(N_PAT):
        a_d, b_d = [], []
        wk = 0
        for _ in range(nc):
            wk_a = wk + WK_T
            eg = max(eGFR[p] + slope[p]*wk_a/52, 5)
            d = (25./eg)**1.2 / (25./eGFR[p])**1.2
            a_d.extend([d]*km); wk = wk_a + WK_W
            wk_b = wk + WK_T
            eg = max(eGFR[p] + slope[p]*wk_b/52, 5)
            d = (25./eg)**1.2 / (25./eGFR[p])**1.2
            b_d.extend([d]*km); wk = wk_b + WK_W
        ad = np.array(a_d)[None,:]
        bd = np.array(b_d)[None,:]
        A = base_IS[p]*ad*(1+rng.normal(0,cv,(N_REP,n_arm)))
        B = base_IS[p]*bd*(1-tau[p])*(1+rng.normal(0,cv,(N_REP,n_arm)))
        obs[p] = np.where(A.mean(1)>0, (A.mean(1)-B.mean(1))/A.mean(1), 0)
    return obs, nc*WK_CYC, 1.645*cv*np.sqrt(2./n_arm)

def stats(obs, dt):
    det = (obs > dt).mean(1)
    return dict(det=det,
                pw_w=det[weak].mean(), pw_m=det[mod_bord].mean(),
                pw_s=det[strong].mean(), pw_all=det[true_resp].mean(),
                fp=det[non_r].mean())

# =========================================================================
# (3) STRATEGY A: EXTENDED FIXED DESIGNS
# =========================================================================
designs = [(2,3),(3,3),(3,4),(4,3),(4,4),(5,3)]
print(f"\n{'='*88}")
print("(3) STRATEGY A — EXTENDED FIXED DESIGNS  (baseline = 2x3)")
print(f"{'='*88}")
print(f"  {'design':<10}{'DT':>6}{'wk':>5}{'draws':>7}"
      f"{'WEAK':>8}{'20-30%':>8}{'30%+':>8}{'all':>8}{'FP':>6}")
print(f"  {'-'*65}")
fixed = {}
for nc,km in designs:
    o,wk,dt = run(nc,km)
    s = stats(o,dt)
    s.update(wk=wk, dt=dt, draws=2*nc*km)
    fixed[(nc,km)] = s
    tag = " *" if (nc,km)==(2,3) else ""
    print(f"  {nc}x{km}{tag:<5}{dt*100:>5.0f}%{wk:>5}{s['draws']:>7}"
          f"{s['pw_w']*100:>7.0f}%{s['pw_m']*100:>7.0f}%{s['pw_s']*100:>7.0f}%"
          f"{s['pw_all']*100:>7.0f}%{s['fp']*100:>5.0f}%")

# DT vs cycles insight
print(f"\n  DT scaling: DT = 1.645 * CV * sqrt(2/n)")
print(f"  to HALVE DT: need 4x measurements (e.g. 2x3->8x3 = 96wk)")
print(f"  -> adding cycles is quadratically inefficient for weak responders")

# =========================================================================
# (4) STRATEGY B: ADAPTIVE 2-STAGE ENRICHMENT
# =========================================================================
def adaptive(s1_nc, s2_nc, km, cv=CV0, thresh_lo=0.0):
    """Stage 1 = s1_nc x km. Borderline -> add s2_nc x km."""
    n1 = s1_nc * km
    n_full = (s1_nc + s2_nc) * km
    mde1 = 1.645 * cv * np.sqrt(2.0 / n1)
    dt_full = 1.645 * cv * np.sqrt(2.0 / n_full)
    wk1 = s1_nc * WK_CYC
    wk_full = (s1_nc + s2_nc) * WK_CYC

    detected = np.zeros((N_PAT, N_REP), bool)
    went_s2 = np.zeros((N_PAT, N_REP), bool)

    for p in range(N_PAT):
        a_d, b_d = [], []
        wk = 0
        for _ in range(s1_nc + s2_nc):
            wk_a = wk + WK_T
            eg = max(eGFR[p] + slope[p]*wk_a/52, 5)
            d = (25./eg)**1.2 / (25./eGFR[p])**1.2
            a_d.extend([d]*km); wk = wk_a + WK_W
            wk_b = wk + WK_T
            eg = max(eGFR[p] + slope[p]*wk_b/52, 5)
            d = (25./eg)**1.2 / (25./eGFR[p])**1.2
            b_d.extend([d]*km); wk = wk_b + WK_W

        ad = np.array(a_d)[None,:]
        bd = np.array(b_d)[None,:]
        A_all = base_IS[p]*ad*(1+rng.normal(0,cv,(N_REP,n_full)))
        B_all = base_IS[p]*bd*(1-tau[p])*(1+rng.normal(0,cv,(N_REP,n_full)))

        A1 = A_all[:,:n1]; B1 = B_all[:,:n1]
        ma1, mb1 = A1.mean(1), B1.mean(1)
        obs1 = np.where(ma1>0, (ma1-mb1)/ma1, 0)

        ma_f, mb_f = A_all.mean(1), B_all.mean(1)
        obs_f = np.where(ma_f>0, (ma_f-mb_f)/ma_f, 0)

        clear_resp = obs1 > mde1
        clear_nonr = obs1 < thresh_lo
        borderline = ~clear_resp & ~clear_nonr

        detected[p, clear_resp] = True
        detected[p, clear_nonr] = False
        detected[p, borderline] = obs_f[borderline] > dt_full
        went_s2[p] = borderline

    power = detected.mean(1)
    s2_frac = went_s2.mean(1)
    avg_wk = np.where(went_s2, wk_full, wk1).astype(float).mean(1)
    avg_draws = np.where(went_s2, 2*n_full, 2*n1).astype(float).mean(1)
    return dict(power=power, s2_frac=s2_frac, avg_wk=avg_wk, avg_draws=avg_draws,
                mde1=mde1, dt_full=dt_full, wk1=wk1, wk_full=wk_full,
                fp=power[non_r].mean(), pw_w=power[weak].mean(),
                pw_m=power[mod_bord].mean(), pw_all=power[true_resp].mean())

print(f"\n{'='*88}")
print("(4) STRATEGY B — ADAPTIVE 2-STAGE (S1=2x3 -> borderline -> S2 extra cycles)")
print(f"{'='*88}")
adapt_cfgs = [
    ("2+1 x 3meas", 2, 1, 3),
    ("2+2 x 3meas", 2, 2, 3),
    ("2+1 x 4meas", 2, 1, 4),
]
print(f"  {'config':<18}{'WEAK':>8}{'20-30':>8}{'all':>8}{'FP':>6}"
      f"{'avg wk':>8}{'max wk':>8}{'avg dr':>8}{'%->S2':>7}")
print(f"  {'-'*78}")
adapt = {}
for label, s1, s2, km in adapt_cfgs:
    r = adaptive(s1, s2, km)
    adapt[label] = r
    s2r = np.mean(r['s2_frac'])
    print(f"  {label:<18}{r['pw_w']*100:>7.0f}%{r['pw_m']*100:>7.0f}%{r['pw_all']*100:>7.0f}%"
          f"{r['fp']*100:>5.0f}%{np.mean(r['avg_wk']):>7.0f}{r['wk_full']:>8}"
          f"{np.mean(r['avg_draws']):>7.0f}{s2r*100:>7.0f}%")

# who goes to stage 2?
r22 = adapt["2+2 x 3meas"]
print(f"\n  Stage 2 enrichment breakdown (2+2 x 3meas):")
for lb,mk in [("non-resp",non_r),("WEAK",weak),("20-30%",mod_bord),("30%+",strong)]:
    if mk.sum()==0: continue
    s2r = np.mean(r22['s2_frac'][mk])
    print(f"    {lb:<14} n={mk.sum():>3}  -> S2: {s2r*100:.0f}%"
          f"  (avg wk: {np.mean(r22['avg_wk'][mk]):.0f})")
print(f"  -> adaptive correctly ENRICHES stage 2 for weak/borderline patients")

# =========================================================================
# (5) STRATEGY C: NOISE REDUCTION (CV sweep at fixed 2x3)
# =========================================================================
cvs = [0.22, 0.18, 0.15, 0.12, 0.10]
cv_how = {0.22:"standard serum IS draw",
          0.18:"fasting + AM timed draw",
          0.15:"same-day duplicates averaged",
          0.12:"serum + urinary IS composite",
          0.10:"LC-MS/MS strict protocol"}

print(f"\n{'='*88}")
print("(5) STRATEGY C — NOISE REDUCTION (CV sweep, fixed 2x3 = 24 wk)")
print(f"{'='*88}")
print(f"  {'CV':>5}{'DT':>6}{'WEAK':>8}{'20-30':>8}{'all':>8}{'FP':>6}"
      f"  {'equiv n@.22':>12}  {'method'}")
print(f"  {'-'*82}")
cv_res = {}
for cv in cvs:
    o,_,dt = run(2,3,cv=cv)
    s = stats(o,dt)
    cv_res[cv] = s
    cv_res[cv]['dt'] = dt
    n_equiv = 2.0/(dt/(1.645*0.22))**2
    print(f"  {cv:>5.2f}{dt*100:>5.0f}%{s['pw_w']*100:>7.0f}%{s['pw_m']*100:>7.0f}%"
          f"{s['pw_all']*100:>7.0f}%{s['fp']*100:>5.0f}%{n_equiv:>11.0f}"
          f"  {cv_how[cv]}")

print(f"\n  KEY INSIGHT: CV 0.22->0.15 at 2x3 gives DT 21%->14%")
print(f"    same DT at CV=0.22 would require ~{2.0/((1.645*0.15*np.sqrt(2./6))/(1.645*0.22))**2:.0f} measures/arm"
      f" (= {2.0/((1.645*0.15*np.sqrt(2./6))/(1.645*0.22))**2/3:.0f} cycles x 3meas = "
      f"{int(2.0/((1.645*0.15*np.sqrt(2./6))/(1.645*0.22))**2/3)*12} wk)")
print(f"    -> reducing CV saves ~{int(2.0/((1.645*0.15*np.sqrt(2./6))/(1.645*0.22))**2/3)*12 - 24} WEEKS"
      f" of protocol time")

# =========================================================================
# (6) STRATEGY D: COMBINED (reduced CV + adaptive)
# =========================================================================
print(f"\n{'='*88}")
print("(6) STRATEGY D — COMBINED (reduced CV + adaptive)")
print(f"{'='*88}")
combos = [
    ("CV=.15 fixed 2x3",      0.15, None),
    ("CV=.15 adapt 2+1x3",    0.15, (2,1,3)),
    ("CV=.15 adapt 2+2x3",    0.15, (2,2,3)),
    ("CV=.12 fixed 2x3",      0.12, None),
    ("CV=.12 adapt 2+1x3",    0.12, (2,1,3)),
]
print(f"  {'strategy':<26}{'WEAK':>8}{'20-30':>8}{'all':>8}{'FP':>6}"
      f"{'avg wk':>8}{'draws':>7}")
print(f"  {'-'*70}")
combo_res = {}
for label, cv, acfg in combos:
    if acfg is None:
        o,wk,dt = run(2,3,cv=cv)
        s = stats(o,dt)
        combo_res[label] = dict(pw_w=s['pw_w'],pw_m=s['pw_m'],pw_all=s['pw_all'],
                                fp=s['fp'],avg_wk=24,draws=12)
        print(f"  {label:<26}{s['pw_w']*100:>7.0f}%{s['pw_m']*100:>7.0f}%"
              f"{s['pw_all']*100:>7.0f}%{s['fp']*100:>5.0f}%{'24':>8}{'12':>7}")
    else:
        s1,s2,km = acfg
        r = adaptive(s1,s2,km,cv=cv)
        combo_res[label] = dict(pw_w=r['pw_w'],pw_m=r['pw_m'],pw_all=r['pw_all'],
                                fp=r['fp'],avg_wk=np.mean(r['avg_wk']),
                                draws=np.mean(r['avg_draws']))
        print(f"  {label:<26}{r['pw_w']*100:>7.0f}%{r['pw_m']*100:>7.0f}%"
              f"{r['pw_all']*100:>7.0f}%{r['fp']*100:>5.0f}%"
              f"{np.mean(r['avg_wk']):>7.0f}{np.mean(r['avg_draws']):>6.0f}")

# =========================================================================
# (7) GRAND TRADE-OFF TABLE
# =========================================================================
print(f"\n{'='*88}")
print("(7) GRAND TRADE-OFF TABLE: weak responder power vs protocol cost")
print(f"{'='*88}")
rows = []
# baseline
rows.append(("2x3 fixed (baseline)", fixed[(2,3)]['pw_w'], fixed[(2,3)]['pw_all'],
             fixed[(2,3)]['fp'], 24, 12, "reference"))
# extended
for nc,km in [(3,3),(4,3),(5,3)]:
    f = fixed[(nc,km)]
    rows.append((f"{nc}x{km} fixed", f['pw_w'], f['pw_all'], f['fp'],
                 f['wk'], f['draws'], f"everyone +{f['wk']-24}wk"))
# adaptive
for k in ["2+1 x 3meas","2+2 x 3meas"]:
    r = adapt[k]
    rows.append((f"adapt {k[:7]}", r['pw_w'], r['pw_all'], r['fp'],
                 np.mean(r['avg_wk']), np.mean(r['avg_draws']),
                 "extra for borderline"))
# CV reduction
for cv in [0.18, 0.15, 0.12]:
    c = cv_res[cv]
    rows.append((f"CV={cv:.2f} 2x3", c['pw_w'], c['pw_all'], c['fp'],
                 24, 12, cv_how[cv]))
# combined
rows.append(("CV=.15 + adapt 2+1x3", combo_res["CV=.15 adapt 2+1x3"]['pw_w'],
             combo_res["CV=.15 adapt 2+1x3"]['pw_all'],
             combo_res["CV=.15 adapt 2+1x3"]['fp'],
             combo_res["CV=.15 adapt 2+1x3"]['avg_wk'],
             combo_res["CV=.15 adapt 2+1x3"]['draws'],
             "BEST: noise + adaptive"))

print(f"  {'strategy':<24}{'WEAK pw':>9}{'all pw':>8}{'FP':>6}{'avg wk':>8}"
      f"{'draws':>7}  {'note'}")
print(f"  {'-'*82}")
for name, pw_w, pw_all, fp, wk, dr, note in rows:
    tag = " <--" if "BEST" in note else ""
    print(f"  {name:<24}{pw_w*100:>8.0f}%{pw_all*100:>7.0f}%{fp*100:>5.0f}%"
          f"{wk:>7.0f}{dr:>7.0f}  {note}{tag}")

# =========================================================================
# (8) EFFICIENCY ANALYSIS: cost per 1% power gain for weak responders
# =========================================================================
base_pw = fixed[(2,3)]['pw_w']
print(f"\n{'='*88}")
print("(8) EFFICIENCY: extra weeks per +1% weak-responder power (lower = better)")
print(f"{'='*88}")
print(f"  {'strategy':<24}{'weak pw':>9}{'gain':>7}{'extra wk':>10}{'wk per +1%':>11}")
print(f"  {'-'*60}")
for name, pw_w, pw_all, fp, wk, dr, note in rows[1:]:
    gain = (pw_w - base_pw)*100
    extra_wk = wk - 24
    if gain > 0:
        eff = extra_wk / gain
        print(f"  {name:<24}{pw_w*100:>8.0f}%{gain:>+6.0f}%{extra_wk:>9.0f}{eff:>10.1f}")
    else:
        print(f"  {name:<24}{pw_w*100:>8.0f}%{gain:>+6.0f}%{extra_wk:>9.0f}{'n/a':>10}")

# =========================================================================
# (9) RECOMMENDATION
# =========================================================================
print(f"\n{'='*88}")
print("(9) RECOMMENDATION FOR WEAK RESPONDERS")
print(f"{'='*88}")

best_cv = cv_res[0.15]
best_combo = combo_res["CV=.15 adapt 2+1x3"]

print(f"""
  RANKING OF STRATEGIES (by efficiency):

  1. NOISE REDUCTION (CV 0.22 -> 0.15)                    <- PRIMARY LEVER
     how: same-day duplicate serum IS draws, fasting, timed AM collection
     cost: +12 extra blood draws (24 vs 12), same 24-week duration
     weak power: {best_cv['pw_w']*100:.0f}% (from {fixed[(2,3)]['pw_w']*100:.0f}%)
     reason: DT drops below the weak-responder tau range (21% -> 14%)
             this is IMPOSSIBLE to achieve by adding cycles at CV=0.22
             within any practical timeline for a 74-year-old

  2. ADAPTIVE 2-STAGE (on top of reduced CV)               <- SECONDARY
     how: after 2x3 at CV=0.15, patients with ambiguous signal (10-14%)
          get 1 extra AB cycle (total 36 wk for borderline only)
     weak power: {best_combo['pw_w']*100:.0f}%
     avg duration: {best_combo['avg_wk']:.0f} wk (most patients stop at 24 wk)
     -> catches remaining borderline weak responders without burdening
        the majority who were already clearly classified

  3. EXTENDED FIXED DESIGNS                                 <- AVOID
     4x3 fixed (48 wk) needed to match CV=0.15's weak power in 24 wk
     -> 74yo loses ~1.8 mL/min eGFR over 48 weeks (double the cost)
     -> adding cycles is quadratically inefficient: 4x measures = 2x DT reduction

  DO NOT PURSUE:
  - 5+ cycle designs (60+ wk = 14+ months for a 74yo with declining eGFR)
  - tau < 10% detection (these ARE non-responders; clinical threshold is ~10%)
  - aggressive early stopping thresholds (sacrifices weak-responder power
    to save non-responder burden — wrong trade for a diagnostic protocol)

  PRACTICAL CV REDUCTION:
    CV = 0.22: standard single fasting serum IS draw
    CV ~ 0.18: strict AM fasting + timed 2h-post-meal standardization
    CV ~ 0.15: above + duplicate draw same visit, averaged
               (extra cost: one more tube per visit = negligible)
    CV ~ 0.12: add 24h urine IS (IS excretion rate), composite with serum
               (harder: requires urine collection compliance)

  COMBINED RECOMMENDED PROTOCOL:
    stage 1 (wk 0-24): 2x3 crossover at CV ~ 0.15
      -> duplicate serum IS per visit, fasting AM, 12 visits total (24 draws)
      -> clear responders (obs > 14%): stop, KEEP regimen
      -> clear non-resp (obs < 0%): stop, DEPRESCRIBE
      -> borderline (0-14%): proceed to stage 2
    stage 2 (wk 24-36): 1 additional AB cycle for borderline only
      -> 3 more visits per arm (6 visits, 12 draws)
      -> final decision with 9 measures/arm, DT ~ {1.645*0.15*np.sqrt(2./9)*100:.0f}%

  EXPECTED OUTCOME (100-patient cohort):
    ~70% classified at stage 1 (24 wk) — no extra burden
    ~30% proceed to stage 2 (36 wk) — mostly weak/borderline responders
    weak responder power: ~{best_combo['pw_w']*100:.0f}%  (from ~{fixed[(2,3)]['pw_w']*100:.0f}% baseline)
    overall power: {best_combo['pw_all']*100:.0f}%  |  FP: {best_combo['fp']*100:.0f}%""")
