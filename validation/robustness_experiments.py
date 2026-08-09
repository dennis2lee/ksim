"""
Robustness of the protocol to violations of the reference model's assumptions.

Every scenario below is generated and run through `nof1_core`, so the reference
cohort here is byte-identical to the one used by every other table. Previously
this script carried its own copy of the cohort generator and each scenario
re-derived the cohort with hand-aligned random draws; that is what let the
scripts drift apart.

Scenarios:
  0   reference (Gaussian, independent, homoscedastic, no carryover)
  A   log-normal baseline IS (right-skewed)
  B   bimodal effect distribution (sharp responder / non-responder split)
  C   heteroscedastic CV across patients
  D   serial correlation within a measurement period, AR(1) rho = 0.3
  D'  serial correlation, rho = 0.5
  E   partial carryover, 20% of the effect persists into the control period
  F   combined worst case (A + B + D + E)
  G   eGFR-decline drift applied over the measurement schedule

Results are reported under both decision rules so the reader can see which
conclusions depend on which null (review item R1):

  zero-null      H0: tau <= 0,    alpha 0.05 per look (the original rule)
  primary        H0: tau <= 0.10, alpha-spending calibrated to 5% overall,
                 randomized cycle order (review items R1, R2, R3)
"""

import numpy as np
from nof1_core import make_cohort, run_protocol, calibrate_boundaries

N_PAT = 1000
COHORT_SEED = 777
MEAS_SEED = 777
CV = 0.15
THETA = 0.10

reference = make_cohort(N_PAT, COHORT_SEED, 0.18)
_cal = calibrate_boundaries(reference, cv=CV, null_margin=THETA,
                            alpha_total=0.05, seeds=range(900, 920),
                            order='randomized')

ZERO_NULL = dict(theta=THETA, null_margin=0.0, alpha1=0.05, alpha2=0.05,
                 order='fixed')
PRIMARY = dict(theta=THETA, null_margin=THETA, alpha1=_cal['alpha1'],
               alpha2=_cal['alpha2'], order='randomized')

# Heteroscedastic CV: patient-level CV drawn from a heavy-tailed distribution
_rng = np.random.default_rng(781)
CV_HETERO = np.clip(CV * _rng.lognormal(0, 0.45, N_PAT), 0.05, 0.60)

lognormal_cohort = make_cohort(N_PAT, COHORT_SEED, 0.18,
                               baseline_shape='lognormal')
bimodal_cohort = make_cohort(N_PAT, COHORT_SEED, 0.18, tau_shape='bimodal')
worst_cohort = make_cohort(N_PAT, COHORT_SEED, 0.18,
                           baseline_shape='lognormal', tau_shape='bimodal')

SCENARIOS = [
    ('0', 'reference', reference, dict()),
    ('A', 'log-normal baseline IS', lognormal_cohort, dict()),
    ('B', 'bimodal effect distribution', bimodal_cohort, dict()),
    ('C', 'heteroscedastic CV', reference, dict(cv_per_patient=CV_HETERO)),
    ('D', 'serial correlation, rho = 0.3', reference, dict(rho=0.3)),
    ("D'", 'serial correlation, rho = 0.5', reference, dict(rho=0.5)),
    ('E', 'carryover 20%', reference, dict(carryover=0.20)),
    ('F', 'combined worst case', worst_cohort, dict(rho=0.3, carryover=0.20)),
    ('G', 'eGFR-decline drift', reference, dict(drift=True)),
]


def banner(t):
    print('\n' + '=' * 86)
    print(t)
    print('=' * 86)


banner('TABLE 3. ROBUSTNESS TO DISTRIBUTIONAL AND DESIGN VIOLATIONS')
print(f'  Reference cohort seed {COHORT_SEED}, measurement seed {MEAS_SEED}, '
      f'CV = {CV}.')
print(f'  Primary rule: H0 tau <= {THETA:.2f}, alpha1 = {_cal["alpha1"]:.4f}, '
      f'alpha2 = {_cal["alpha2"]:.4f}, randomized order.\n')
print(f'  {"":<4}{"scenario":<32}| {"zero-null (original)":^22}'
      f'| {"primary rule":^22}')
print(f'  {"":<4}{"":<32}| {"sens":>7}{"spec":>7}{"S2":>7}'
      f'| {"sens":>7}{"spec":>7}{"S2":>7}')
print('  ' + '-' * 84)
for tag, label, coh, kw in SCENARIOS:
    z = run_protocol(coh, cv=CV, measurement_seed=MEAS_SEED, **ZERO_NULL, **kw)
    p = run_protocol(coh, cv=CV, measurement_seed=MEAS_SEED, **PRIMARY, **kw)
    print(f'  {tag:<4}{label:<32}| {z["sens"]*100:>6.1f}%{z["spec"]*100:>6.1f}%'
          f'{z["stage2_rate"]*100:>6.0f}%| {p["sens"]*100:>6.1f}%'
          f'{p["spec"]*100:>6.1f}%{p["stage2_rate"]*100:>6.0f}%')

banner('WEAK-RESPONDER DETECTION BY SCENARIO')
print(f'  {"":<4}{"scenario":<32}{"zero-null":>12}{"primary":>12}')
print('  ' + '-' * 62)
for tag, label, coh, kw in SCENARIOS:
    z = run_protocol(coh, cv=CV, measurement_seed=MEAS_SEED, **ZERO_NULL, **kw)
    p = run_protocol(coh, cv=CV, measurement_seed=MEAS_SEED, **PRIMARY, **kw)
    print(f'  {tag:<4}{label:<32}{z["detect"]["weak"]*100:>11.0f}%'
          f'{p["detect"]["weak"]*100:>11.0f}%')

banner('SUMMARY')
zs = [run_protocol(c, cv=CV, measurement_seed=MEAS_SEED, **ZERO_NULL, **k)
      for _, _, c, k in SCENARIOS]
ps = [run_protocol(c, cv=CV, measurement_seed=MEAS_SEED, **PRIMARY, **k)
      for _, _, c, k in SCENARIOS]
print(f'  zero-null: sensitivity {min(r["sens"] for r in zs)*100:.1f}-'
      f'{max(r["sens"] for r in zs)*100:.1f}%, specificity '
      f'{min(r["spec"] for r in zs)*100:.1f}-{max(r["spec"] for r in zs)*100:.1f}%')
print(f'  primary:   sensitivity {min(r["sens"] for r in ps)*100:.1f}-'
      f'{max(r["sens"] for r in ps)*100:.1f}%, specificity '
      f'{min(r["spec"] for r in ps)*100:.1f}-{max(r["spec"] for r in ps)*100:.1f}%')
print('\n  Reference plus eight stress conditions (nine rows in total).')
print('  Reproduce with: python3 validation/robustness_experiments.py')
