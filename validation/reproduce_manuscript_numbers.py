"""
Reproduce every number quoted in the manuscript body.

This script used to carry its own copy of the cohort generator plus a block of
dummy draws whose only purpose was to push the random stream into the same
state as `robustness_experiments.py`. One of those dummy lines reused the
variable `_` for both the eGFR slope and the eGFR itself and then raised a
negative number to a fractional power, producing a RuntimeWarning on every run.
All of that is gone: the cohort and the protocol now come from `nof1_core`, so
there is nothing to keep in sync.

Both decision rules are printed. The zero-null column is what earlier versions
of the manuscript reported; the primary column is the rule whose null matches
the 10% responder definition, with alpha-spending and randomized cycle order.
"""

import numpy as np
from nof1_core import (make_cohort, run_protocol, critical_value,
                       calibrate_boundaries, baseline_moments)

N_PAT = 1000
COHORT_SEED = 777
MEAS_SEED = 777
CV = 0.15
THETA = 0.10

cohort = make_cohort(N_PAT, COHORT_SEED, 0.18)
cal = calibrate_boundaries(cohort, cv=CV, null_margin=THETA,
                           alpha_total=0.05, seeds=range(900, 920))

ZERO_NULL = dict(theta=THETA, null_margin=0.0, alpha1=0.05, alpha2=0.05,
                 order='fixed')
PRIMARY = dict(theta=THETA, null_margin=THETA, alpha1=cal['alpha1'],
               alpha2=cal['alpha2'], order='randomized')


def banner(t):
    print('\n' + '=' * 74)
    print(t)
    print('=' * 74)


banner('COHORT CHARACTERISTICS')
m, s = baseline_moments(cohort)
print(f'  baseline IS: realized mean {m:.3f} ug/mL, SD {s:.3f} ug/mL')
print('  (the standardization imposes 5.4 / 3.6 before the 0.5 floor; the')
print('   floor moves 32 of 1000 values and shifts the moments slightly)')
st = cohort.strata(THETA)
for k in ('non_responder', 'weak', 'moderate', 'strong'):
    print(f'  {k:<16}{int(st[k].sum()):>6}')
print(f'  {"total":<16}{N_PAT:>6}')
print(f'  mean tau among true responders: '
      f'{cohort.tau[cohort.true_responder(THETA)].mean()*100:.0f}%')

banner('TABLE 1. CRITICAL VALUE BY DESIGN, CV AND NULL')
print(f'  {"design":<8}{"CV":>6}{"n/arm":>7}{"weeks":>7}'
      f'{"H0 tau<=0":>12}{"H0 tau<=0.10":>15}')
print('  ' + '-' * 60)
for design, cv, n, wk in (('2x3', 0.25, 6, 24), ('2x3', 0.15, 6, 24),
                          ('3x3', 0.25, 9, 36), ('4x3', 0.25, 12, 48),
                          ('3x3', 0.15, 9, 36)):
    print(f'  {design:<8}{cv:>6.2f}{n:>7}{wk:>7}'
          f'{critical_value(cv, n, 0.05, 0.0)*100:>11.1f}%'
          f'{critical_value(cv, n, 0.05, THETA)*100:>14.1f}%')

banner('TABLE 2. SINGLE-RUN CLASSIFICATION')
for label, kw in (('zero-null (previous manuscript rule)', ZERO_NULL),
                  ('primary rule', PRIMARY)):
    r = run_protocol(cohort, cv=CV, measurement_seed=MEAS_SEED, **kw)
    print(f'\n  {label}')
    print(f'    true responders {r["n_responders"]}, '
          f'true non-responders {r["n_nonresponders"]}')
    print(f'    TP={r["tp"]}  FN={r["fn"]}  FP={r["fp"]}  TN={r["tn"]}')
    print(f'    sensitivity {r["sens"]*100:.1f}%, '
          f'specificity {r["spec"]*100:.1f}%')
    print(f'    Stage 2 rate {r["stage2_rate"]*100:.0f}%, '
          f'mean duration {r["mean_weeks"]:.1f} weeks')
    print(f'    subgroup detection: '
          + ', '.join(f'{k} {v*100:.0f}%'
                      for k, v in r['detect'].items() if k != 'non_responder'))

banner('OVERALL TYPE I ERROR AT THE NULL BOUNDARY')
print('  Measured on a cohort with every tau set to the null margin, so this')
print('  is the actual type I error, not the false-positive rate among true')
print('  non-responders (which mixes in the definition mismatch).\n')
for margin, a1, a2, label in ((0.0, 0.05, 0.05, 'zero-null, uncorrected'),
                              (THETA, 0.05, 0.05, 'matched-null, uncorrected'),
                              (THETA, cal['alpha1'], cal['alpha2'],
                               'matched-null, alpha-spending')):
    nc = make_cohort(N_PAT, COHORT_SEED, 0.18)
    nc.tau = np.full(N_PAT, margin)
    rates = [run_protocol(nc, cv=CV, measurement_seed=sd,
                          theta=margin + 1e-9, null_margin=margin,
                          alpha1=a1, alpha2=a2)['cls']
             for sd in range(900, 950)]
    print(f'  {label:<34}{np.mean([(c == "R").mean() for c in rates])*100:>6.1f}%')

banner('BASELINE DISTRIBUTION IS IRRELEVANT BY CONSTRUCTION')
ln = make_cohort(N_PAT, COHORT_SEED, 0.18, baseline_shape='lognormal')
a = run_protocol(cohort, cv=CV, measurement_seed=MEAS_SEED, **ZERO_NULL)
b = run_protocol(ln, cv=CV, measurement_seed=MEAS_SEED, **ZERO_NULL)
print(f'  Gaussian baseline mean {cohort.base_is.mean():.2f}, '
      f'log-normal baseline mean {ln.base_is.mean():.2f}')
print(f'  classifications identical: {np.array_equal(a["cls"], b["cls"])}')
print(f'  largest difference in the estimate: '
      f'{np.nanmax(np.abs(a["obs_final"] - b["obs_final"])):.1e}')
print('  The estimator is (mean_A - mean_B)/mean_A, in which the patient\'s')
print('  own baseline cancels exactly. No baseline distribution can change a')
print('  classification, so this is an analytic invariance rather than an')
print('  empirical robustness result.')

banner('INLINE NUMBERS')
r = run_protocol(cohort, cv=CV, measurement_seed=MEAS_SEED, **PRIMARY)
print(f'  mean duration under the primary rule: {r["mean_weeks"]:.1f} weeks')
print(f'  critical value, H0 tau<=0.10, CV 0.15, n=6: '
      f'{critical_value(0.15, 6, 0.05, THETA)*100:.1f}%')
print(f'  critical value, H0 tau<=0.10, CV 0.25, n=6: '
      f'{critical_value(0.25, 6, 0.05, THETA)*100:.1f}%')
print(f'  calibrated alphas: {cal["alpha1"]:.4f} / {cal["alpha2"]:.4f}')
print('\n  Reproduce with: python3 validation/reproduce_manuscript_numbers.py')
