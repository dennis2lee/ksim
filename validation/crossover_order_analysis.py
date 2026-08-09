"""
Crossover order and period effects (review item R3).

Every cycle in the original protocol ran control (A) before treatment (B), so a
treatment effect and a time trend enter the arm contrast the same way. Showing
that one monotonic eGFR-decline scenario happened to bias the result
conservatively does not establish that the fixed order is safe: seasonality,
adherence learning, regression to the mean and visit fatigue are not monotonic,
and some of them bias the other way.

This script compares four cycle orders

  fixed           A B A B   (original)
  reverse         B A B A
  counterbalanced A B B A   (both arms centred on week 12)
  randomized      AB or BA drawn independently for each cycle

under five time-trend scenarios and two estimators (the arm-mean ratio, and a
per-patient model carrying an explicit linear week term). The quantity that
matters most is the bias of the effect estimate, which does not depend on the
decision rule; operating characteristics are reported under the primary rule
for context.

Outputs feed manuscript Table 5 and Table S10.
"""

import numpy as np
from nof1_core import (make_cohort, run_protocol, calibrate_boundaries,
                       linear_period_effect, seasonal_period_effect)

N_PAT = 1000
COHORT_SEED = 777
THETA = 0.10
REP_SEEDS = range(700, 720)          # 20 replications; bias converges quickly

cohort = make_cohort(N_PAT, COHORT_SEED, 0.18)

# Primary rule: matched null, alpha-spending calibrated to 5% overall.
_cal = calibrate_boundaries(cohort, cv=0.15, null_margin=THETA,
                            alpha_total=0.05, seeds=range(900, 920),
                            order='randomized')
PRIMARY = dict(theta=THETA, null_margin=THETA,
               alpha1=_cal['alpha1'], alpha2=_cal['alpha2'])

ORDERS = ['fixed', 'reverse', 'counterbalanced', 'randomized']

SCENARIOS = [
    ('no time trend', dict()),
    ('eGFR decline (model drift)', dict(drift=True)),
    ('linear +0.5%/week', dict(period_effect=linear_period_effect(0.005))),
    ('linear -0.5%/week', dict(period_effect=linear_period_effect(-0.005))),
    ('seasonal, 10% amplitude', dict(period_effect=seasonal_period_effect(0.10))),
    ('carryover 20%', dict(carryover=0.20)),
]


def evaluate(order, estimator, scen_kw):
    """Mean bias in the effect estimate plus operating characteristics."""
    bias, sens, spec = [], [], []
    for s in REP_SEEDS:
        r = run_protocol(cohort, cv=0.15, measurement_seed=s, order=order,
                         estimator=estimator, **PRIMARY, **scen_kw)
        ok = ~np.isnan(r['obs_final'])
        bias.append(float(np.mean(r['obs_final'][ok] - cohort.tau[ok])))
        sens.append(r['sens'])
        spec.append(r['spec'])
    return (float(np.mean(bias)), float(np.mean(sens)), float(np.mean(spec)))


def banner(t):
    print('\n' + '=' * 88)
    print(t)
    print('=' * 88)


banner('R3. EFFECT-ESTIMATE BIAS BY CYCLE ORDER AND TIME TREND '
       '(arm-mean ratio estimator)')
print('  Bias = mean(estimated reduction - true tau), percentage points.')
print(f'  {"scenario":<30}' + ''.join(f'{o:>15}' for o in ORDERS))
print('  ' + '-' * 86)
ratio_bias = {}
for label, kw in SCENARIOS:
    cells = []
    for o in ORDERS:
        b, se, sp = evaluate(o, 'ratio', kw)
        ratio_bias[(label, o)] = (b, se, sp)
        cells.append(f'{b*100:>+14.2f}')
    print(f'  {label:<30}' + ''.join(cells))

banner('R3. THE SAME, WITH AN EXPLICIT PER-PATIENT LINEAR WEEK TERM')
print('  Estimator regresses log(IS) on treatment and week within each patient.')
print(f'  {"scenario":<30}' + ''.join(f'{o:>15}' for o in ORDERS))
print('  ' + '-' * 86)
ols_bias = {}
for label, kw in SCENARIOS:
    cells = []
    for o in ORDERS:
        b, se, sp = evaluate(o, 'ols_period', kw)
        ols_bias[(label, o)] = (b, se, sp)
        cells.append(f'{b*100:>+14.2f}')
    print(f'  {label:<30}' + ''.join(cells))

banner('R3. OPERATING CHARACTERISTICS UNDER THE PRIMARY RULE '
       '(matched null, alpha-spending)')
print(f'  {"scenario":<30}{"order":<17}{"sens":>8}{"spec":>8}'
      f'{"bias(pp)":>10}')
print('  ' + '-' * 73)
for label, kw in SCENARIOS:
    for o in ORDERS:
        b, se, sp = ratio_bias[(label, o)]
        print(f'  {label:<30}{o:<17}{se*100:>7.1f}%{sp*100:>7.1f}%'
              f'{b*100:>+10.2f}')

banner('SUMMARY')
worst = {}
for o in ORDERS:
    mags = [abs(ratio_bias[(l, o)][0]) for l, _ in SCENARIOS]
    worst[o] = max(mags) * 100
print('  Worst-case absolute bias across the six scenarios, ratio estimator:')
for o in ORDERS:
    print(f'    {o:<17}{worst[o]:>6.2f} pp')
worst_ols = {o: max(abs(ols_bias[(l, o)][0]) for l, _ in SCENARIOS) * 100
             for o in ORDERS}
print('  Worst-case absolute bias, with an explicit week term:')
for o in ORDERS:
    print(f'    {o:<17}{worst_ols[o]:>6.2f} pp')
print('\n  Reproduce with: python3 validation/crossover_order_analysis.py')
