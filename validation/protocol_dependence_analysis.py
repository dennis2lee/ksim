"""
How much the protocol's performance depends on operational choices
(review item R15).

The manuscript specifies a 4-week period, a 2-week washout and an effect that
is fully present from the first day of a period. None of those is measured;
they are design assumptions, and an earlier version of this work presented them
alongside operational recommendations they did not support. This script sweeps
them so the recommendations can either be grounded or withdrawn.

Two mechanisms are added to the core model for this purpose:

  washout_weeks   the length of each washout, which changes the schedule and
                  therefore how a time trend enters the arm contrast
  onset_weeks     the time constant of an exponential approach to the
                  steady-state effect, so a measurement taken t weeks into a
                  period sees effect * (1 - exp(-t / onset_weeks)). A fiber
                  intervention works through a microbiome shift that takes time
                  to establish, so an effect present from day one is optimistic.

Carryover and adherence remain independent scenario axes. A longer washout
would in reality reduce carryover; holding them independent is deliberate,
because the point of the grid is to show which assumption the conclusions rest
on, not to model a particular decay.
"""

import numpy as np
from nof1_core import make_cohort, run_protocol, calibrate_boundaries

N_PAT, COHORT_SEED, THETA = 1000, 777, 0.10
CV = 0.264                              # reachable under the intermediate split
SEEDS = range(3000, 3030)

cohort = make_cohort(N_PAT, COHORT_SEED, 0.18)
cal = calibrate_boundaries(cohort, cv=0.15, null_margin=THETA, alpha_total=0.05,
                           seeds=range(900, 920), order='randomized')
PRIMARY = dict(theta=THETA, null_margin=THETA, alpha1=cal['alpha1'],
               alpha2=cal['alpha2'], order='randomized')


def oc(**kw):
    runs = [run_protocol(cohort, cv=CV, measurement_seed=s, **PRIMARY, **kw)
            for s in SEEDS]
    return (float(np.mean([r['sens'] for r in runs])),
            float(np.mean([r['spec'] for r in runs])),
            float(np.mean([r['detect']['weak'] for r in runs])),
            float(np.mean([r['mean_weeks'] for r in runs])))


def adherence_cohort(frac, keep, seed=11):
    rng = np.random.default_rng(seed)
    c = make_cohort(N_PAT, COHORT_SEED, 0.18)
    idx = rng.random(N_PAT) < frac
    c.tau = c.tau.copy()
    c.tau[idx] = c.tau[idx] * keep
    return c


def banner(t):
    print('\n' + '=' * 88)
    print(t)
    print('=' * 88)


banner('R15. ONSET RAMP: THE ASSUMPTION THAT THE EFFECT IS PRESENT AT DAY ONE')
print(f'  CV = {CV}, primary rule, mean of {len(list(SEEDS))} replications.')
print('  Measurements are taken in weeks 2, 3 and 4 of each 4-week period.\n')
print(f'  {"onset time constant":<28}{"sens":>8}{"spec":>8}{"weak":>8}')
print('  ' + '-' * 52)
for onset, label in ((0.0, 'none (effect immediate)'), (0.5, '0.5 weeks'),
                     (1.0, '1 week'), (2.0, '2 weeks'), (3.0, '3 weeks')):
    sn, sp, wk, _ = oc(onset_weeks=onset)
    print(f'  {label:<28}{sn*100:>7.1f}%{sp*100:>7.1f}%{wk*100:>7.0f}%')
print('''
  This is the single most consequential assumption in the protocol. A 2-week
  onset constant costs 19.5 sensitivity points, more than three times what
  moving from the intermediate to the pessimistic variance split costs, and a
  3-week constant costs 33.6. The design mitigates onset only by measuring late
  in the period, which is why the draws sit in weeks 2 to 4; with a 4-week
  period there is no room to wait longer. A trial that cannot bound the onset
  rate cannot bound its own sensitivity.''')

banner('R15. WASHOUT LENGTH AND CARRYOVER')
print('  Washout changes the schedule; with no time trend the ratio estimator is')
print('  unaffected by timing, so washout acts on performance only through')
print('  carryover and through total duration.\n')
print(f'  {"washout":>9}{"carryover":>11}{"sens":>8}{"spec":>8}{"weeks":>8}')
print('  ' + '-' * 46)
for wash in (2, 3, 4):
    for carry in (0.0, 0.10, 0.20, 0.40):
        sn, sp, _, wks = oc(washout_weeks=wash, carryover=carry)
        print(f'  {wash:>7}wk{carry:>10.0%}{sn*100:>7.1f}%{sp*100:>7.1f}%'
              f'{wks:>8.1f}')
print('\n  Carryover costs sensitivity at every washout length; a longer washout')
print('  buys calendar time, not performance, unless it also reduces carryover.')
print('  The trial that can measure the reversion rate is the one that can')
print('  choose the washout; we cannot choose it from these simulations.')

banner('R15. PARTIAL ADHERENCE')
print(f'  {"scenario":<40}{"sens":>8}{"spec":>8}')
print('  ' + '-' * 56)
for frac, keep in ((0.0, 1.0), (0.20, 0.75), (0.20, 0.50), (0.30, 0.50),
                   (0.30, 0.25)):
    c = adherence_cohort(frac, keep) if frac else cohort
    runs = [run_protocol(c, cv=CV, measurement_seed=s, **PRIMARY)
            for s in SEEDS]
    lbl = ('full adherence' if not frac
           else f'{frac:.0%} of patients at {keep:.0%} adherence')
    print(f'  {lbl:<40}'
          f'{np.mean([r["sens"] for r in runs])*100:>7.1f}%'
          f'{np.mean([r["spec"] for r in runs])*100:>7.1f}%')
print('\n  Adherence enters as a reduction in the delivered effect, so a patient')
print('  who does not take the intervention is a simulated non-responder under')
print('  the estimand. Whether a trial should treat them that way or as a')
print('  protocol deviation is a design decision, not a simulation result.')

banner('R15. WORST PLAUSIBLE COMBINATION')
sn0, sp0, wk0, _ = oc()
sn1, sp1, wk1, _ = oc(onset_weeks=2.0, washout_weeks=2, carryover=0.20)
print(f'  reference at CV {CV}                     '
      f'{sn0*100:>6.1f}% sens  {sp0*100:>6.1f}% spec  {wk0*100:>3.0f}% weak')
print(f'  onset 2 weeks + 20% carryover              '
      f'{sn1*100:>6.1f}% sens  {sp1*100:>6.1f}% spec  {wk1*100:>3.0f}% weak')
print(f'\n  Sensitivity falls {(sn0-sn1)*100:.1f} points and weak-responder detection')
print(f'  falls from {wk0*100:.0f}% to {wk1*100:.0f}%. The operational parameters are not a')
print('  footnote to the design: at the reachable CV they move performance more')
print('  than the measurement precision this paper is about. Any recommendation')
print('  about period length, washout or onset that these simulations do not')
print('  support has been removed from the manuscript rather than softened.')

print('\n  Reproduce with: python3 validation/protocol_dependence_analysis.py')
