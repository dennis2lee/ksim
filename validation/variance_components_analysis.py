"""
Variance components and protocol attrition (review item R4).

The manuscript modelled measurement error as one lumped coefficient of
variation and then described a "package" of dietary standardization, fixed
sampling time, fasting and duplicate assays as lowering it from 0.25 to 0.15.
That conflates components that respond to different interventions. Repeating an
assay on the same tube averages down only the analytical component; the
day-to-day biological variation of the analyte and the pre-analytical handling
variation are properties of the specimen and are untouched by re-assaying it:

    cv_total^2 = cv_bio^2 + cv_pre^2 + cv_analytical^2 / n_replicates

The consequence is a hard floor. Pretorius et al. estimated the within-person
biological CV of total serum indoxyl sulfate at 35.9% in healthy volunteers. If
even a third of that variance is irreducible day-to-day biology, no assay
protocol reaches a total CV of 0.15. This script makes that floor explicit and
sweeps the assumptions that determine where it sits.

None of the component values below are measured in CKD. They are scenario
assumptions, and the manuscript now describes CV 0.15 as a hypothetical
scenario rather than an achievable target.

Also reports protocol attrition: missing measurements, dropout, partial
adherence and an intercurrent-event shock.

Outputs feed manuscript Table 4 and Tables S11, S12.
"""

import numpy as np
from nof1_core import (make_cohort, run_protocol, total_cv,
                       calibrate_boundaries)

N_PAT = 1000
COHORT_SEED = 777
THETA = 0.10
REP_SEEDS = range(700, 730)

cohort = make_cohort(N_PAT, COHORT_SEED, 0.18)
_cal = calibrate_boundaries(cohort, cv=0.15, null_margin=THETA,
                            alpha_total=0.05, seeds=range(900, 920),
                            order='randomized')
PRIMARY = dict(theta=THETA, null_margin=THETA,
               alpha1=_cal['alpha1'], alpha2=_cal['alpha2'])


def banner(t):
    print('\n' + '=' * 84)
    print(t)
    print('=' * 84)


# ---------------------------------------------------------------------------
# Component decomposition
# ---------------------------------------------------------------------------
banner('R4. WHAT TOTAL CV IS REACHABLE, BY VARIANCE COMPONENT')
print('  cv_bio is split into a part that sampling standardization can remove')
print('  (diet, time of day, fasting state) and an irreducible day-to-day part.')
print('  Duplicate assay divides only the analytical variance.\n')

# Scenarios span how much of the 35.9% within-person biological CV reported by
# Pretorius et al. is attributable to standardizable sources.
COMPONENTS = [
    # label, cv_bio_irreducible, cv_bio_standardizable, cv_pre, cv_analytical
    ('optimistic: 2/3 of biological CV is standardizable',
     0.207, 0.293, 0.08, 0.06),
    ('intermediate: 1/2 standardizable',
     0.254, 0.254, 0.08, 0.06),
    ('pessimistic: 1/3 standardizable',
     0.293, 0.207, 0.08, 0.06),
]

print(f'  {"scenario":<52}{"none":>8}{"std.":>8}{"std.+dup":>10}')
print('  ' + '-' * 78)
reachable = {}
for label, bio_irr, bio_std, pre, ana in COMPONENTS:
    no_std = total_cv(cv_bio=np.hypot(bio_irr, bio_std), cv_pre=pre,
                      cv_analytical=ana, n_replicates=1)
    std = total_cv(cv_bio=bio_irr, cv_pre=pre * 0.5,
                   cv_analytical=ana, n_replicates=1)
    std_dup = total_cv(cv_bio=bio_irr, cv_pre=pre * 0.5,
                       cv_analytical=ana, n_replicates=2)
    reachable[label] = (no_std, std, std_dup)
    print(f'  {label:<52}{no_std:>8.3f}{std:>8.3f}{std_dup:>10.3f}')

print('\n  The duplicate assay column differs from the column to its left by at')
print('  most {:.4f}. Re-assaying the same specimen cannot deliver the'
      .format(max(abs(v[1] - v[2]) for v in reachable.values())))
print('  0.25 -> 0.15 reduction the manuscript previously attributed to a')
print('  package that included it. Reaching 0.15 requires the irreducible')
print('  biological component alone to be below 0.15.')

banner('R4. WHAT IRREDUCIBLE BIOLOGICAL CV WOULD 0.15 TOTAL REQUIRE?')
for pre, ana, k in ((0.04, 0.06, 2), (0.08, 0.06, 2), (0.08, 0.06, 1)):
    resid = 0.15 ** 2 - pre ** 2 - ana ** 2 / k
    val = np.sqrt(resid) if resid > 0 else float('nan')
    print(f'  cv_pre={pre:.2f}, cv_analytical={ana:.2f}, replicates={k}: '
          f'cv_bio must be <= {val:.3f}')
print('  Against a reported within-person biological CV of 0.359, that is a')
print('  reduction of more than half, for which there is no CKD evidence.')

# ---------------------------------------------------------------------------
# Operating characteristics across the plausible CV range
# ---------------------------------------------------------------------------
banner('R4. OPERATING CHARACTERISTICS ACROSS THE PLAUSIBLE CV RANGE')
print('  Primary rule (matched null, alpha-spending). CV 0.15 is shown as one')
print('  hypothetical scenario among several, not as a representative value.\n')
print(f'  {"CV":>6}{"basis":<40}{"sens":>8}{"spec":>8}{"weak":>8}{"S2":>6}')
print('  ' + '-' * 76)
CV_BASIS = [
    (0.15, 'hypothetical best case (see above)'),
    (0.20, 'strong standardization'),
    (0.25, 'reference simulation assumption'),
    (0.30, 'modest standardization'),
    (0.36, 'Pretorius within-person biological CV'),
]
for cv, basis in CV_BASIS:
    runs = [run_protocol(cohort, cv=cv, measurement_seed=s, order='randomized',
                         **PRIMARY) for s in REP_SEEDS]
    print(f'  {cv:>5.2f}  {basis:<38}'
          f'{np.mean([r["sens"] for r in runs])*100:>7.1f}%'
          f'{np.mean([r["spec"] for r in runs])*100:>7.1f}%'
          f'{np.mean([r["detect"]["weak"] for r in runs])*100:>7.0f}%'
          f'{np.mean([r["stage2_rate"] for r in runs])*100:>5.0f}%')

# ---------------------------------------------------------------------------
# Attrition
# ---------------------------------------------------------------------------
banner('R4. PROTOCOL ATTRITION AND INTERCURRENT EVENTS')
print('  Primary rule, randomized order, CV = 0.25 (reference assumption).\n')
print(f'  {"scenario":<44}{"analysed":>10}{"sens":>8}{"spec":>8}')
print('  ' + '-' * 70)

ATTRITION = [
    ('complete data', dict()),
    ('10% of measurements missing', dict(p_missing=0.10)),
    ('20% of measurements missing', dict(p_missing=0.20)),
    ('10% dropout', dict(p_dropout=0.10)),
    ('20% dropout', dict(p_dropout=0.20)),
    ('20% missing + 20% dropout', dict(p_missing=0.20, p_dropout=0.20)),
]
for label, kw in ATTRITION:
    runs = [run_protocol(cohort, cv=0.25, measurement_seed=s,
                         order='randomized', **PRIMARY, **kw)
            for s in REP_SEEDS]
    print(f'  {label:<44}{np.mean([r["n_analysed"] for r in runs]):>10.0f}'
          f'{np.mean([r["sens"] for r in runs])*100:>7.1f}%'
          f'{np.mean([r["spec"] for r in runs])*100:>7.1f}%')

print('\n  Partial adherence is modelled as an attenuated treatment effect and')
print('  an intercurrent event as an inflated CV in a subset of patients.')
print(f'\n  {"scenario":<44}{"sens":>8}{"spec":>8}')
print('  ' + '-' * 60)

# partial adherence: a fraction of patients receive an attenuated effect
for frac, keep in ((0.20, 0.5), (0.30, 0.5), (0.20, 0.0)):
    rng = np.random.default_rng(11)
    c2 = make_cohort(N_PAT, COHORT_SEED, 0.18)
    idx = rng.random(N_PAT) < frac
    c2.tau = c2.tau.copy()
    c2.tau[idx] = c2.tau[idx] * keep
    runs = [run_protocol(c2, cv=0.25, measurement_seed=s, order='randomized',
                         **PRIMARY) for s in REP_SEEDS]
    lbl = (f'{frac:.0%} of patients at {keep:.0%} adherence')
    print(f'  {lbl:<44}'
          f'{np.mean([r["sens"] for r in runs])*100:>7.1f}%'
          f'{np.mean([r["spec"] for r in runs])*100:>7.1f}%')
    # note: true responder status is judged on the delivered effect, so this
    # measures the protocol's behaviour when the delivered effect differs from
    # the one the patient would have had under full adherence

# intercurrent event: CV inflated in a subset
for frac, mult in ((0.10, 2.0), (0.20, 2.0)):
    rng = np.random.default_rng(12)
    cvp = np.full(N_PAT, 0.25)
    hit = rng.random(N_PAT) < frac
    cvp[hit] *= mult
    runs = [run_protocol(cohort, cv=0.25, cv_per_patient=cvp,
                         measurement_seed=s, order='randomized', **PRIMARY)
            for s in REP_SEEDS]
    lbl = f'{frac:.0%} with intercurrent event (CV x{mult:.0f})'
    print(f'  {lbl:<44}'
          f'{np.mean([r["sens"] for r in runs])*100:>7.1f}%'
          f'{np.mean([r["spec"] for r in runs])*100:>7.1f}%')

print('\n  Reproduce with: python3 validation/variance_components_analysis.py')
