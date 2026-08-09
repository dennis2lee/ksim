"""
Cohorts parameterized to published trials (review item R6).

The previous version of this analysis contained three source errors that
changed how the results should be read:

  1. A row labelled "EPPIC AST-120 arm (Schulman 2015)" used N = 460. EPPIC-1
     and EPPIC-2 randomized 2,035 patients (PMID 25349205). 460 is the sample
     size of CAP-KD, a different trial (Akizawa 2009, PMID 19615804).
  2. EPPIC does not report serum indoxyl sulfate at all; its endpoint was a
     composite of dialysis initiation, transplantation and creatinine
     doubling. The dose-dependent IS lowering attributed to it belongs to the
     Schulman 2006 dose-ranging study (PMID 16564934).
  3. A row labelled "2025 meta-analysis pooled (11 RCTs)" used N = 200. The
     indoxyl sulfate analysis of that meta-analysis pooled 11 trials and 398
     patients; 200 is the size of the intervention arm alone
     (PMC11860371, Table 2).

Every input is now listed explicitly with its source and its status. Only the
baseline IS mean and SD and the mean effect are taken from publications. The
SD of the individual effect distribution and the non-responder fraction cannot
be recovered from group-level reports and are author assumptions; they are
labelled as such in the output and in Table S3.

Where a Monte Carlo cohort size is chosen for simulation convenience rather
than taken from a trial, the `n_source` column says so. It is not presented as
a published cohort parameter.
"""

import numpy as np
from nof1_core import make_cohort, run_protocol, calibrate_boundaries

THETA = 0.10
PUB_SEEDS = range(700, 750)

_ref = make_cohort(1000, 777, 0.18)
_cal = calibrate_boundaries(_ref, cv=0.15, null_margin=THETA,
                            alpha_total=0.05, seeds=range(900, 920),
                            order='randomized')
PRIMARY = dict(theta=THETA, null_margin=THETA,
               alpha1=_cal['alpha1'], alpha2=_cal['alpha2'], order='randomized')

# ---------------------------------------------------------------------------
# Inputs. Fields:
#   label, n, n_source, tau_mean, tau_source, tau_sd, base_mean, base_sd,
#   base_source, nr_frac
# tau_sd and nr_frac are ASSUMPTIONS in every row; no published report gives
# individual-level effect SDs for serum IS.
# ---------------------------------------------------------------------------
COHORTS = [
    dict(label='Schulman 2006 dose-ranging (AST-120)',
         n=200, n_source='Monte Carlo size (trial N not extracted)',
         tau_mean=0.35,
         tau_source='dose-dependent IS lowering, AJKD 2006;47:565-577',
         tau_sd=0.15, base_mean=5.4, base_sd=3.6,
         base_source='Lin 2011 CKD stage 4 (assumed, not trial-reported)',
         nr_frac=0.15),
    dict(label='CAP-KD (Akizawa 2009)',
         n=460, n_source='trial N, AJKD 2009;54:459-467',
         tau_mean=0.35,
         tau_source='ASSUMED equal to Schulman 2006; CAP-KD reports no IS',
         tau_sd=0.15, base_mean=5.4, base_sd=3.6,
         base_source='Lin 2011 CKD stage 4 (assumed, not trial-reported)',
         nr_frac=0.15),
    dict(label='EPPIC-1 and EPPIC-2 (Schulman 2015)',
         n=2035, n_source='trial N, JASN 2015;26:1732-1746',
         tau_mean=0.35,
         tau_source='ASSUMED; EPPIC reports no IS outcome at all',
         tau_sd=0.15, base_mean=5.4, base_sd=3.6,
         base_source='Lin 2011 CKD stage 4 (assumed, not trial-reported)',
         nr_frac=0.15),
    dict(label='Rossi SYNERGY synbiotic (2016)',
         n=37, n_source='trial N, CJASN 2016;11:223-231',
         tau_mean=0.14, tau_source='reported IS reduction ~12-14%',
         tau_sd=0.10, base_mean=4.0, base_sd=2.5,
         base_source='trial baseline (approximate)',
         nr_frac=0.30),
    dict(label='Esgalhado resistant starch crossover (2020)',
         n=26, n_source='trial N, Food Funct 2020;11:2617-2625',
         tau_mean=0.25, tau_source='reported IS reduction ~25%',
         tau_sd=0.12, base_mean=6.0, base_sd=3.0,
         base_source='trial baseline (approximate), hemodialysis',
         nr_frac=0.20),
    dict(label='Sirich fiber parallel (2014)',
         n=40, n_source='trial N, CJASN 2014;9:1603-1610',
         tau_mean=0.25, tau_source='reported IS reduction ~25%',
         tau_sd=0.13, base_mean=5.0, base_sd=2.8,
         base_source='trial baseline (approximate), hemodialysis',
         nr_frac=0.20),
    dict(label='2025 fiber meta-analysis, IS outcome',
         n=398, n_source='pooled N for the IS outcome (11 RCTs, 200 + 198)',
         tau_mean=0.25,
         tau_source='SMD -0.34 does not convert to a percentage; 0.25 assumed',
         tau_sd=0.14, base_mean=5.4, base_sd=3.6,
         base_source='Lin 2011 CKD stage 4 (assumed, not meta-reported)',
         nr_frac=0.18),
]


def build(spec, seed):
    """Cohort with the trial's baseline moments and effect distribution."""
    rng = np.random.default_rng(seed)
    n = spec['n']
    base = np.clip(rng.normal(spec['base_mean'], spec['base_sd'], n), 0.5, None)
    tau = np.clip(rng.normal(spec['tau_mean'], spec['tau_sd'], n), 0, 0.70)
    nr = rng.random(n) < spec['nr_frac']
    tau[nr] = np.clip(rng.normal(0.03, 0.02, nr.sum()), 0, 0.08)
    c = make_cohort(n, seed, spec['nr_frac'])
    c.base_is = base
    c.tau = tau
    return c


def banner(t):
    print('\n' + '=' * 100)
    print(t)
    print('=' * 100)


banner('R6. INPUT PARAMETERS AND THEIR PROVENANCE')
print(f'  {"study":<42}{"N":>6}  {"N provenance":<48}')
print('  ' + '-' * 98)
for s in COHORTS:
    print(f'  {s["label"]:<42}{s["n"]:>6}  {s["n_source"]:<48}')
print()
print(f'  {"study":<42}{"tau":>6}  {"effect provenance":<48}')
print('  ' + '-' * 98)
for s in COHORTS:
    print(f'  {s["label"]:<42}{s["tau_mean"]:>6.2f}  {s["tau_source"]:<48}')
print()
print(f'  {"study":<42}{"tau SD":>8}{"NR frac":>9}   status')
print('  ' + '-' * 98)
for s in COHORTS:
    print(f'  {s["label"]:<42}{s["tau_sd"]:>8.2f}{s["nr_frac"]:>9.2f}'
          f'   BOTH ASSUMED (no source)')

banner('R6. OPERATING CHARACTERISTICS (primary rule, CV = 0.25, '
       'mean of 50 replications, seeds 700-749)')
print('  Read as parameter sensitivity, not external validation: no')
print('  individual-patient data from any of these trials were available.\n')
print(f'  {"study":<42}{"N":>6}{"resp":>7}{"sens":>8}{"spec":>8}{"weak":>7}')
print('  ' + '-' * 78)
for s in COHORTS:
    sens, spec, weak, resp = [], [], [], []
    for sd in PUB_SEEDS:
        c = build(s, sd)
        r = run_protocol(c, cv=0.25, measurement_seed=sd, **PRIMARY)
        sens.append(r['sens']); spec.append(r['spec'])
        weak.append(r['detect']['weak']); resp.append(r['n_responders'])
    print(f'  {s["label"]:<42}{s["n"]:>6}{np.mean(resp):>7.0f}'
          f'{np.nanmean(sens)*100:>7.1f}%{np.nanmean(spec)*100:>7.1f}%'
          f'{np.nanmean(weak)*100:>6.0f}%')

print('\n  The three AST-120 rows share one assumed effect size because only')
print('  Schulman 2006 measured indoxyl sulfate. Their differing operating')
print('  characteristics reflect cohort size alone, which is the point: an')
print('  effect attributed to EPPIC would be an attribution error, not a')
print('  finding.')
print('\n  Reproduce with: python3 validation/published_cohorts.py')
