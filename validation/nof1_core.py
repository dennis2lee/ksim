"""
nof1_core.py -- the single generative model and protocol runner used by every
table and figure script in this repository.

Before this module existed each analysis script carried its own copy of the
cohort generator, and the copies had drifted: `evsi_analysis.py` used seed 800
and omitted the eGFR-slope draw, while `robustness_experiments.py` used seed 777
and included it, so two tables that both claimed to describe "the reference
model" were in fact built on different virtual cohorts. Everything now routes
through `make_cohort` and `run_protocol` here.

Two seeds are kept separate and must be passed separately:

  cohort_seed       fixes the virtual patients (baseline IS, tau, eGFR, slope).
  measurement_seed  fixes the measurement noise drawn during a protocol run.

Reusing one cohort across analyses is therefore explicit: pass the same
`cohort_seed` and vary `measurement_seed`. The reference configuration
(REFERENCE below) reproduces the manuscript's Table 2 confusion matrix exactly.

Decision rule
-------------
The protocol estimates a within-patient reduction

    obs = (mean_A - mean_B) / mean_A

whose standard error is approximately cv * sqrt(2/n). A patient is called a
responder when obs exceeds a critical value

    crit = null_margin + z(alpha) * cv * sqrt(2/n)

`null_margin` is the effect size under the null hypothesis. Setting it to 0
tests H0: tau <= 0 ("any detectable reduction"); setting it to the clinical
responder definition theta tests H0: tau <= theta, which is the rule that
matches how true responders are labelled. Both are supported so the two can be
compared directly rather than conflated.
"""

import numpy as np
from statistics import NormalDist

_N = NormalDist()          # stdlib; avoids a scipy dependency for two functions
_z = _N.inv_cdf            # standard normal quantile
_phi = _N.cdf              # standard normal CDF

# ---------------------------------------------------------------------------
# Reference configuration
# ---------------------------------------------------------------------------

REFERENCE = dict(
    n_pat=1000,
    cohort_seed=777,
    measurement_seed=777,
    cv=0.15,
    n_s1=6,             # measurements per arm, Stage 1 (2 cycles x 3)
    n_s2=3,             # additional measurements per arm, Stage 2
    wk_treat=4,
    wk_wash=2,
    theta=0.10,         # clinical responder definition: tau >= theta
    null_margin=0.0,    # H0: tau <= null_margin
    alpha1=0.05,        # one-sided, Stage 1
    alpha2=0.05,        # one-sided, Stage 2
    nr_frac=0.18,
    order='fixed',
)

# One crossover cycle occupies 12 weeks: a 4-week period, a 2-week washout, a
# second 4-week period, a second 2-week washout. Measurements are taken in the
# final `meas_per_period` weeks of each period, so that the microbiome shift has
# time to establish before the first draw.
CYCLE_WEEKS = 12
PERIOD_WEEKS = 4
WASHOUT_WEEKS = 2


def cycle_measurement_weeks(cycle, meas_per_period):
    """Measurement weeks for the two periods of one 0-indexed cycle."""
    base = cycle * CYCLE_WEEKS
    first_start = PERIOD_WEEKS - meas_per_period + 1
    p1 = tuple(base + first_start + i for i in range(meas_per_period))
    offset = PERIOD_WEEKS + WASHOUT_WEEKS
    p2 = tuple(w + offset for w in p1)
    return p1, p2


def make_schedule(order, n_cycles=2, meas_per_period=3, rng=None,
                  start_cycle=0):
    """Build the (arm, measurement weeks) schedule for a run.

    order : 'fixed'            control first in every cycle (original protocol)
            'reverse'          treatment first in every cycle
            'counterbalanced'  alternating AB, BA, ... (ABBA for two cycles),
                               which centres both arms on the same mean week
            'randomized'       AB or BA drawn independently for each cycle

    With order='fixed', n_cycles=2 and meas_per_period=3 this returns the
    original manuscript schedule: A(2,3,4) B(8,9,10) A(14,15,16) B(20,21,22).
    """
    sched = []
    for c in range(start_cycle, start_cycle + n_cycles):
        p1, p2 = cycle_measurement_weeks(c, meas_per_period)
        if order == 'fixed':
            arms = ('A', 'B')
        elif order == 'reverse':
            arms = ('B', 'A')
        elif order == 'counterbalanced':
            arms = ('A', 'B') if c % 2 == 0 else ('B', 'A')
        elif order == 'randomized':
            arms = ('A', 'B') if rng.random() < 0.5 else ('B', 'A')
        else:
            raise ValueError(f'unknown order {order!r}')
        sched.append((arms[0], p1))
        sched.append((arms[1], p2))
    return sched


# ---------------------------------------------------------------------------
# Cohort generation
# ---------------------------------------------------------------------------

class Cohort:
    """A virtual cohort. Attributes are per-patient arrays of length n_pat."""

    __slots__ = ('base_is', 'tau', 'egfr', 'slope', 'gut', 'n_pat',
                 'seed', 'nr_mask', 'baseline_shape', 'tau_shape')

    def __init__(self, **kw):
        for k, v in kw.items():
            setattr(self, k, v)

    def true_responder(self, theta):
        return self.tau >= theta

    def strata(self, theta=0.10):
        """Weak / moderate / strong responder masks, as used in the manuscript."""
        r = self.true_responder(theta)
        return dict(
            non_responder=~r,
            weak=r & (self.tau < 0.20),
            moderate=r & (self.tau >= 0.20) & (self.tau < 0.35),
            strong=self.tau >= 0.35,
        )


def make_cohort(n_pat=1000, seed=777, nr_frac=0.18,
                baseline_shape='gaussian', tau_shape='gaussian'):
    """Generate a virtual CKD stage 3b-4 cohort.

    The draw order (eGFR, gut capacity, eGFR slope, baseline IS, tau,
    non-responder mask) is fixed and must not be reordered: it defines the
    random stream, and the manuscript's reported numbers are tied to it.

    baseline_shape : 'gaussian'   affine-standardized eGFR-dependent draw
                     'lognormal'  right-skewed alternative (robustness A)
    tau_shape      : 'gaussian'   clipped normal around 0.30 * gut capacity
                     'bimodal'    sharp responder / non-responder split (B)
    """
    rng = np.random.default_rng(seed)

    egfr = np.clip(rng.normal(22, 6, n_pat), 10, 35)
    gut = np.clip(rng.normal(1.0, 0.35, n_pat), 0.3, 1.8)
    slope = np.clip(rng.normal(-2.0, 0.8, n_pat), -5.0, -0.5)

    if baseline_shape == 'gaussian':
        raw = (25.0 / egfr) ** 1.2 * np.clip(rng.normal(1.0, 0.30, n_pat), 0.3, 2.5)
        # Affine-standardize to the CKD stage 4 cohort mean 5.4 / SD 3.6 ug/mL
        # (Lin 2011), then floor at 0.5 to keep concentrations positive. The
        # floor perturbs the moments slightly; realized values are reported by
        # `baseline_moments` rather than asserted to be exactly 5.4 / 3.6.
        base_is = np.clip(5.4 + 3.6 * (raw - raw.mean()) / raw.std(), 0.5, None)
    elif baseline_shape == 'lognormal':
        ln_sigma = np.sqrt(np.log(1 + (3.6 / 5.4) ** 2))
        ln_mu = np.log(5.4) - ln_sigma ** 2 / 2
        base_is = (np.exp(rng.normal(ln_mu, ln_sigma, n_pat))
                   * (25.0 / egfr) ** 1.2 / ((25.0 / 22) ** 1.2))
    else:
        raise ValueError(f'unknown baseline_shape {baseline_shape!r}')

    if tau_shape == 'gaussian':
        tau = np.clip(rng.normal(0.30 * gut, 0.14), 0, 0.70)
        nr = rng.random(n_pat) < nr_frac
        tau[nr] = np.clip(rng.normal(0.03, 0.02, nr.sum()), 0, 0.08)
    elif tau_shape == 'bimodal':
        nr = rng.random(n_pat) < nr_frac
        tau = np.where(nr,
                       np.clip(rng.normal(0.02, 0.015, n_pat), 0, 0.06),
                       np.clip(rng.normal(0.38, 0.08, n_pat), 0, 0.70))
    else:
        raise ValueError(f'unknown tau_shape {tau_shape!r}')

    return Cohort(base_is=base_is, tau=tau, egfr=egfr, slope=slope, gut=gut,
                  n_pat=n_pat, seed=seed, nr_mask=nr,
                  baseline_shape=baseline_shape, tau_shape=tau_shape)


def baseline_moments(cohort):
    """Realized mean and SD of baseline IS, after flooring.

    The standardization step imposes 5.4 / 3.6 on the pre-floor values; the
    0.5 floor then shifts them. Report these, do not claim the targets are met
    exactly.
    """
    return float(cohort.base_is.mean()), float(cohort.base_is.std())


# ---------------------------------------------------------------------------
# Measurement model
# ---------------------------------------------------------------------------

def total_cv(cv_bio=None, cv_pre=None, cv_analytical=None, n_replicates=1,
             cv=None):
    """Combine independent variance components into one total CV.

    Replicate assays on the same specimen average down only the analytical
    component; day-to-day biological variation and pre-analytical variation
    are properties of the specimen and are untouched by re-assaying it.

        cv_total^2 = cv_bio^2 + cv_pre^2 + cv_analytical^2 / n_replicates

    Passing `cv` instead returns it unchanged, so callers that work with a
    single lumped error term keep working.
    """
    if cv is not None:
        return cv
    return float(np.sqrt(cv_bio ** 2 + cv_pre ** 2
                         + cv_analytical ** 2 / n_replicates))


def _draw_period(rng, base, effect, cv, n, rho, weeks, drift, egfr0, slope,
                 period_effect, p_missing):
    """Draw n measurements for one treatment period of one patient."""
    if rho > 0:
        raw = np.empty(n)
        raw[0] = rng.normal(0, cv)
        for j in range(1, n):
            raw[j] = rho * raw[j - 1] + np.sqrt(1 - rho ** 2) * rng.normal(0, cv)
    else:
        raw = rng.normal(0, cv, n)

    meas = base * (1 - effect) * (1 + raw)

    weeks = np.asarray(weeks, dtype=float)
    if drift:
        ew = egfr0 + slope * (weeks / 52.0)
        meas = meas * (egfr0 / ew) ** 1.2
    if period_effect is not None:
        meas = meas * period_effect(weeks)

    if p_missing > 0:
        keep = rng.random(n) >= p_missing
        if keep.any():
            meas = meas[keep]
        else:
            meas = meas[:1]          # never lose a period entirely
    return meas


def linear_period_effect(pct_per_week):
    """Multiplicative linear time trend, e.g. 0.005 = +0.5% per week."""
    return lambda weeks: 1.0 + pct_per_week * weeks


def seasonal_period_effect(amplitude, period_weeks=52.0):
    """Non-linear (sinusoidal) time trend, e.g. seasonal dietary variation."""
    return lambda weeks: 1.0 + amplitude * np.sin(2 * np.pi * weeks / period_weeks)


# ---------------------------------------------------------------------------
# Estimators
# ---------------------------------------------------------------------------

def ratio_estimate(a_vals, b_vals):
    """Simple arm-mean ratio, the manuscript's original estimator."""
    ma = a_vals.mean()
    return (ma - b_vals.mean()) / ma if ma > 0 else 0.0


def ols_period_estimate(a_vals, a_weeks, b_vals, b_weeks):
    """Reduction estimated with a per-patient linear time term.

    Regresses log(measurement) on a treatment indicator and week, then converts
    the treatment coefficient back to a fractional reduction. This separates a
    treatment effect from a monotonic drift that the arm-mean ratio absorbs.
    Falls back to the ratio estimator if the design matrix is rank-deficient.
    """
    y = np.log(np.concatenate([a_vals, b_vals]))
    trt = np.concatenate([np.zeros(len(a_vals)), np.ones(len(b_vals))])
    wk = np.concatenate([np.asarray(a_weeks, float), np.asarray(b_weeks, float)])
    X = np.column_stack([np.ones_like(y), trt, wk])
    if np.linalg.matrix_rank(X) < X.shape[1]:
        return ratio_estimate(a_vals, b_vals)
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return float(1 - np.exp(beta[1]))


# ---------------------------------------------------------------------------
# Critical values
# ---------------------------------------------------------------------------

def critical_value(cv, n, alpha=0.05, null_margin=0.0):
    """One-sided critical difference for H0: tau <= null_margin.

    This is a critical difference, not a power-based minimum detectable effect:
    there is no z(1-beta) term. At null_margin = 0, cv = 0.15 and n = 6 it
    returns 0.1425.
    """
    return null_margin + _z(1 - alpha) * cv * np.sqrt(2.0 / n)


def obf_alphas(alpha_total, info_fractions):
    """O'Brien-Fleming style one-sided alpha-spending.

    Returns the cumulative alpha spent at each look under the spending function
    alpha(t) = 2 * (1 - Phi(z(alpha_total/2) / sqrt(t))), the one-sided analogue
    of the standard OBF boundary. These are nominal per-look levels; the
    achieved overall error under this protocol's ratio estimator is calibrated
    by simulation in `calibrate_boundaries`, because the ratio statistic is not
    exactly normal and the futility rule also removes patients at Stage 1.
    """
    z = _z(1 - alpha_total)
    return [float(1 - _phi(z / np.sqrt(t))) for t in info_fractions]


# ---------------------------------------------------------------------------
# Protocol
# ---------------------------------------------------------------------------

def run_protocol(cohort, cv=0.15, measurement_seed=777,
                 n_cycles=2, meas_per_period=3, theta=0.10, null_margin=0.0,
                 alpha1=0.05, alpha2=0.05, order='fixed',
                 rho=0.0, carryover=0.0, drift=False,
                 period_effect=None, estimator='ratio',
                 p_missing=0.0, p_dropout=0.0,
                 cv_per_patient=None,
                 n_s1=None, n_s2=None):
    """Run the full two-stage adaptive protocol over a cohort.

    Stage 1 collects n_s1 measurements per arm over two crossover cycles and
    applies a three-way rule: call responder above the Stage 1 critical value,
    call non-responder below the null margin, otherwise continue to Stage 2.
    Stage 2 adds n_s2 measurements per arm and applies a single final rule to
    the pooled data.

    Returns a dict with per-patient classifications and summary operating
    characteristics. Dropouts are reported separately and excluded from
    sensitivity and specificity denominators rather than silently called
    non-responders.
    """
    rng = np.random.default_rng(measurement_seed)
    n = cohort.n_pat
    cvp = np.full(n, cv) if cv_per_patient is None else np.asarray(cv_per_patient)

    # n_s1 / n_s2 are derived from the design; they remain accepted as explicit
    # overrides so older call sites keep working.
    per_arm = meas_per_period
    if n_s1 is None:
        n_s1 = n_cycles * meas_per_period
    if n_s2 is None:
        n_s2 = meas_per_period

    static_sched = (None if order == 'randomized'
                    else make_schedule(order, n_cycles, meas_per_period))
    stage2_sched = None if order == 'randomized' else make_schedule(
        order, 1, meas_per_period, start_cycle=n_cycles)

    cls = np.full(n, 'N', dtype='U1')
    went_s2 = np.zeros(n, bool)
    dropped = np.zeros(n, bool)
    obs1_all = np.full(n, np.nan)
    obs_final = np.full(n, np.nan)

    for p in range(n):
        b = cohort.base_is[p]
        t = cohort.tau[p]
        c = cvp[p]
        sched = (make_schedule('randomized', n_cycles, meas_per_period, rng)
                 if order == 'randomized' else static_sched)

        a_vals, a_weeks, b_vals, b_weeks = [], [], [], []
        prev_on = False
        for arm, weeks in sched:
            on = (arm == 'B')
            effect = t if on else (t * carryover if prev_on else 0.0)
            vals = _draw_period(rng, b, effect, c, per_arm, rho, weeks, drift,
                                cohort.egfr[p], cohort.slope[p], period_effect,
                                p_missing)
            wk = np.asarray(weeks, float)[:len(vals)]
            (b_vals if on else a_vals).append(vals)
            (b_weeks if on else a_weeks).append(wk)
            prev_on = on

        A = np.concatenate(a_vals); AW = np.concatenate(a_weeks)
        B = np.concatenate(b_vals); BW = np.concatenate(b_weeks)

        if p_dropout > 0 and rng.random() < p_dropout:
            dropped[p] = True
            continue

        obs1 = (ratio_estimate(A, B) if estimator == 'ratio'
                else ols_period_estimate(A, AW, B, BW))
        obs1_all[p] = obs1

        crit1 = critical_value(c, n_s1, alpha1, null_margin)
        if obs1 > crit1:
            cls[p] = 'R'
            obs_final[p] = obs1
            continue
        if obs1 < null_margin:
            cls[p] = 'N'
            obs_final[p] = obs1
            continue

        # Stage 2
        went_s2[p] = True
        s2 = (make_schedule('randomized', 1, meas_per_period, rng,
                            start_cycle=n_cycles)
              if order == 'randomized' else stage2_sched)
        prev_on = (sched[-1][0] == 'B')
        for arm, weeks in s2:
            on = (arm == 'B')
            effect = t if on else (t * carryover if prev_on else 0.0)
            prev_on = on
            vals = _draw_period(rng, b, effect, c, n_s2, rho, weeks, drift,
                                cohort.egfr[p], cohort.slope[p], period_effect,
                                p_missing)
            wk = np.asarray(weeks, float)[:len(vals)]
            if on:
                B = np.concatenate([B, vals]); BW = np.concatenate([BW, wk])
            else:
                A = np.concatenate([A, vals]); AW = np.concatenate([AW, wk])

        obs2 = (ratio_estimate(A, B) if estimator == 'ratio'
                else ols_period_estimate(A, AW, B, BW))
        obs_final[p] = obs2
        crit2 = critical_value(c, n_s1 + n_s2, alpha2, null_margin)
        cls[p] = 'R' if obs2 > crit2 else 'N'

    return summarize(cohort, cls, went_s2, dropped, theta,
                     obs1=obs1_all, obs_final=obs_final)


def summarize(cohort, cls, went_s2, dropped, theta=0.10, obs1=None,
              obs_final=None):
    """Operating characteristics, with dropouts excluded from the denominators."""
    analysed = ~dropped
    resp = cohort.true_responder(theta) & analysed
    nonr = (~cohort.true_responder(theta)) & analysed
    called = (cls == 'R') & analysed

    tp = int((called & resp).sum())
    fn = int((~called & resp).sum())
    fp = int((called & nonr).sum())
    tn = int((~called & nonr).sum())

    strata = cohort.strata(theta)
    detect = {k: (float(((cls == 'R') & m & analysed).sum()
                        / max((m & analysed).sum(), 1)))
              for k, m in strata.items()}

    n_s2 = int((went_s2 & analysed).sum())
    n_an = int(analysed.sum())
    return dict(
        cls=cls, went_s2=went_s2, dropped=dropped,
        obs1=obs1, obs_final=obs_final,
        n_analysed=n_an, n_dropped=int(dropped.sum()),
        tp=tp, fn=fn, fp=fp, tn=tn,
        n_responders=tp + fn, n_nonresponders=fp + tn,
        sens=tp / (tp + fn) if tp + fn else float('nan'),
        spec=tn / (tn + fp) if tn + fp else float('nan'),
        fp_rate=fp / (fp + tn) if fp + tn else float('nan'),
        ncc=tp + tn,
        stage2_rate=n_s2 / n_an if n_an else float('nan'),
        mean_weeks=24 + 12 * (n_s2 / n_an) if n_an else float('nan'),
        detect=detect,
    )


def weighted_accuracy(res, fp_cost=1.0, fn_cost=1.0):
    """Accuracy with asymmetric misclassification costs.

    Net correct classifications (TP + TN) weights a false positive and a false
    negative equally, which is a strong and unstated assumption. This returns
    the expected loss per patient so conclusions can be checked across cost
    ratios.
    """
    n = res['tp'] + res['fn'] + res['fp'] + res['tn']
    return (fp_cost * res['fp'] + fn_cost * res['fn']) / n if n else float('nan')


def calibrate_boundaries(cohort, cv, null_margin, alpha_total=0.05,
                         n_s1=6, n_s2=3, seeds=range(900, 950), **kw):
    """Find per-look alphas whose achieved overall type I error is alpha_total.

    Runs the protocol on a null cohort (every tau set to the null margin) and
    scales the O'Brien-Fleming nominal levels by a common factor until the
    observed false-positive rate matches the target. Calibrating by simulation
    rather than asserting the asymptotic boundary is necessary here because the
    statistic is a ratio and the Stage 1 futility rule removes patients before
    the second look.
    """
    null_cohort = make_cohort(cohort.n_pat, cohort.seed)
    null_cohort.tau = np.full(cohort.n_pat, null_margin)

    info = [n_s1 / (n_s1 + n_s2), 1.0]
    a1_nom, a2_nom = obf_alphas(alpha_total, info)

    def achieved(scale):
        rates = []
        for s in seeds:
            r = run_protocol(null_cohort, cv=cv, measurement_seed=s,
                             n_s1=n_s1, n_s2=n_s2, theta=null_margin + 1e-9,
                             null_margin=null_margin,
                             alpha1=min(a1_nom * scale, 0.5),
                             alpha2=min(a2_nom * scale, 0.5), **kw)
            an = ~r['dropped']
            rates.append(((r['cls'] == 'R') & an).sum() / max(an.sum(), 1))
        return float(np.mean(rates))

    lo, hi = 0.01, 20.0
    for _ in range(40):
        mid = np.sqrt(lo * hi)
        if achieved(mid) > alpha_total:
            hi = mid
        else:
            lo = mid
    scale = np.sqrt(lo * hi)
    return dict(alpha1=min(a1_nom * scale, 0.5), alpha2=min(a2_nom * scale, 0.5),
                nominal=(a1_nom, a2_nom), scale=scale,
                achieved=achieved(scale))
