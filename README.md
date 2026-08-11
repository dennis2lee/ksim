# ksim: In-Silico N-of-1 Trial Design for Uremic Toxin Clearance in CKD

> **This is a research simulation framework and a candidate trial design, not a
> clinical decision tool.** Every operating characteristic here describes how a
> protocol behaves under a model, not in a clinic. A responder or non-responder
> call is a statement about what a simulated protocol demonstrated; it is not a
> recommendation to start, continue or stop any treatment. Do not use this code
> or these thresholds for patient care.

Simulation framework for designing and validating personalized n-of-1 crossover
trials that detect individual indoxyl sulfate (IS) reduction from gut-clearance
interventions in chronic kidney disease (CKD stage 3b–4).

**Paper**: Lee P, Lee TJ. "Reachable Measurement Precision and the Cost of a
Mismatched Null in N-of-1 Biomarker Trials: An In Silico Study in Chronic Kidney
Disease." (2026, under review)

## The Problem

Gut-targeted interventions (dietary fiber, probiotics, oral sorbents) reduce
indoxyl sulfate in CKD — a 2025 meta-analysis of 11 RCTs confirms a pooled
effect (SMD −0.34). But these trials report only population averages. A
clinician cannot tell whether a specific patient responds or not. The
within-person biological variability of IS (CV 35.9%, Pretorius 2013) makes single pre-post
comparisons unreliable.

## What This Code Does

1. **Generates virtual CKD cohorts** from published clinical data (IS levels,
   eGFR distributions, treatment effect sizes, biological variability)
2. **Simulates an n-of-1 crossover protocol** (2-cycle AB design with adaptive
   second stage for borderline patients)
3. **Quantifies the classification-efficiency frontier**: the relationship between
   measurement noise (CV), the one-sided decision threshold (DT), and the number of patients
   correctly classified as responders or non-responders
4. **Tests robustness** under 8 distributional violations (log-normal IS, bimodal
   treatment effects, correlated noise, partial carryover, combined worst case)

The central finding: sampling standardization is a more efficient design lever
than extending the protocol. It buys 15.9 sensitivity points with no additional
weeks and no additional venipunctures, against 11.8 points for doubling the
number of crossover cycles at a cost of 22 additional weeks. That comparison is
made at CV values a standardization package can actually reach, which is 0.22 to
0.30 rather than the 0.15 assumed in earlier versions of this work.

## Repository Structure

```
ksim/
├── simulations/          14 scripts — toxin models and intervention engineering
│   ├── gut_clearance_model.py          Baseline IS/PCS/urea generation & clearance
│   ├── drug_dev_simulation.py          Monte Carlo (N=20K), adherence, dose-response
│   ├── iter1_safety_feasibility.py     Side effects, feasibility gates, dialysis delay
│   ├── iter2_creative_bypass.py        Replace infeasible arms with deployable alternatives
│   ├── iter3_competing_risk.py         74-year-old competing risk model
│   ├── iter4_what_moves_needle.py      Toxin route vs AKI/CV prevention
│   ├── iter5_toxin_mortality_coupling.py  f_tox parameter sweep
│   ├── devA_sachet_formulation.py      Optimal fiber blend (RS/inulin/acacia)
│   ├── devC_collision_timing.py        Sorbent–probiotic timing separation
│   ├── dev_ketoacid_protein.py         Protein target + ketoacid optimization
│   ├── dev_integration_validation.py   Naive vs engineered stack comparison
│   ├── final_outcome_and_redteam.py    AST-120 drug interaction analysis
│   ├── redteam_loop.py                 15-round adversarial refinement (R: 48→82)
│   └── why_experts_failed_revalidation.py  10-pass expert failure escape analysis
│
├── protocol/             3 scripts — n-of-1 design and candidate classification rule
│   ├── nof1_virtual_cohort.py          N=100 virtual cohort, crossover power analysis
│   ├── nof1_weak_rescue.py             Weak responder strategies (adaptive + CV reduction)
│   └── clinical_protocol.py            Step-by-step protocol + classify_patient() function
│
├── validation/           one generative model, plus the analyses that call it
│   ├── nof1_core.py                    THE model: cohort, measurement, protocol, decision rules
│   ├── reproduce_all.py                One command; regenerates every published number
│   ├── requirements.txt                Pinned environment
│   ├── reproduce_manuscript_numbers.py Cohort, Table 1, Table 2, inline numbers
│   ├── decision_rule_analysis.py       Decision rules, CV sweep, cost-ratio sensitivity
│   ├── sequential_error_validation.py  Held-out boundary validation, Monte Carlo error
│   ├── correlation_boundary_analysis.py Correlation-adjusted decision boundaries
│   ├── protocol_dependence_analysis.py Onset, washout, carryover, adherence
│   ├── robustness_experiments.py       Reference plus 8 distributional stress conditions
│   ├── crossover_order_analysis.py     Cycle order and period effects
│   ├── variance_components_analysis.py Variance components, attrition, adherence
│   ├── efficiency_analysis.py          Joint design efficiency on a common burden scale
│   ├── published_cohorts.py            Cohorts parameterized to published trials
│   ├── threshold_sensitivity.py        Responder-definition threshold sweep
│   ├── nonresponder_fraction_sensitivity.py
│   ├── threshold_approximation_check.py Critical-value accuracy vs the empirical null
│   ├── sensitivity_analysis.py         Tornado sweep across literature-calibrated params
│   ├── literature_recalibration.py     Before/after comparison with published data
│   └── novelty_and_gap_analysis.py     RCT vs n-of-1 gap quantification
│
├── LIMITATIONS.md        Known limitations and anticipated reviewer responses
├── FUTURE_WORK.md        Research roadmap (14 directions, prioritized)
└── README.md             This file
```

## Quick Start

**Requirements**: Python 3.8+, NumPy. Matplotlib needed only for figure generation.

```bash
pip install numpy matplotlib

# Reproduce every number in the manuscript and the supplement, in one command.
# Any Python warning is promoted to an error, so a clean run means a clean run.
# Writes a full transcript to validation/reproduce_all.log (about 4 minutes).
python validation/reproduce_all.py

# Or run a single analysis
python validation/decision_rule_analysis.py
python validation/variance_components_analysis.py
```

## Key Results

The primary decision rule tests H0: τ ≤ 0.10, so the statistical null matches
the 10% definition of a true responder. Numbers under the zero-effect null
(H0: τ ≤ 0), which earlier versions of this work reported as primary, are kept
as a clearly labelled secondary analysis of a different question.

| Metric | Primary rule | Zero-effect null | Source script |
|--------|--------------|------------------|---------------|
| Single-run sensitivity (CV 0.15) | 83.4% | 95.1% | `validation/reproduce_manuscript_numbers.py` |
| Single-run specificity (CV 0.15) | 99.2% | 85.7% | same |
| Overall type I error, held-out seeds | 4.98% ± 0.05% | — | `validation/sequential_error_validation.py` |
| Type I error at within-period correlation 0.5 | 11.81% | — | `validation/correlation_boundary_analysis.py` |
| Model-implied reachable CV | 0.22 to 0.30 | same | `validation/variance_components_analysis.py` |
| Sensitivity at the reachable CV 0.26 | 65.2% | — | `validation/efficiency_analysis.py` |
| Standardization gain, three variance splits | +21.9 / +15.9 / +10.0 points | — | same |
| Doubling cycles gain | +11.7 to +11.8 points, +22 weeks | — | same |
| Cost of a 2-week onset time constant | −19.5 points | — | `validation/protocol_dependence_analysis.py` |

Three claims made in earlier versions of this work do not survive and are
withdrawn:

- **There is no optimal CV.** Under the matched null, net correct
  classifications rise monotonically as CV falls (636 at CV 0.36 to 919 at
  CV 0.10), at every false-positive-to-false-negative cost ratio examined. The
  apparent optimum near CV 0.12 was an artifact of testing against a
  zero-effect null while labelling responders at 10%.
- **CV 0.15 is not reachable.** Duplicate assays divide only the analytical
  variance, contributing at most 0.004 to the total. Reaching 0.15 would
  require an irreducible biological CV of 0.112 to 0.138 against a reported
  within-person biological CV of 0.359.

- **Standardization is not always the better lever.** It beats doubling the
  cycles under the optimistic and intermediate variance splits and loses to it
  under the pessimistic one (+10.0 against +11.7 sensitivity points). Reporting
  the intermediate split alone concealed a sign change.

One robustness scenario was also found to be vacuous: the estimator is a
within-patient ratio, so the baseline cancels exactly and a log-normal baseline
leaves every classification bit-identical. It is now reported as an analytic
invariance, not an empirical result.

## Parameter Sources

| Parameter | Value | Literature source |
|-----------|-------|-------------------|
| IS biological CV | 0.25 (conservative modeling baseline) | Pretorius et al. 2013, *Clin Chim Acta*, reports 35.9% for total serum IS in healthy volunteers. We simulate from 0.25, which understates rather than inflates the gain from standardization. |
| IS baseline CKD4 | 5.4 ± 3.6 µg/mL | Lin et al. 2011 |
| Fiber IS reduction (pooled) | SMD −0.34 | Wathanavasin et al. 2025, *Toxins*, 11 RCTs, N=398 |
| CKD4 eGFR slope | −2.0 mL/min/1.73 m²/yr | CRIC Study / MDRD |
| Non-responder fraction | 18% (scenario assumption) | Not estimable from published data; see Table S5 sweep |
| Gut microbiome variability | SD 0.35 | Wu et al. 2011, *Science* |

## Clinical Protocol Summary

The protocol classifies individual CKD patients as IS responders or non-responders:

1. **Eligibility**: CKD 3b–4, stable eGFR, IS assay available
2. **Standardize**: Fasting AM draws and timed collection, reaching CV ≈ 0.26
3. **Stage 1**: 2 × 3 crossover (24 weeks, 12 visits), period order randomized
   independently within each cycle
4. **Classify** against H0: τ ≤ 10%, so the test matches the responder
   definition: obs_red > 35.1% → Responder; obs_red < 10% → Non-responder;
   in between → Borderline
5. **Stage 2**: Borderline patients get 1 extra cycle (wk 24–36), threshold 30.5%
6. **After classification**: The call records what the protocol demonstrated. It is not a recommendation to continue or stop treatment.

The classification rule is implemented in `protocol/clinical_protocol.py` as the
`classify_patient()` function — input IS measurements, output classification.

## License

Released under the [MIT License](LICENSE). Copyright (c) 2026 Pyeongwoo Lee and Timothy Juheon Lee.

This project is for research and educational purposes. Simulation parameters
are illustrative and calibrated to published data — not for clinical dosing
decisions without empirical validation.

## Citation

```
Lee P, Lee TJ. Reachable Measurement Precision and the Cost of a Mismatched Null
in N-of-1 Biomarker Trials: An In Silico Study in Chronic Kidney Disease. 2026.


## Author

Pyeongwoo Lee — Independent Researcher, Sunnyvale, CA, USA
Timothy Juheon Lee — Independent Researcher, Sunnyvale, CA, USA
Contact: dennis2.lee@gmail.com
