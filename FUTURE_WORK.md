# Future Research Roadmap

Based on the current project state: in-silico n-of-1 protocol design validated on N=1,000
virtual cohort, literature-calibrated parameters, CV-as-design-lever insight,
six-step clinical protocol with decision engine, and manuscript targeting CPT:PSP.

---

## 1. Validation Pathway — from simulation to real measurement

### 1A. Retrospective IS variability characterization (Priority: ★★★)

**Goal**: Confirm the CV=0.15 assumption by measuring actual within-person IS variability
under standardized conditions in CKD patients.

**Design**: Single-center observational study. Recruit 20–30 CKD stage 3b–4 patients.
Collect duplicate fasting AM serum IS draws at 3 visits over 2 weeks (no intervention).
Compute within-person CV under standardized versus unstandardized conditions.

**Why first**: The entire protocol's efficiency rests on achievable CV. If CV=0.15 is
not reachable, the DT stays at 24% and the design reverts to a 4-cycle requirement.
Conversely, if CV turns out to be 0.12, the protocol becomes even more powerful.

**What's needed**: IRB approval, LC-MS/MS IS assay access (or partnership with a reference
lab), 20–30 consenting CKD patients, ~2 months of data collection.

**Difficulty**: Low (observational, no intervention, standard blood draws).

**Output**: Short report or letter confirming achievable CV; feeds directly into
manuscript revision or follow-up paper.

### 1B. Single-patient proof-of-concept n-of-1 (Priority: ★★★)

**Goal**: Run the full 2×3 protocol on one willing CKD patient to demonstrate
operational feasibility — can the timeline, washouts, measurement schedule, and
fiber titration be followed in practice?

**Design**: One CKD 3b–4 patient, the 24-week protocol as specified in
`08_clinical_protocol/clinical_protocol.py`. Measure IS (duplicate, fasting AM)
at all 12 visits. Apply the decision rule. Report observed reduction,
classification, adherence, tolerability, and protocol deviations.

**Why early**: Demonstrates the protocol is executable, not just simulable.
Even a single-patient report — if it shows the measurement schedule is adherent
and the IS signal is interpretable — carries weight with reviewers and clinicians.

**What's needed**: One motivated patient, nephrologist collaborator, IS assay,
24 weeks, minimal funding (blood draws only).

**Difficulty**: Low-moderate (regulatory: may qualify as quality improvement, not
formal trial, depending on jurisdiction; fiber/probiotic are OTC).

**Output**: Case report or n-of-1 methods paper (e.g., BMJ Open, J Personalized Medicine).

### 1C. Retrospective RCT individual-level reanalysis (Priority: ★★)

**Goal**: Obtain individual-patient data (IPD) from a published fiber or AST-120
trial and apply the classification algorithm retrospectively. Confirm that the
simulated operating characteristics match what we see in real data.

**Design**: Contact authors of Sirich 2014, Rossi SYNERGY 2016, or Esgalhado 2020
for IPD. Apply our decision threshold to each patient's observed IS change.
Compare the classification rate with our simulation predictions.

**What's needed**: Data-sharing agreement with original investigators, or access to
a trial repository. Statistical reanalysis.

**Difficulty**: Moderate (depends on data availability; IPD sharing is improving
but not guaranteed).

**Output**: Validation paper demonstrating that simulation predictions hold on real
data. Strong evidence for or against the non-responder fraction estimate.

---

## 2. Methodological Extensions

### 2A. Multi-endpoint: IS + PCS + urea composite (Priority: ★★)

**Goal**: Extend the protocol to detect individual effects on PCS and urea
alongside IS. Some patients may respond on PCS but not IS (fiber is more
PCS-specific in some studies). A composite endpoint could rescue patients who
are IS non-responders but PCS responders.

**Design**: Simulate 2×3 crossover measuring IS, PCS, and urea simultaneously.
Define responder as significant reduction in any one toxin. Compute the gain in
sensitivity from a composite versus IS-only.

**What's needed**: Pretorius 2013 reports a within-person biological CV of 50.5%
for total pCS, higher than the 35.9% for total IS, so a composite endpoint would
carry the noisier marker. Also needed: the correlation structure between IS and
PCS within patients.

**Difficulty**: Low (simulation extension, same framework).

**Output**: Extension paper or supplementary analysis for the primary manuscript.

### 2B. Bayesian adaptive design with formal alpha-spending (Priority: ★★)

**Goal**: Replace the current threshold-based two-stage design with a Bayesian
sequential framework that formally controls the false-positive rate and allows
for continuous monitoring rather than fixed interim analysis.

**Design**: Implement a group-sequential or Bayesian predictive probability
design. At each measurement, update the posterior on the patient's true τ.
Stop for efficacy or futility when posterior probability exceeds pre-specified
boundaries. Compute operating characteristics via simulation.

**What's needed**: Statistical programming (Stan or PyMC), familiarity with
alpha-spending functions (O'Brien-Fleming, Pocock).

**Difficulty**: Moderate (statistical methodology).

**Output**: Methodological paper on Bayesian n-of-1 design for biomarkers with
high biological variability. Strong fit for CPT:PSP or Statistics in Medicine.

### 2C. Mixture-model empirical Bayes for individual effect estimation (Priority: ★)

**Goal**: Improve individual τ estimates beyond simple thresholding by fitting
a two-component mixture model (responders + non-responders) and computing
posterior probabilities of membership for each patient.

**Design**: After running the protocol, fit a mixture of N(μ_resp, σ²) and
N(0, σ²_null) to the observed reductions. Each patient gets a posterior
probability of being a responder. Compare classification accuracy with the
threshold-based approach.

**What's needed**: EM algorithm implementation, simulation validation.

**Difficulty**: Moderate.

**Output**: Methodological extension; may be incorporated into the primary
paper as a sensitivity analysis or published separately.

### 2D. Measurement noise model refinement (Priority: ★)

**Goal**: Replace the i.i.d. Gaussian noise model with a more realistic model
incorporating (a) serial correlation within a treatment period, (b) heteroscedasticity
(variance that depends on IS level), and (c) assay-specific error structure.

**What's needed**: Access to repeated IS measurements from the same patients to
estimate autocorrelation and variance structure. Possibly from Study 1A above.

**Difficulty**: Low-moderate.

**Output**: Supplementary analysis showing robustness (or sensitivity) to
noise model assumptions.

---

## 3. Scope Expansion

### 3A. Other CKD populations (Priority: ★★)

**Goal**: Validate protocol performance in diabetic CKD, younger patients,
CKD stage 5 pre-dialysis, and post-transplant patients. Each population has
different baseline IS levels, eGFR trajectories, and microbiome compositions.

**Design**: Generate population-specific virtual cohorts from published data
for each subgroup. Re-run the full protocol simulation. Identify populations
where the protocol needs modification (e.g., different decision thresholds,
longer washouts).

**What's needed**: Published IS data stratified by CKD etiology and demographics.

**Difficulty**: Low (simulation, already have the framework).

**Output**: Extension paper or supplementary table for the primary manuscript.

### 3B. Hemodialysis patients — different kinetics (Priority: ★)

**Goal**: Adapt the protocol for hemodialysis patients, where IS kinetics are
dominated by dialysis clearance rather than residual renal function. The
measurement timing relative to dialysis sessions becomes critical.

**What's needed**: Published IS kinetics during and between dialysis sessions;
rebound curves post-dialysis.

**Difficulty**: Moderate (different pharmacokinetic model needed).

### 3C. Other interventions beyond fiber/probiotic (Priority: ★★)

**Goal**: Apply the framework to pharmacological interventions: tryptophanase
inhibitors (if they reach clinical development), targeted sorbents, engineered
bacteria, or dietary protein modulation. Each has different onset kinetics,
dose-response, and washout characteristics.

**What's needed**: Intervention-specific PK/PD parameters from preclinical or
early clinical data.

**Difficulty**: Low per intervention (parameter swap in existing framework).

### 3D. Digital decision support tool (Priority: ★★)

**Goal**: Package the `classify_patient()` function from `08_clinical_protocol/`
into a standalone web application or spreadsheet template that a clinician can
use at the bedside. Input: 12 IS values. Output: classification + confidence.

**What's needed**: Web development (Flask/Streamlit) or Excel VBA; UX design
for clinical users; potentially CE marking or FDA 510(k) if marketed as a
medical device (may be exempt as clinical decision support).

**Difficulty**: Low for prototype, moderate for regulatory.

**Output**: Software tool paper; app or web link.

---

## 4. Clinical Translation

### 4A. Prospective multi-patient n-of-1 series (Priority: ★★★)

**Goal**: Run the protocol on 10–20 patients as a prospective n-of-1 series.
Aggregate individual results to estimate the population-level responder fraction
and validate the simulated operating characteristics.

**Design**: Multi-center or single-center. Each patient runs the full 2×3
protocol (with optional Stage 2). Pre-register on ClinicalTrials.gov.
Primary outcome: proportion of patients classified as responders.
Secondary: observed IS reduction distribution, adherence, tolerability,
comparison with simulation predictions.

**What's needed**: IRB approval, 2–3 sites, IS assay standardization across
sites (or central lab), 10–20 patients, ~12 months enrollment + 6 months
follow-up. Funding: modest (dietary intervention, standard labs).

**Difficulty**: Moderate (multi-patient coordination, but low-risk intervention).

**Output**: Clinical validation paper in Clinical Kidney Journal or CJASN.
This is the study that converts simulation into evidence.

### 4B. Regulatory and guideline pathway (Priority: ★)

**Goal**: Map the regulatory requirements for (a) using n-of-1 as a prescribing
tool (not a marketing authorization tool), and (b) the decision support software
as a clinical tool.

**Key questions**: Does n-of-1 evidence qualify under precision medicine
frameworks? Can it support individual prescribing decisions without a
population-level RCT? How do KDIGO or ERA guidelines accommodate n-of-1?

**What's needed**: Regulatory affairs consultation; literature review of n-of-1
in guideline frameworks; engagement with KDIGO working groups.

**Difficulty**: High (regulatory/policy, slow timeline).

### 4C. Health economic evaluation (Priority: ★)

**Goal**: Quantify the cost-effectiveness of the n-of-1 protocol versus
treat-all or treat-none strategies. The protocol costs ~12 visits and 24–36
blood draws; the savings are in averted burden for non-responders and
potentially in delayed dialysis for responders (if f_tox > 0).

**What's needed**: Health economic modeling (Markov or microsimulation),
cost data for IS assays, clinic visits, fiber/probiotic regimens, and
dialysis.

**Difficulty**: Moderate.

**Output**: Cost-effectiveness paper.

---

## Summary Table

| ID | Direction | Priority | Difficulty | Timeline | Output |
|----|-----------|----------|------------|----------|--------|
| **1A** | **IS variability characterization** | **★★★** | **Low** | **3–6 mo** | **Letter/short report** |
| **1B** | **Single-patient proof-of-concept** | **★★★** | **Low-mod** | **6–9 mo** | **Case report** |
| 1C | Retrospective IPD reanalysis | ★★ | Moderate | 6–12 mo | Validation paper |
| 2A | Multi-endpoint (IS+PCS+urea) | ★★ | Low | 1–2 mo | Extension/suppl |
| **2B** | **Bayesian adaptive design** | **★★** | **Moderate** | **3–4 mo** | **Methods paper (CPT:PSP)** |
| 2C | Mixture-model EB | ★ | Moderate | 2–3 mo | Methods extension |
| 2D | Noise model refinement | ★ | Low-mod | 1–2 mo | Supplementary |
| 3A | Other CKD populations | ★★ | Low | 1–2 mo | Extension paper |
| 3B | Hemodialysis adaptation | ★ | Moderate | 3–6 mo | New paper |
| 3C | Other interventions | ★★ | Low | 1–2 mo/each | Extension |
| 3D | Digital decision tool | ★★ | Low-mod | 2–4 mo | Software paper |
| **4A** | **Prospective n-of-1 series** | **★★★** | **Moderate** | **12–18 mo** | **Clinical validation (CKJ/CJASN)** |
| 4B | Regulatory pathway | ★ | High | 6–12 mo | Policy paper |
| 4C | Health economics | ★ | Moderate | 3–6 mo | HE paper |

---

## Recommended Sequence (Critical Path)

```
NOW (in-silico, current paper)
  │
  ├── 1A: IS variability study (3-6 mo, confirms CV assumption)
  │     │
  │     └── feeds back to refine protocol parameters
  │
  ├── 2B: Bayesian adaptive design (3-4 mo, parallel with 1A)
  │     │
  │     └── second methods paper for CPT:PSP
  │
  └── 1B: Single-patient proof-of-concept (6-9 mo, can start during 1A)
        │
        └── 4A: Prospective 10-20 patient series (12-18 mo)
              │
              └── Clinical validation → guideline consideration
```

The first two steps (1A and 1B) can proceed with minimal funding and regulatory
burden. They produce the empirical evidence that the current purely in-silico
framework lacks. Everything downstream — the multi-patient series, the digital
tool, the regulatory engagement — becomes credible once a real IS measurement
confirms the simulation's prediction.
