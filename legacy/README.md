# Archived work, outside the scope of the paper

Nothing in this directory is run by `validation/reproduce_all.py`, and no number
in the manuscript or its supplement depends on any of it. It is kept as a record
of how the project got to its current design, and it is not maintained.

**Do not read these files as the current implementation.** Several of them state
conclusions the paper reverses. Each file carries a header naming which of its
statements is superseded; the ones that matter are collected here.

| File | What it says that the paper contradicts |
|---|---|
| `protocol/clinical_protocol.py` | Sets `CV_TARGET = 0.15` and calls it achievable; derives thresholds from a zero-effect null with no alpha-spending; prints treatment actions including deprescribing |
| `protocol/nof1_weak_rescue.py` | Treats a reduction of CV from 0.22 to 0.15 as an available design lever |
| `protocol/nof1_virtual_cohort.py` | Prints deprescribing as an outcome of classification |
| `sensitivity_analysis.py` | Contains an earlier feasibility argument for CV 0.15 |
| `novelty_and_gap_analysis.py` | Reports non-responders as correctly deprescribed |
| `simulations/redteam_loop.py` | Scores a quarterly deprescribing review as a design option |
| `simulations/why_experts_failed_revalidation.py` | Argues for deprescribing as a clinical objective |

What the paper concludes instead:

- A total CV of 0.15 is **not** reachable. Duplicate assays divide only the
  analytical variance, which contributes at most 0.004 to the total. The
  model-implied range under the assumed variance decomposition is 0.22 to 0.30.
- The decision rule tests H0: τ ≤ 0.10, matching the null to the definition of a
  responder, with alpha-spending boundaries calibrated by simulation. Testing
  against a zero-effect null while labelling responders at 10% is the mismatch
  the paper exists to quantify.
- No clinical decision rule was validated. A classification records what a
  simulated protocol demonstrated. It is not a recommendation to start, continue
  or stop any treatment.

See `../LIMITATIONS.md` sections (d) and (g) for the full withdrawals.

The remaining files under `simulations/` are exploratory work on toxin models and
intervention engineering. They informed no number in the paper.
