"""
Regenerate every number in the manuscript and the supplement with one command:

    python3 validation/reproduce_all.py

Runs each analysis script in dependency order and reports pass/fail plus
runtime. Any Python warning is promoted to an error, so a clean run means a
clean run. Writes the combined transcript to reproduce_all.log next to this
file, which is the log referenced in the Data availability statement.

Requires only Python 3.8+ and numpy; see requirements.txt for the pinned
version used to produce the submitted figures.
"""

import io
import contextlib
import sys
import time
import runpy
import warnings
from pathlib import Path

HERE = Path(__file__).resolve().parent

SCRIPTS = [
    ('reproduce_manuscript_numbers.py', 'Cohort, Table 1, Table 2, inline numbers'),
    ('decision_rule_analysis.py', 'R1/R2: decision rules, CV sweep, cost ratios'),
    ('robustness_experiments.py', 'Table 3: distributional stress tests'),
    ('crossover_order_analysis.py', 'R3: cycle order and period effects'),
    ('variance_components_analysis.py', 'R4: variance components, attrition'),
    ('efficiency_analysis.py', 'R5: joint design efficiency'),
    ('published_cohorts.py', 'R6: cohorts parameterized to published trials'),
    ('threshold_sensitivity.py', 'Table S4: responder-definition sensitivity'),
    ('nonresponder_fraction_sensitivity.py', 'Table S5: non-responder fraction'),
    ('threshold_approximation_check.py', 'Table S6: critical-value accuracy'),
]


def main():
    log = io.StringIO()
    failures = []
    print(f'Reproducing all manuscript numbers ({len(SCRIPTS)} scripts)\n')
    for name, blurb in SCRIPTS:
        path = HERE / name
        if not path.exists():
            print(f'  MISSING  {name}')
            failures.append(name)
            continue
        buf = io.StringIO()
        t0 = time.time()
        try:
            with warnings.catch_warnings():
                warnings.simplefilter('error')
                with contextlib.redirect_stdout(buf):
                    runpy.run_path(str(path), run_name='__main__')
            status = 'ok'
        except Exception as exc:                       # noqa: BLE001
            status = f'FAILED: {type(exc).__name__}: {exc}'
            failures.append(name)
        dt = time.time() - t0
        print(f'  {status:<10}{dt:>7.1f}s  {name:<40}{blurb}')
        log.write(f'\n{"#" * 78}\n# {name}\n{"#" * 78}\n')
        log.write(buf.getvalue())

    out = HERE / 'reproduce_all.log'
    out.write_text(log.getvalue())
    print(f'\nTranscript written to {out}')
    if failures:
        print(f'FAILURES: {", ".join(failures)}')
        return 1
    print('All scripts completed with no warnings.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
