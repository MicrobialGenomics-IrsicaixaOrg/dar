# Scientific validation

This directory contains development-only Monte Carlo validation for the centralized statistical models in `dar`. It is excluded from package builds and does not add user-facing functions.

## Quick local run

From the repository root, run:

```bash
Rscript --vanilla dev/scientific-validation/run-validation.R \
  --profile=quick \
  --engine=all \
  --scenario=all \
  --output_dir=validation-results/quick/all
```

Optional engines whose dependencies are unavailable are recorded as `dependency_missing`. A compatible engine that starts but cannot fit the model is recorded as `engine_error` and makes a strict run fail.

Use comma-separated engine or scenario names to restrict a run. `--engine=none` performs simulation and artifact smoke checks without fitting a DA method. The full statistical profile is selected with `--profile=full`.

## Scenarios

The suite covers a global null, balanced positive and negative effects, independent zero inflation, fixed confounding, condition-by-time effects, repeated subjects with random effects, and an explicit complete-case cohort for missing and unbalanced metadata.

Counts are generated from centered latent log abundances, converted to
compositions with a softmax, and sampled with variable multinomial library
sizes. The scheduled workflow varies the base seed between runs, so regression
detection does not depend on identical Monte Carlo draws. Each simulation
produces an exact taxon-contrast truth table and a manifest.

## Metrics and estimands

Power, type-I error, FDR, direction recovery, finite-result rate, execution rate, bias, and interval coverage are kept separate by engine, scenario, replicate, contrast, and effect metric. Bias and coverage are marked unavailable when an engine effect is not on a scientifically comparable scale. The suite never requires different engines to return equal p-values or effect magnitudes.

The quick profile uses two replicates and treats Monte Carlo thresholds as
informational. The full profile uses the replicate counts in `profiles.csv`.
Gates in `thresholds.csv` match engines, scenarios and contrasts, so power and
direction requirements apply only to deliberately strong effects. LinDA bias
and coverage are gated because its log2 fold-change scale is comparable with
the simulation truth. The direction-reversal sentinel must always fail its
gate; this protects numerator-minus-denominator orientation.

`baselines/accepted-full.csv` stores reviewed full-run intervals and
metric-specific tolerances. Use `promote-baseline.R REPORT_DIR RUN_ID` only
after reviewing a complete run, commit the CSV and metadata together, and then
confirm them with a different seed. Threshold and tolerance changes require
the scientific justification in the same pull request.

## Artifacts

Every run writes CSV tables, a recursive dependency-version manifest, one RDS
bundle and session information. The workflow splits ANCOM-BC2 by scenario,
detects missing jobs, renders `scientific-validation.html` and writes key status
counts to the Actions summary. Artifacts are retained for 90 days. Complete
`full/all` runs from `devel` additionally archive compact reports and
machine-readable tables on the `validation-results` branch. Quick and partial
runs never update that durable history.
