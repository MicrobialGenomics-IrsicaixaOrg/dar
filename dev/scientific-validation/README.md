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

Counts are generated from centered latent log abundances, converted to compositions with a softmax, and sampled with variable multinomial library sizes. Seeds are derived only from the base seed, scenario, and replicate. Each simulation produces the microbiome data, an exact taxon-contrast truth table, and a manifest.

## Metrics and estimands

Power, type-I error, FDR, direction recovery, finite-result rate, execution rate, bias, and interval coverage are kept separate by engine, scenario, replicate, contrast, and effect metric. Bias and coverage are marked unavailable when an engine effect is not on a scientifically comparable scale. The suite never requires different engines to return equal p-values or effect magnitudes.

The quick profile uses two replicates and treats Monte Carlo thresholds as informational. The full profile uses the replicate counts in `profiles.csv`. Gates in `thresholds.csv` compare confidence intervals across independent replicates with deliberately tolerant scientific limits. The direction-reversal sentinel must always fail its gate; this protects numerator-minus-denominator orientation.

## Artifacts

Every run writes CSV tables, one RDS bundle, and session information. The dedicated GitHub Actions workflow runs each engine independently, uploads per-engine artifacts, aggregates them, and renders `scientific-validation.html`. It runs only when manually dispatched or on the monthly schedule, so pull-request checks are unaffected.
