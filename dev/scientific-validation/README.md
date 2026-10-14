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

Power, type-I error, FDR, direction recovery, finite-result rate, execution rate,
signed mean error (`bias`), relative mean absolute error (`relative_mae`), and
interval coverage are kept separate by engine, scenario, replicate, contrast,
and effect metric. Estimation errors and coverage are unavailable when an engine
effect is not on a scientifically comparable scale. The suite never requires
different engines to return equal p-values or effect magnitudes.

The scoring universe is the truth table crossed with attempted compatible runs,
not the rows returned by an engine. Failed fits remain in that universe;
capability skips and unavailable dependencies do not. Missing taxon-contrast
rows are not discoveries, and missing or non-finite effects are not recovered
directions. Power and direction recovery therefore use all true signals, while
`finite_rate` uses all expected hypotheses. Native output filtering can reduce
these end-to-end rates; a missing row is never silently assumed to be a valid
filtering decision. Missing effects, p-values and standard errors remain missing.

Per-replicate metrics include `n_expected`, `n_returned`, `n_missing`, `n_signals`,
`n_tested`, `n_finite`, `n_signal_finite` and `n_intervals`, making their
denominators explicit. `bias` and `relative_mae` are conditional on finite
estimates of true signals; interval coverage additionally requires a finite,
non-negative standard error. These conditional diagnostics must be interpreted
together with availability, not as evidence that missing estimates were correct.

`relative_mae` is `mean(abs(estimate - scaled_truth)) /
mean(abs(scaled_truth))` over the same finite signal estimates. Taking absolute
errors before averaging prevents balanced positive and negative errors from
cancelling. It measures estimation error, not signed bias, and is unbounded
above. Null contrasts without true signals have no relative-error metric.

The quick profile uses two replicates and treats Monte Carlo thresholds as
informational. The full profile uses the replicate counts in `profiles.csv`.
Gates in `thresholds.csv` match engines, scenarios and contrasts, so power and
direction requirements apply only to deliberately strong effects. The new
relative-error gate uses an upper limit of 0.50 for strong condition contrasts
on comparable DESeq2, ANCOM-BC2 and LinDA scales. This is a conservative
gross-error criterion: mean absolute error above half the known effect magnitude
is unacceptable, and complete attenuation to zero gives 1. As with other upper
gates, failure requires the Monte Carlo interval to lie entirely above the
limit. It is not inherited
from the former 0.30 signed-bias gate, an empirical calibration, or a claim of
general statistical accuracy. Weak time and difference-in-differences effects
remain descriptive. DESeq2 and LinDA interval coverage is gated because their
log2 fold-change scales are comparable with the simulation truth. The
direction-reversal sentinel must always fail its gate; this protects
numerator-minus-denominator orientation.

`baselines/accepted-full.csv` stores reviewed full-run intervals and
metric-specific tolerances. Use `promote-baseline.R REPORT_DIR RUN_ID` only
after reviewing a complete run, commit the CSV and metadata together, and then
confirm them with a different seed. Threshold and tolerance changes require
the scientific justification in the same pull request.

The historical baseline predates complete-universe scoring and contains the
retired `abs_relative_bias` metric. Its values must not be renamed to
`relative_mae`. Review new full-profile results and calibrate a new baseline
tolerance before promotion; promotion rejects metrics without a reviewed
tolerance. Existing historical baseline files are preserved, not silently
rewritten or treated as a calibrated reference for the new metric.

## Artifacts

Every run writes CSV tables, a recursive dependency-version manifest, one RDS
bundle and session information. The workflow splits ANCOM-BC2 by scenario,
detects missing jobs, renders `scientific-validation.html` and writes key status
counts to the Actions summary. Artifacts are retained for 90 days. Complete
`full/all` runs from `devel` additionally archive compact reports and
machine-readable tables on the `validation-results` branch. Quick and partial
runs never update that durable history.
