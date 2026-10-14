# Consensus benchmark (development only)

The [pilot protocol](pilot-protocol.md) records the scientific design and its
prospective amendments. The current version is 0.2, dated 2026-10-14. This
directory is excluded from package builds and adds no public API or dependencies
to the installed `dar` package. Generated runs belong in ignored `results/`.

## Scope and safety

The CLI defaults to **preflight**: two strong-signal datasets, all five engines,
two identical-seed attempts each (20 fits). That mode never dispatches the
nine-cell, 180-dataset pilot or a cloud job. Documentation and a runnable command
are not authorization to launch either mode. The current amendment prepares
method changes only; its preflight and pilot have not been executed.

Every run uses a new output directory. The current checkout is installed into
`OUTPUT/runtime/library`; it does not replace the user's installed dar. Missing
native dependencies cause a stop; none are downloaded or installed automatically.
Workers run with clean R profiles, one thread, isolated native fits and a fixed
600-second timeout. The first implementation is sequential, within the protocol's
maximum two concurrent workers. Peak RSS is sampled externally, not confused
with R's vector-heap limit. No results are silently overwritten or resumed.

## Local preflight

From the repository root, in a terminal:

```sh
Rscript --vanilla dev/benchmark/run-benchmark.R \
  --mode=preflight \
  --output_dir=dev/benchmark/results/preflight-NEW-RUN-ID
```

Do not reuse an existing nonempty output directory. Exit status zero means the
technical preflight passed; nonzero means failures remain. Inspect `report.md`,
`checks.csv`, `runs.csv` and per-fit logs before deciding the next action.
Repeated preflight attempts are a reproducibility check, not independent
statistical replicates; their metrics are not passed to the scientific decision.

## Tests without native statistical fits

```sh
Rscript --vanilla dev/benchmark/tests/testthat.R
```

Tests use generated inputs and result-contract fixtures. They test the design,
truth, scoring, decision gates, actual bake/cool behavior and CLI lock. They do
not execute the 20 preflight fits or any full-pilot fits.

## Full pilot remains locked

Implementing a pilot runner is not authorization to execute it. It requires both
an explicit `confirm_pilot=true` argument and a passing preflight receipt supplied
via `preflight_dir`. The receipt must match hashes of package/benchmark source,
protocol, R/dependency versions and the empirical reservoir. A failed or stale
receipt cannot unlock it; there is no automatic expansion or best-rule search.

The runner keeps the four-engine voting denominator fixed. Native significance
is authoritative. A missing engine makes both consensuses unavailable, not empty
successful results. Native outputs, effective arguments and RNG states are
captured by observation in isolated processes. The benchmark does not patch
native functions or rethreshold package output; the package selects its effective
engine implementation, including the documented local LinDA correction.

The ANCOM vote uses native `diff_robust` with `pseudo_sens = TRUE`, including its
pseudocount sensitivity safeguard. The preflight must confirm that the package's
significance equals this native decision, not merely that the option was supplied.
LinDA retains adaptive zero handling and winsorization. With MicrobiomeStat 1.4,
dar selects a local copy of the native function with the single adaptive
`"Imputation"` assignment corrected to `"imputation"`; the installed dependency
is unchanged. Evidence must identify the effective backend and its verified
correction, not only the requested `adaptive = TRUE` argument.
MaAsLin3 requests centered abundance contrast coefficients. A native
incompatibility with a frozen method contract is a preflight failure, not a reason
to silently change settings or exclude that motor.

## Amendment history and existing runs

Protocol 0.1 was executed in `results/pilot-20261014-01`. Its ANCOM calls are
basic (`pseudo_sens = FALSE`) and its MicrobiomeStat 1.4 backend predates the local
LinDA correction. Those results, receipts and manifests remain historical evidence
of that contract; they are not relabelled as robust or corrected results.

Protocol 0.2 changes those two method contracts only. Seeds, scenarios, truth,
metrics, panel, consensus rules and decision gates are unchanged. Before any
later pilot execution, a newly authorized preflight must pass in a fresh output
directory with matching source/protocol/runtime identity. An old receipt or fit
cannot be reused under the changed contract. Neither a passing preflight nor this
amendment authorizes the pilot automatically.

## Files and evidence

- `fit-ledger.csv`: planned budget and attempt statuses, persisted before fits.
- `contract.rds/json`, `versions.csv`, `provenance.json`, `session-info.txt`:
  source/runtime identity and the run-local package installation.
- `reservoir/`: frozen bundled-data reference, donor-independence caveat and IDs.
- `datasets/`: count/metadata/truth/RNG bundles and actual consensus artifacts.
- `fits/`: job specifications, captured native output, prepared object and logs.
- `runs.csv`, `metrics.csv`, `checks.csv`, `report.md`: successes, missing results,
  failures, resource observations and explicit checks.
- `preflight-receipt.rds/json`: technical result only; not a publication claim.

The maintenance baseline in `dev/scientific-validation/` is neither read nor
promoted by this benchmark.
