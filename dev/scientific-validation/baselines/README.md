# Accepted scientific-validation baseline

`accepted-full.csv` is generated only from reviewed full-profile results with
`promote-baseline.R`. It contains the metrics covered by active scientific
gates, their Monte Carlo intervals, comparison direction, tolerance, and source
workflow run. The initial baseline combines the completed engine artifacts from
runs `30824448541` and `30815107597`; the composite source is recorded in every
row and in `accepted-full.json`. Changing it requires review together with
`thresholds.csv` and must be followed by a second full run with a different
seed.
