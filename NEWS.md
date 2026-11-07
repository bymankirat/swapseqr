# swapseqr 0.2.0

* Uses the SwapSeq 1.1.0 engine. The power figures in `attr(result, "power")`
  are now measured: SwapSeq hides a few test swaps and relabels in copies of
  your labels, reruns the check and reports how many it caught (`swap`,
  `single`, `swapWithCheck`, `singleWithCheck`, plus `tests` and
  `falseFlagsPer1000`). The 1.0 estimate could read near 0% on clean tables
  that SwapSeq actually checks well.
* `power` is `NA` when too few people have three or more samples to hide tests in.

# swapseqr 0.1.0

* First release. `check_swaps()` runs the SwapSeq engine (the same JavaScript
  core as the web tool, via V8) on a phyloseq object or a count matrix.
* Tests for a planted swap, clean data, output format and input errors.
* Automatic R CMD check on Mac, Windows and Linux.
