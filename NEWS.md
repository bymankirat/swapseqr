# swapseqr 0.3.0

* `type = "markers"`: strain-level fingerprints for shotgun data. Give a
  presence table of strain markers, genes or SNVs (for example MetaPhlAn
  marker presence) and samples are compared on the rare markers they share.
* `host`: optional host-DNA read counts at common human SNPs (columns
  `sample`, `site`, `alt`, `total`). Each sample's reads are scored against
  every person's pooled reads, so samples can be matched to people by their own
  DNA where microbes change too fast to identify anyone (newborns, antibiotics).
  New columns `host_verdict`, `host_looks_like`, `host_llr`, `host_sites` and
  `flagged_by_host`. Tested on simulated genotypes only so far: samples need at least
  2,000 covered SNP sites, where it caught 95% of swaps (all of them from 2,500).
* Counts mode gives exactly the same results as 0.2.1 (SwapSeq engine 1.2.0).

# swapseqr 0.2.1

* SwapSeq 1.1.1 engine: measured power no longer reads low when a table
  already holds errors (near-certain flags are set aside before the hidden
  tests). Flags and probabilities are unchanged.

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
