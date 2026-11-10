# swapseqr 0.3.2

* Small studies: a clear swap in a table of 24 to 40 samples is now a strong
  flag (60-70% caught on subsamples of real tables, up from 1-30%), with no
  more flags on clean labels. People need three or more samples for this.
* Missing subject labels (`NA` or empty) are left out of the check and returned
  with `scorable = FALSE`; before, they were treated as one person and could be
  flagged.
* Clear errors for duplicate sample IDs, a non-numeric table and a `prior`
  outside (0, 0.5).
* Relative abundances are turned into pseudo-counts instead of being swamped
  by the pseudocount.
* Two samples with no reads are no longer reported as duplicates; the
  duplicate check reports only exact copies in tables too small to judge.
* The "more swap pairs than the prior explains" message needs three or more
  pairs.

# swapseqr 0.3.1

* Several errors in one table: swap pairs where both samples point at each
  other's person are now strong flags when planted swaps stand out well in the
  table (`attr(result, "separation")` 0.7 or more). With 3% of samples swapped
  in 7 public shotgun tables, strong flags caught 41% of hidden swaps instead of
  5% (other samples flagged: 2.8 per 1,000 instead of 0.1). A message suggests
  `prior = 0.02` or `0.05` when a table has more swap pairs than the prior
  explains; at 3% that caught about 80%.
* Duplicate check: two samples labelled as different people whose communities
  are near-identical (one tube, extract or file under two names) are listed in
  `attr(result, "duplicates")`. With planted technical duplicates in 15 public
  shotgun tables: 97% found in cross-sectional tables, 63% in longitudinal ones,
  no unplanted pair reported.

# swapseqr 0.3.0

* `type = "markers"`: strain-level fingerprints for shotgun data. Give a
  presence table of strain markers, genes or SNVs (for example MetaPhlAn
  marker presence) and samples are compared on the rare markers they share.
* `host`: optional host-DNA read counts at common human SNPs (columns
  `sample`, `site`, `alt`, `total`). Each sample's reads are scored against
  every person's pooled reads, so samples can be matched to people by their own
  DNA where microbes change too fast to identify anyone (newborns, antibiotics).
  New columns `host_verdict`, `host_looks_like`, `host_llr`, `host_sites` and
  `flagged_by_host`. `host_verdict` is "same person", "other person",
  "identical DNA" (identical twins, or one person under two IDs), "unclear" or
  "unknown". Tested on simulated genotypes only so far, including a stress test
  with relatives, twins, uneven DNA and contamination: almost no false flags;
  samples need at least 1,500 covered SNP sites, and it works best with 6 or more
  samples per person.
* Counts mode gives the same results as 0.2.1 (SwapSeq 1.1.1 engine, November update).

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
