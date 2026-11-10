# swapseqr

<!-- badges: start -->
[![R-CMD-check](https://github.com/bymankirat/swapseqr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bymankirat/swapseqr/actions/workflows/R-CMD-check.yaml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23963612.svg)](https://doi.org/10.5281/zenodo.23963612)
<!-- badges: end -->

R interface to [SwapSeq](https://bymankirat.github.io/swapseq/): it finds swapped tubes and mislabelled samples in longitudinal microbiome studies, before they reach your results.

When the same people are sampled more than once, each person's microbiome works as a fingerprint. `check_swaps()` checks whether every sample matches the person on its label. It reports:

- **swaps**: two samples that each look like the other's person
- **mislabels**: a sample that looks like a different person
- **sample-sheet shifts**: a run of samples on a plate where every label is off by the same number of wells

The method is the same JavaScript core used by the SwapSeq web tool. It runs inside R through [V8](https://cran.r-project.org/package=V8), so the R package and the web tool give identical results. Your data never leave your computer.

## Install

```r
# install.packages("remotes")
remotes::install_github("bymankirat/swapseqr")
```

## Use

With a phyloseq object, give the name of the sample-data column that holds the subject (person) ID:

```r
library(swapseqr)
res <- check_swaps(ps, subject = "subject_id")
res[res$level != "ok", ]
```

With a count matrix (features in rows, samples in columns), give the subject IDs as a vector:

```r
res <- check_swaps(counts, subject = meta$subject_id)
```

Optional extras:

- `well`, `plate`: plate positions, used to detect sample-sheet shifts
- `household`: household IDs, so that people who live together are not over-flagged

### Shotgun data (0.3.0)

Strain-level fingerprints from a marker or gene presence table (for example MetaPhlAn marker presence),
features in rows and samples in columns:

```r
res <- check_swaps(markers, subject = meta$subject_id, type = "markers")
```

Host DNA: shotgun reads carry some of the person's own DNA, which identifies them even when their microbes
change too fast (newborns, antibiotics). Give read counts at common human SNPs as a data frame with columns
`sample`, `site`, `alt` and `total`:

```r
res <- check_swaps(ps, subject = "subject_id", host = host_counts)
res[res$flagged_by_host, ]
```

The host-DNA check is tested in simulation only so far (including relatives, twins and contamination), needs at least 1,500 covered SNP sites per sample and works best with 6 or more samples per person.
How to make `host_counts`, and the privacy rules for human genetic data:
<https://github.com/bymankirat/swapseq/blob/main/docs/host-dna.md>.

## Output

You get one row per sample:

| column | meaning |
|---|---|
| `level` | `flag` (strong evidence), `check` (worth a look), `ok` |
| `kind` | `swap` or `mismatch` |
| `probability` | estimated chance the label is wrong |
| `looks_like` | the subject this sample most resembles |
| `partner` | the other sample in a likely swap |
| `rare_*` | how many rare features it shares with its own subject and with the look-alike |

`attr(res, "power")` gives the measured power: SwapSeq hides a few test swaps and relabels in copies of your labels, reruns the check, and reports the share it caught.

## Limits

- Each subject needs at least two samples. A sample whose subject has no other samples cannot be checked.
- Detection power depends on the body site and the study. In the [SwapSeq Atlas](https://bymankirat.github.io/swapseq/atlas/) of 42 public tables it was high for adult stool, forehead skin and vaginal swabs, and very low for infant stool, patients on heavy antibiotics, lab mice from one donor and genus-level tables.
- A flag is evidence, not proof. Check flagged samples against lab records before you relabel or drop anything.
- This is version 0.3.0 (SwapSeq 1.1.1 engine, method v3). Please report problems at <https://github.com/bymankirat/swapseqr/issues>.

## Cite

Singh M. (2026). *SwapSeq: a sample mix-up checker for longitudinal microbiome studies* (v1.1.1). Zenodo. <https://doi.org/10.5281/zenodo.23963612> (RRID:SCR_029248)

## Author

Mankirat Singh ([ORCID 0009-0005-3375-2123](https://orcid.org/0009-0005-3375-2123))

## License

GPL (>= 3)
