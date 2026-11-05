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

`attr(res, "power")` estimates how likely this dataset is to catch a swap, given how distinct its subjects are.

## Limits

- Each subject needs at least two samples. A sample whose subject has no other samples cannot be checked.
- Detection power depends on the body site and the study. It is high for gut, lower for sites that vary a lot within a person.
- A flag is evidence, not proof. Check flagged samples against lab records before you relabel or drop anything.
- This is version 0.1.0, an early release. Please report problems at <https://github.com/bymankirat/swapseqr/issues>.

## Cite

Singh M. SwapSeq: detecting sample swaps in longitudinal microbiome studies. Zenodo. <https://doi.org/10.5281/zenodo.23963612> (RRID:SCR_029248)

## Author

Mankirat Singh ([ORCID 0009-0005-3375-2123](https://orcid.org/0009-0005-3375-2123))

## License

GPL (>= 3)
