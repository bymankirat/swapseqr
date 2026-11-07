#' Check samples for swaps and wrong labels
#'
#' Runs SwapSeq on a microbiome count table and reports, for every sample,
#' whether it fits the subject on its label.
#'
#' @param x A phyloseq object, or a count matrix (features in rows, samples in columns).
#' @param subject For phyloseq: name of the sample_data column with subject IDs.
#'   For a matrix: a vector with one subject ID per sample.
#' @param prior Expected share of mislabelled samples (default 0.01 = 1%).
#' @param well,plate,household Optional: column names (phyloseq) or vectors (matrix).
#' @param type `"counts"` (default) for taxa, ASV or OTU counts. `"markers"` for
#'   strain-level tables from shotgun data (strain markers, genes or SNVs, for
#'   example MetaPhlAn marker presence): any value above 0 means present, and
#'   samples are compared on the rare markers they share.
#' @param host Optional host-DNA read counts from shotgun data, as a data frame
#'   with columns `sample`, `site`, `alt` and `total` (reads showing the
#'   alternative allele and all reads at a common human SNP). Samples are then
#'   also matched to people by their own DNA, which works where microbes change
#'   too fast to identify anyone (newborns, antibiotics). Samples need at least
#'   1,500 covered SNP sites, and it works best with 6 or more samples per
#'   person; tested on simulated data only. Human genetic data: keep it on your own
#'   computer and follow your ethics approval.
#' @return A data frame with one row per sample. Measured power (share of test
#'   swaps and relabels hidden in copies of your labels that were caught) is in
#'   `attr(result, "power")`. With `host`, the columns `host_verdict`,
#'   `host_looks_like` and `flagged_by_host` show what the host DNA says.
#' @importFrom methods as
#' @export
check_swaps <- function(x, subject, prior = 0.01,
                        well = NULL, plate = NULL, household = NULL,
                        type = c("counts", "markers"), host = NULL) {
  type <- match.arg(type)

  # 1) Count matrix (features x samples) aur labels nikalo
  if (inherits(x, "phyloseq")) {
    if (!requireNamespace("phyloseq", quietly = TRUE))
      stop("Please install the phyloseq package.")
    counts <- as(phyloseq::otu_table(x), "matrix")
    if (!phyloseq::taxa_are_rows(x)) counts <- t(counts)      # features rows mein chahiye
    meta <- as(phyloseq::sample_data(x), "data.frame")
    meta <- meta[colnames(counts), , drop = FALSE]             # samples ka order match karo
    pick <- function(col) {
      if (is.null(col)) return(NULL)
      if (!col %in% names(meta)) stop("Column '", col, "' is not in sample_data.")
      as.character(meta[[col]])
    }
    subject <- pick(subject); well <- pick(well)
    plate <- pick(plate); household <- pick(household)
  } else {
    counts <- as.matrix(x)
  }

  # 2) Basic checks
  if (is.null(colnames(counts))) colnames(counts) <- paste0("sample", seq_len(ncol(counts)))
  if (is.null(rownames(counts))) rownames(counts) <- paste0("feature", seq_len(nrow(counts)))
  if (length(subject) != ncol(counts)) stop("Need exactly one subject ID per sample.")
  if (any(counts < 0, na.rm = TRUE)) stop("Counts must be non-negative.")
  counts[is.na(counts)] <- 0

  # 3) Engine ke format mein data: har sample ke non-zero features (0-based index) + counts
  samples <- lapply(seq_len(ncol(counts)), function(j) {
    nz <- which(counts[, j] > 0)
    list(id = colnames(counts)[j], idx = I(nz - 1L), val = I(unname(counts[nz, j])))
  })
  opt <- function(v) if (is.null(v)) NULL else I(as.character(v))
  input <- list(features = I(rownames(counts)), samples = samples,
                subject = I(as.character(subject)), prior = prior,
                well = opt(well), plate = opt(plate), house = opt(household),
                type = type, host = host_input(host, colnames(counts)))
  json <- jsonlite::toJSON(input, auto_unbox = TRUE, null = "null", digits = NA)

  # 4) JavaScript engine chalao (wahi jo website par chalta hai)
  ctx <- V8::v8()
  ctx$source(system.file("js", "swapseq-core.js", package = "swapseqr"))
  out <- jsonlite::fromJSON(ctx$call("swapseq_run", V8::JS(json)))

  # 5) Result ko saaf data frame banao
  res <- out$rows
  res$level <- factor(res$level, levels = c("flag", "check", "ok"))
  attr(res, "power") <- out$power
  attr(res, "shifts") <- out$shifts
  swap_power <- if (is.null(out$power$swap)) NA_real_ else out$power$swap
  attr(res, "host") <- out$host
  message(sprintf("SwapSeq: %d strong flags, %d on the check list, out of %d samples. Measured power to catch swaps: %.0f%%.",
                  sum(res$level == "flag"), sum(res$level == "check"), nrow(res),
                  100 * swap_power))
  if (!is.null(out$host)) {
    if (isTRUE(out$host$ok)) {
      message(sprintf("Host DNA: %d samples with enough host reads; %d flagged by host DNA.",
                      out$host$samples_with_host_dna, sum(res$flagged_by_host)))
    } else message("Host DNA: not used (", out$host$reason, ").")
  }
  res
}

# host read counts (long data frame) -> one list per sample: site -> c(alt, total)
host_input <- function(host, sample_ids) {
  if (is.null(host)) return(NULL)
  host <- as.data.frame(host)
  need <- c("sample", "site", "alt", "total")
  if (!all(need %in% names(host)))
    stop("`host` needs the columns sample, site, alt and total.")
  host <- host[host$total > 0 & host$sample %in% sample_ids, need]
  if (any(host$alt > host$total | host$alt < 0))
    stop("In `host`, alt must be between 0 and total.")
  by <- split(host, factor(host$sample, levels = sample_ids))
  unname(lapply(by, function(d) {
    if (!nrow(d)) return(NULL)
    stats::setNames(lapply(seq_len(nrow(d)), function(k) c(d$alt[k], d$total[k])),
                    as.character(d$site))
  }))
}
