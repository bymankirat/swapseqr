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
#' @return A data frame with one row per sample. Power estimates are in `attr(result, "power")`.
#' @importFrom methods as
#' @export
check_swaps <- function(x, subject, prior = 0.01,
                        well = NULL, plate = NULL, household = NULL) {

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

  # 3) Engine ke format mein data: har sample ke non-zero features (0-based index) + counts
  samples <- lapply(seq_len(ncol(counts)), function(j) {
    nz <- which(counts[, j] > 0)
    list(id = colnames(counts)[j], idx = I(nz - 1L), val = I(unname(counts[nz, j])))
  })
  opt <- function(v) if (is.null(v)) NULL else I(as.character(v))
  input <- list(features = I(rownames(counts)), samples = samples,
                subject = I(as.character(subject)), prior = prior,
                well = opt(well), plate = opt(plate), house = opt(household))
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
  message(sprintf("SwapSeq: %d strong flags, %d on the check list, out of %d samples. Power to catch swaps: %.0f%%.",
                  sum(res$level == "flag"), sum(res$level == "check"), nrow(res),
                  100 * swap_power))
  res
}
