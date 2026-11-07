test_that("a planted swap is flagged and the two samples point at each other", {
  st <- make_study()
  subject <- st$subject
  a <- 1; b <- st$per + 1                       # P1_1 and P2_1
  subject[c(a, b)] <- subject[c(b, a)]          # swap their labels

  res <- suppressMessages(check_swaps(st$counts, subject))
  flagged <- as.character(res$sample[res$level == "flag"])

  expect_setequal(flagged, c("P1_1", "P2_1"))
  expect_equal(res$looks_like[res$sample == "P1_1"], "P1")
  expect_equal(res$looks_like[res$sample == "P2_1"], "P2")
  expect_equal(res$partner[res$sample == "P1_1"], "P2_1")
  expect_equal(res$partner[res$sample == "P2_1"], "P1_1")
})

test_that("clean data gives no strong flags", {
  st <- make_study()
  res <- suppressMessages(check_swaps(st$counts, st$subject))
  expect_equal(sum(res$level == "flag"), 0)
})

test_that("output has one row per sample and the expected columns", {
  st <- make_study()
  res <- suppressMessages(check_swaps(st$counts, st$subject))
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), ncol(st$counts))
  expect_true(all(c("sample", "subject", "level", "probability",
                    "looks_like", "partner") %in% names(res)))
  expect_s3_class(res$level, "factor")
  expect_equal(levels(res$level), c("flag", "check", "ok"))
  expect_true(all(res$probability >= 0 & res$probability <= 1, na.rm = TRUE))
  expect_type(attr(res, "power"), "list")
})

test_that("check_swaps prints a one-line summary", {
  st <- make_study()
  expect_message(check_swaps(st$counts, st$subject), "SwapSeq:")
})

test_that("bad input gives a clear error", {
  st <- make_study()
  expect_error(check_swaps(st$counts, st$subject[-1]), "one subject ID per sample")
  bad <- st$counts
  bad[1, 1] <- -1
  expect_error(check_swaps(bad, st$subject), "non-negative")
})

test_that("phyloseq input gives the same answer as a matrix", {
  skip_if_not_installed("phyloseq")
  st <- make_study()
  subject <- st$subject
  subject[c(1, st$per + 1)] <- subject[c(st$per + 1, 1)]
  meta <- data.frame(person = subject, row.names = colnames(st$counts))
  ps <- phyloseq::phyloseq(phyloseq::otu_table(st$counts, taxa_are_rows = TRUE),
                           phyloseq::sample_data(meta))

  from_ps <- suppressMessages(check_swaps(ps, subject = "person"))
  from_mx <- suppressMessages(check_swaps(st$counts, subject))
  expect_equal(as.character(from_ps$level), as.character(from_mx$level))
  expect_equal(from_ps$probability, from_mx$probability)
})

test_that("power is measured on hidden test errors", {
  st <- make_study()
  res <- suppressMessages(check_swaps(st$counts, st$subject))
  pw <- attr(res, "power")
  expect_true(isTRUE(pw$measured))
  expect_gt(pw$tests$swapped, 0)
  # a toy study with strong individual signal: most hidden swaps are caught
  expect_gt(pw$swap, 0.7)
})
