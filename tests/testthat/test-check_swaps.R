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
  expect_gt(pw$swapWithCheck, 0.7)
})

test_that("markers mode finds a planted swap in a presence table", {
  st <- make_study(own = 100)
  pres <- (st$counts > 0) * 1L
  subject <- st$subject
  subject[c(1, st$per + 1)] <- subject[c(st$per + 1, 1)]
  res <- suppressMessages(check_swaps(pres, subject, type = "markers"))
  expect_setequal(as.character(res$sample[res$level == "flag"]), c("P1_1", "P2_1"))
})

test_that("host DNA catches a swap that the microbes cannot see", {
  st <- make_study(n_subjects = 12, per = 5)
  # no personal microbes at all: everyone looks alike, as in newborns
  flat <- st$counts
  flat[51:nrow(flat), ] <- 0
  host <- make_host(colnames(flat), st$subject, sites = 4000)
  subject <- st$subject
  subject[c(1, st$per + 1)] <- subject[c(st$per + 1, 1)]
  res <- suppressMessages(check_swaps(flat, subject, host = host))
  expect_true(isTRUE(attr(res, "host")$ok))
  expect_setequal(as.character(res$sample[res$flagged_by_host]), c("P1_1", "P2_1"))
  expect_equal(res$host_looks_like[res$sample == "P1_1"], "P1")
  expect_equal(sum(res$level == "flag" & !res$sample %in% c("P1_1", "P2_1")), 0)
})

test_that("host input is checked", {
  st <- make_study()
  expect_error(check_swaps(st$counts, st$subject, host = data.frame(sample = "x")),
               "sample, site, alt and total")
})

test_that("a near-identical sample under another person is reported as a possible duplicate", {
  st <- make_study()
  x <- st$counts
  src <- which(colnames(x) == "P6_2")
  dst <- which(colnames(x) == "P5_3")
  set.seed(9)
  x[, dst] <- ifelse(x[, src] > 0, pmax(1L, round(x[, src] * runif(nrow(x), 0.85, 1.15))), 0L)
  res <- suppressMessages(check_swaps(x, st$subject))
  d <- attr(res, "duplicates")
  expect_equal(nrow(d), 1)
  expect_setequal(c(d$sample_a, d$sample_b), c("P5_3", "P6_2"))
  clean <- suppressMessages(check_swaps(st$counts, st$subject))
  expect_equal(nrow(attr(clean, "duplicates")), 0)
})

test_that("missing subject labels and bad input are handled", {
  st <- make_study()
  lab <- st$subject; lab[c(4, 17)] <- NA
  expect_message(res <- check_swaps(st$counts, lab), "no subject label")
  expect_equal(nrow(res), ncol(st$counts))
  expect_false(any(res$scorable[c(4, 17)]))
  expect_equal(sum(res$level == "flag"), 0)
  x <- st$counts; colnames(x)[2] <- colnames(x)[1]
  expect_error(check_swaps(x, st$subject), "more than once")
  expect_error(check_swaps(st$counts, st$subject, prior = 5), "prior must be")
})
