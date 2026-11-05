# Toy study: 30 people, 5 samples each. Every person has 30 microbes of
# their own plus 50 that everyone shares, so each person has a clear
# fingerprint.
make_study <- function(seed = 1, n_subjects = 30, per = 5, own = 30, shared = 50) {
  set.seed(seed)
  n_features <- shared + n_subjects * own
  counts <- matrix(0L, n_features, n_subjects * per)
  rownames(counts) <- paste0("f", seq_len(n_features))
  colnames(counts) <- paste0("P", rep(seq_len(n_subjects), each = per), "_",
                             rep(seq_len(per), n_subjects))
  j <- 0
  for (s in seq_len(n_subjects)) for (k in seq_len(per)) {
    j <- j + 1
    counts[1:shared, j] <- ifelse(runif(shared) < 0.8, sample(1:50, shared, TRUE), 0)
    mine <- shared + (s - 1) * own + seq_len(own)
    counts[mine, j] <- ifelse(runif(own) < 0.7, sample(1:30, own, TRUE), 0)
  }
  subject <- paste0("P", rep(seq_len(n_subjects), each = per))
  list(counts = counts, subject = subject, per = per)
}
