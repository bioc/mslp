load("example_expr.rda")
load("example_corr_mut.rda")
set.seed(42)
res <- corr_slp(example_expr, example_corr_mut)
expect_equal(head(res$slp_entrez), c("751071", "10395", "57864", "54629", "54898", "1778"))
