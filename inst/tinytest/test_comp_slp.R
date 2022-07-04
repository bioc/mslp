load("example_z.rda")
load("example_comp_mut.rda")
res <- comp_slp(example_z, example_comp_mut)
expect_equal(head(res$slp_entrez), c("25879", "25879", "79075", "51001", "79075", "6045"))
