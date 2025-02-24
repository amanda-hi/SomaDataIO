# Setup ----
ranks <- SomaDataIO::t_tests[,c("SeqId", "rank")]
rank_vec <- SomaDataIO::t_tests$rank
names(rank_vec) <- SomaDataIO::t_tests$SeqId


# Testing -----

test_that("`run_gsea()` errors when input data frame contains >2 columns", {
  rank_df <- SomaDataIO::t_tests[,c("SeqId", "t_test", "p.value", "rank")]
  expect_error(
    run_gsea(ranks = rank_df, pathIds = go2seqid$biological_process,
             adat = example_data),
            "`ranks` must be a 2-column data frame.")
})

test_that("`run_gsea()` errors when input vector is not correctly formatted", {
  rank_vec <- SomaDataIO::t_tests$SeqId
  names(rank_vec) <- SomaDataIO::t_tests$rank
  expect_error(
    run_gsea(ranks = rank_vec, pathIds = go2seqid$biological_process,
             adat = example_data),
    "Names of `ranks` must be SomaLogic SeqIds or AptNames.")
})
