# Setup ----
ranks <- SomaDataIO::t_tests[,c("SeqId", "rank")]
rank_vec <- SomaDataIO::t_tests$rank
names(rank_vec) <- SomaDataIO::t_tests$SeqId


# Testing -----
test_that("`run_gsea()` produces the expected output", {
  out <- run_gsea(ranks = ranks, pathIds = go2seqid$molecular_function,
                     adat = example_data)
  expect_s3_class(out, "data.frame")
  expect_length(colnames(out), 8L) # can change to dim when function is more finalized
})

test_that("`run_gsea()` errors when input data frame contains >2 columns", {
  big_df <- SomaDataIO::t_tests[,c("SeqId", "t_test", "p.value", "rank")]
  expect_error(
    run_gsea(ranks = big_df, pathIds = go2seqid$molecular_function,
             adat = example_data),
            "`ranks` must be a 2-column data frame.")
})

test_that("`run_gsea()` can handle SeqIds OR AptNames in the `ranks` object", {
  apt_vec <- SomaDataIO::t_tests$rank |>
    setNames(SomaDataIO::t_tests$AptName)
  out_vec <- run_gsea(ranks = apt_vec, pathIds = go2seqid$molecular_function,
                      adat = example_data)

  out_df <- run_gsea(ranks = SomaDataIO::t_tests[,c("AptName", "rank")],
                     pathIds = go2seqid$molecular_function,
                     adat = example_data)

  expect_true(all(is.SeqId(unlist(out_vec$leadingEdge))))
  # If a seed is not set, exact stats will differ slightly every time, but
  # these should be the same
  expect_identical(out_vec$pathway, out_df$pathway)
  expect_identical(out_vec$leadingEdge, out_df$leadingEdge)
})

test_that("`run_gsea()` errors when `ranks` vector is not correctly formatted", {
  rank_vec <- SomaDataIO::t_tests$SeqId
  names(rank_vec) <- SomaDataIO::t_tests$rank
  expect_error(
    run_gsea(ranks = rank_vec, pathIds = go2seqid$molecular_function,
             adat = example_data),
    "Names of `ranks` must be SomaLogic SeqIds or AptNames.")
})

test_that("`run_gsea()` errors when no analyes in `ranks` are found in `pathIds`", {
  ranks <- data.frame(SeqId = c("12345-67", "98765-43", "01234-56", "3579-01", "2468-02"),
                      Rank = seq(1:5L))
  expect_error(
    run_gsea(ranks = ranks, pathIds = go2seqid$molecular_function,
             adat = example_data)
  )
})

test_that("`selectAnalytes()` drops lower-ranked analytes targeting the same proteins", {
  res <- selectAnalytes(example_data, rank_vec)
  expect_length(res$filtered_results, 4744L)
  expect_type(res, "list")
  expect_type(res$removed_analytes, "character")
  expect_type(res$filtered_results, "integer")
})

test_that("`selectAnalytes()` is silent when `verbose=FALSE`", {
  expect_message(
    selectAnalytes(example_data, rank_vec),
                 "multi-analyte targets were dropped from the"
  )
  expect_silent(
    selectAnalytes(example_data, rank_vec, verbose = FALSE),
    "multi-analyte targets were dropped from the"
  )
})




