devtools::load_all()
library(dplyr)
library(tidyr)
library(purrr)

# center/scale
cs <- function(.x) {    # .x = numeric vector
  out <- .x - mean(.x)  # center
  out / sd(out)         # scale
}

# prepare data set for analysis
cleanData <- example_data |>
  filter(SampleType == "Sample") |>          # rm control samples
  drop_na(Sex) |>                            # rm NAs if present
  log10() |>                                 # log10-transform (Math Generic)
  modify_at(getAnalytes(example_data), cs)   # center/scale analytes

t_tests <- getAnalyteInfo(cleanData) |>
  select(AptName, SeqId, Target = TargetFullName, EntrezGeneSymbol, UniProt)

t_tests <- t_tests |>
  mutate(
    formula = map(AptName, ~ as.formula(paste(.x, "~ Sex"))), # create formula
    t_test  = map(formula, ~ stats::t.test(.x, data = cleanData)),  # fit t-tests
    t_stat  = map_dbl(t_test, "statistic"),            # pull out t-statistic
    p.value = map_dbl(t_test, "p.value"),              # pull out p-values
    fdr     = p.adjust(p.value, method = "BH")         # FDR for multiple testing
  ) |>
  arrange(p.value) |>            # re-order by `p-value`
  mutate(rank = row_number())    # add numeric ranks

save(t_tests, file = "data/t_tests.rda", compress = "xz")
