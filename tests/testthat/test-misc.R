test_that("rand_id works", {
  set.seed(123)
  expect_equal(rand_id(), "step_EoYnc")
})

test_that("get_comparisons works", {
  phy <- dar::metaHIV_phy
  var <- "RiskGroup2"

  is_list <- get_comparisons(var, phy)
  is_tibb <- get_comparisons(var, phy, as_list = FALSE)
  filtered <- get_comparisons(var, phy, n_cut = 30)

  expect_equal(class(is_list), "list")
  expect_s3_class(is_tibb, "tbl_df")
  expect_length(is_list, 3)
  expect_equal(nrow(is_tibb), 3)
  expect_length(filtered, 1)
})

test_that("to_tibble works", {
  phy <- dar::metaHIV_phy
  df <- phyloseq::otu_table(phy)

  test_1 <- to_tibble(df)
  test_2 <- to_tibble(df, id_name = "random_id")

  expect_length(test_1, 157)
  expect_s3_class(test_1, "tbl_df")
  expect_s3_class(test_2, "tbl_df")
  expect_equal(names(test_1)[1], "otu_id")
  expect_equal(names(test_2)[1], "random_id")
})

test_that("str_to_expr", {
  rec <-
    recipe(metaHIV_phy, "RiskGroup2", "Species") %>%
    step_subset_taxa(expr = 'Kingdom %in% c("Bacteria", "Archaea")') %>%
    step_filter_taxa(.f = "function(x) sum(x > 0) >= (0.02 * length(x))") %>%
    step_ancom()


  exprs <-
    rec@steps %>%
    purrr::map_chr(step_to_expr)

  expect_length(exprs, 3)
  expect_true(all(stringr::str_detect(exprs, "dar:::run_")))

})
