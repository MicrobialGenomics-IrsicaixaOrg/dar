#' lefse analysis
#'
#' Lefser is metagenomic biomarker discovery tool that is based on LEfSe tool
#' and is published by Huttenhower et al. 2011. Lefser is the R implementation
#' of the LEfSe method. Using statistical analyses, lefser compares microbial
#' populations of healthy and diseased subjects to discover differencially
#' expressed microorganisms. Lefser than computes effect size, which estimates
#' magnitude of differential expression between the populations for each
#' differentially expressed microorganism. Subclasses of classes can also be
#' assigned and used within the analysis.
#'
#' @param rec A Recipe object. The step will be added to the sequence of
#'   operations for this Recipe.
#' @param kruskal.threshold numeric(1) The p-value for the Kruskal-Wallis Rank
#'   Sum Test (default 0.05).
#' @param wilcox.threshold numeric(1) The p-value for the Wilcoxon Rank-Sum Test
#'   when 'blockCol' is present (default 0.05).
#' @param lda.threshold numeric(1) The effect size threshold (default 2.0).
#' @param blockCol character(1) Optional column name in 'colData(expr)'
#'   indicating the blocks, usually a factor with two levels (e.g., 'c("adult",
#'   "senior")'; default NULL).
#' @param assay The i-th assay matrix in the ‘SummarizedExperiment' (’expr';
#'   default 1).
#' @param trim.names If 'TRUE' extracts the most specific taxonomic rank of
#'   organism.
#' @param rarefy Boolean indicating if OTU counts must be rarefyed. This
#'   rarefaction uses the standard R sample function to resample from the
#'   abundance values in the otu_table component of the first argument, physeq.
#'   Often one of the major goals of this procedure is to achieve parity in
#'   total number of counts between samples, as an alternative to other formal
#'   normalization procedures, which is why a single value for the sample.size
#'   is expected. If 'no_seed', rarefaction is performed without a set seed. 
#' @param id A character string that is unique to this step to identify it.
#'
#' @include recipe-class.R
#' @family Diff taxa steps
#' @aliases step_lefse
#' @return An object of class `Recipe`
#' @export
#' @autoglobal
#' @tests
#' data(metaHIV_phy)
#' 
#' test <-
#'  recipe(metaHIV_phy, "RiskGroup2", "Phylum") |>
#'  step_subset_taxa(tax_level = "Kingdom", taxa = c("Bacteria", "Archaea")) |>
#'  step_filter_by_prevalence() |> 
#'  step_lefse() |> 
#'  step_lefse(rarefy = FALSE) 
#'  
#' expect_match(
#'   capture_warnings(prep(test)),
#'   "Convert counts to relative abundances with 'relativeAb\\(\\)'",
#'   all = TRUE
#' ) |>
#'   expect_snapshot()
#' 
#' data(test_prep_rec)
#' expect_error(step_lefse(test_prep_rec))
#' @examples 
#' data(metaHIV_phy)
#' 
#' ## Init Recipe
#' rec <- 
#'   recipe(metaHIV_phy, "RiskGroup2", "Phylum") |>
#'   step_subset_taxa(tax_level = "Kingdom", taxa = c("Bacteria", "Archaea")) |>
#'   step_filter_taxa(.f = "function(x) sum(x > 0) >= (0.3 * length(x))")
#' 
#' rec
#' 
#' ## Define step with default parameters and prep
#' rec <- 
#'   step_lefse(rec) |>
#'   prep(parallel = FALSE)
#'   
#' rec
#' 
#' ## Running lefse without rarefaction (not recommended)
#' rec <- 
#'   recipe(metaHIV_phy, "RiskGroup2", "Species") |>
#'   step_lefse(rarefy = FALSE)
#'   
#' rec
methods::setGeneric(
  name = "step_lefse",
  def = function(rec,
                 kruskal.threshold = 0.05,
                 wilcox.threshold = 0.05,
                 lda.threshold = 2,
                 blockCol = NULL,
                 assay = 1L,
                 trim.names = FALSE,
                 rarefy = TRUE,
                 id = rand_id("lefse")) {
    standardGeneric("step_lefse")
  }
)

#' @rdname step_lefse
#' @export
#' @autoglobal
methods::setMethod(
  f = "step_lefse",
  signature = c(rec = "Recipe"),
  definition = function(rec,
                        kruskal.threshold,
                        wilcox.threshold,
                        lda.threshold,
                        blockCol,
                        assay,
                        trim.names,
                        rarefy,
                        id) {

    recipes_pkg_check(required_pkgs_lefse(), "step_lefser()")
    if (!rarefy & !contains_rarefaction(rec)) {
      rlang::inform(c(
        "!" = glue::glue(
          "Run lefse without rarefaction is not recommended", 
          " ({crayon::blue(paste0('id = ', id))})"
        )
      ))
    }

    add_step(
      rec,
      step_lefse_new(
        kruskal.threshold = kruskal.threshold,
        wilcox.threshold = wilcox.threshold,
        lda.threshold = lda.threshold,
        blockCol = blockCol,
        assay = assay,
        trim.names = trim.names,
        rarefy = rarefy,
        id = id
      )
    )
  }
)

#' @rdname step_lefse
#' @export
#' @autoglobal
methods::setMethod(
  f = "step_lefse",
  signature = c(rec = "PrepRecipe"),
  definition = function(rec,
                        kruskal.threshold,
                        wilcox.threshold,
                        lda.threshold,
                        blockCol,
                        assay,
                        trim.names,
                        rarefy,
                        id) {
    rlang::abort("This function needs a non-PrepRecipe!")
  }
)

#' @noRd
#' @keywords internal
#' @autoglobal
step_lefse_new <-
  function(kruskal.threshold,
           wilcox.threshold,
           lda.threshold,
           blockCol,
           assay,
           trim.names,
           rarefy,
           id) {
    step(
      subclass = "lefse",
      kruskal.threshold = kruskal.threshold,
      wilcox.threshold = wilcox.threshold,
      lda.threshold = lda.threshold,
      blockCol = blockCol,
      assay = assay,
      trim.names = trim.names,
      rarefy = rarefy,
      id = id
    )
  }

#' @noRd
#' @keywords internal
#' @autoglobal
required_pkgs_lefse <- function(x, ...) { 
  c("bioc::lefser", "bioc::SummarizedExperiment") 
}

#' @noRd
#' @keywords internal
#' @autoglobal
run_lefse <-
  function(rec,
           kruskal.threshold = kruskal.threshold,
           wilcox.threshold = wilcox.threshold,
           lda.threshold = lda.threshold,
           blockCol = blockCol,
           assay = assay,
           trim.names = trim.names,
           rarefy = rarefy ) {

  lefse_mat <- prepro_lefse(rec, rarefy)
  vars <- get_var(rec)
  vars %>%
    purrr::set_names() %>%
    purrr::map(function(var) {
      get_comparisons(var, get_phy(rec), as_list = TRUE, n_cut = 1) %>%
        purrr::map_dfr(function(comparison) {
          sample_data <- 
            dplyr::filter(sample_data(rec), !!dplyr::sym(var) %in% comparison)
         
          se <- SummarizedExperiment::SummarizedExperiment(
            assays = list(counts = lefse_mat[, sample_data$sample_id]),
            colData = sample_data
          )

          # set.seed(1234)
          lefse_res <- lefser::lefser(
            se,
            groupCol = var,
            kruskal.threshold = 1,
            wilcox.threshold = 1,
            lda.threshold = 0,
            blockCol = blockCol,
            assay = assay,
            trim.names = trim.names
          )
          
          adjpval <-
            ifelse(is.null(blockCol), kruskal.threshold, wilcox.threshold)
          
          lefse_res %>%
            tibble::as_tibble() %>%
            dplyr::rename(otu = Names) %>%
            dplyr::mutate(otu = stringr::str_remove_all(otu, "\`")) %>%
            dplyr::left_join(
              kruskal_test(se, sample_data[[var]]), by = "otu"
            ) %>%
            dplyr::arrange(pvalue, dplyr::desc(abs(scores))) %>%
            dplyr::mutate(
              comparison = stringr::str_c(comparison, collapse = "_"),
              var = var,
              taxa = stringr::str_remove_all(otu, ".*[|]"), 
              effect = scores,
              signif = ifelse(
                adjp < adjpval & abs(scores) >= lda.threshold,
                TRUE,
                FALSE
              )
            ) %>%
            dplyr::left_join(tax_table(rec), by = "taxa") %>%
            dplyr::rename(lefse_id = otu)
        })
    })
}

#' @noRd
#' @keywords internal
#' @autoglobal
prepro_lefse <- function(rec, rarefy) {

  tax_level <- get_tax(rec)
  phy <- 
    get_phy(rec) %>% 
    use_rarefy(rarefy)

  # Defining iterating tax levels -------------------------------------------
  tax_otp <- 
    c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species") %>% 
    .[. %in% colnames(phy@tax_table)]
  
  tax_lev_f <- tax_otp[seq_len(match(tax_level, tax_otp))]

  # Computing output table --------------------------------------------------
  lefse_input <-
    purrr::map2_dfr(tax_lev_f, seq_along(tax_lev_f), function(tax, it) {
      dat <- phyloseq::tax_glom(phy, taxrank = tax, NArm = FALSE)
      abundance_df <- phyloseq::otu_table(dat) %>% to_tibble("RTC")
      tax_lev_names <-
        dat %>%
        phyloseq::tax_table(object = .) %>%
        data.frame() %>%
        tibble::as_tibble() %>%
        dplyr::select(dplyr::all_of(seq_len(it))) %>%
        dplyr::mutate(dplyr::across(
          dplyr::everything(),
          .fns = function(x) stringr::str_replace_all(x, " |[.]|-", "_")
        )) 
        # dplyr::mutate(dplyr::across(
        #   .fns = function(x)
        #     stringr::str_replace_all(x, "\\[|\\]", "")
        # )) %>% 

      if (it == 1) {
        tax_lev_names <- tax_lev_names %>% dplyr::pull(!!tax)
      } else {
        tax_lev_names <-
          tax_lev_names %>%
          tidyr::unite('name', Kingdom:(!!tax), sep = "|") %>%
          dplyr::pull()
      }
      abundance_df <- abundance_df %>% dplyr::mutate(RTC = tax_lev_names)
    }) %>%
    dplyr::filter(
      stringr::str_count(RTC, "[|]") == match(tax_level, tax_otp) - 1
    ) %>%
    dplyr::arrange(RTC) %>%
    data.frame(row.names = 1, check.names = FALSE) %>%
    as.matrix()

  lefse_input
}

#' @noRd
#' @keywords internal
#' @autoglobal
kruskal_test <- function(se, levels, assay = 1L) {
  group <- levels %>% as.factor() %>% as.numeric()
  expr <- SummarizedExperiment::assay(se, i = assay)
  apply(expr, 1L, function(x) {
    stats::kruskal.test(x ~ group)[["p.value"]]
  }) %>%
    tibble::as_tibble(rownames = "otu") %>%
    dplyr::rename(pvalue = value) %>%
    dplyr::mutate(adjp = stats::p.adjust(pvalue, method = "BH"))
}
