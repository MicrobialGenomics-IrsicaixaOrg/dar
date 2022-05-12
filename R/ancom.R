#' @export
step_ancom <- function(recipe, ..., skip = FALSE, id = rand_id("ancom")) {
  inputs <- enquos(...)

  add_step(
    recipe,
    step_ancom_new(
      inputs = inputs,
      skip = skip,
      id = id
    )
  )
}

#' @export
step_ancom_new <- function(inputs, skip, id) {
  step(
    subclass = "ancom",
    input = inputs,
    skip = skip,
    id = id
  )
}

#' Multiple comparison ANCOM
#'
#' @param var Categorical variable
#' @param phy
#' @param tax_level
#'
#' @return list
#' @noRd
ancom_multicomp <- function(phy, var, tax_level) {
  get_comparisons(var, phy, as_list = TRUE, n_cut = 1) %>%
    purrr::map_dfr(function(comparison) {

      tax_table <-
        phyloseq::tax_table(phy) %>%
        to_tibble("taxa_id") %>%
        dplyr::select(taxa_id, taxa = !!tax_level)

      meta <-
        phyloseq::sample_data(phy) %>%
        to_tibble("sample_id") %>%
        dplyr::select(sample_id, !!var) %>%
        dplyr::filter(!!sym(var) %in% comparison)

      ancom_df <-
        phyloseq::otu_table(phy) %>%
        .[, meta$sample_id] %>%
        as.data.frame.matrix()

      prepro <- ANCOM::feature_table_pre_process(
        feature_table = ancom_df,
        meta_data = meta,
        sample_var = 'sample_id',
        lib_cut = 0,
        out_cut = 0,
        zero_cut = 1,
        neg_lb = TRUE
      )

      ancom_res <- ANCOM::ANCOM(
        prepro$feature_table,
        prepro$meta_data,
        prepro$structure_zeros,
        main_var = var
      )

      n_taxa <- ifelse(
        is.null(prepro$struc_zero),
        nrow(prepro$feature_table),
        sum(apply(prepro$struc_zero, 1, sum) == 0)
      )

      ancom_res$out %>%
        tibble::as_tibble() %>%
        dplyr::mutate(
          clr_mean_diff = ancom_res$fig$data$x,
          structural_zero = ancom_res$fig$data$zero_ind,
          x_label = min(ancom_res$fig$data$x),
          y_label = 0.7 * (n_taxa - 1),
          comparison = stringr::str_c(comparison, collapse = "_"),
          var = var
        ) %>%
        dplyr::left_join(tax_table, by = "taxa_id")
    })
}

#' Plot ANCOM results
#'
#' @param dar Result of ANCOM
#' @param var Categorical variable
#' @param pct
#' @param no_structural
#'
#' @return ggplot2 object
#' @noRd
plot_ancom <- function(dar, var, pct = 0.7, no_structural = TRUE) {

  plt <- NULL
  res_slot <- dar$ANCOM %>% purrr::pluck(var)
  if (not_null(res_slot)) {
    dat <-
      res_slot %>%
      dplyr::mutate(
        significant = !!sym(glue::glue("detected_{pct}")),
        structural_zero = factor(structural_zero)
      )

    if (no_structural) {
      dat <-
        dat %>%
        dplyr::mutate(significant = dplyr::if_else(structural_zero == "Yes", FALSE, significant))
    }

    p <-
      dat %>%
      ggplot(aes(clr_mean_diff, W, label = taxa)) +
      facet_wrap( ~ comparison) +
      geom_point(aes(alpha = structural_zero, colour = significant)) +
      geom_hline(aes(yintercept = y_label), linetype = "dashed") +
      geom_text(
        data = . %>% dplyr::filter(!duplicated(y_label)),
        aes(min(x_label), y_label, label = "W[0.7]"),
        size = 4,
        vjust = -0.5,
        hjust = -1,
        color = "orange",
        parse = TRUE
      ) +
      scale_colour_manual(values = c("gray", "lightblue")) +
      scale_alpha_manual(values = c(1, 0.2)) +
      labs(x = "CLR mean difference", y = "W statistic") +
      theme(legend.position = "none")

    plt <-
      plotly::ggplotly(p, tooltip = c("x", "colour", "alpha", "taxa")) %>%
      plotly::config(displaylogo = F)
  }
  plt
}
