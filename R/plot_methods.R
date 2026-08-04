#' Plot otuput of the `overlap_df` function as a heatmap. 
#'
#' @param rec A prepped `Recipe` object.
#' @param steps Character vector with step_ids to take in account.
#' @param font_size Size of the axis font.
#' @param type Indicates whether to use all taxa ("all") or only those that are
#'   differentially abundant in at least one method ("da"). Default as "all". 
#' @param target Optional modeled target to include.
#' @param contrast_id Optional modeled contrast identifier to include.
#'
#' @aliases corr_heatmap
#' @importFrom heatmaply heatmaply_cor
#' @return heatmap
#' @export
#' @autoglobal
#' @tests
#' data(test_rec)
#' data(test_prep_rec)
#' expect_error(corr_heatmap(test_rec))
#' expect_s3_class(corr_heatmap(test_prep_rec), "plotly")
#' expect_s3_class(
#'   corr_heatmap(test_prep_rec, steps = steps_ids(test_prep_rec, "da")), 
#'   "plotly"
#' )
#' expect_s3_class(corr_heatmap(test_prep_rec, type = "da"), "plotly")
#' @examples
#' data(test_prep_rec)
#'
#' ## Running the function returns a UpSet plot ordered by frequency.
#' corr_heatmap(test_prep_rec)
#'
#' ## corr_heatmap function needs a PrepRecipe. If you pass a a non-prep
#' ## Recipe the output is an error.
#' data(test_rec)
#' err <- testthat::expect_error(corr_heatmap(test_rec))
#' err
corr_heatmap <- function(rec,
                         steps = steps_ids(rec, "da"),
                         font_size = 15,
                         type = "all",
                         target = NULL,
                         contrast_id = NULL) {
  
  check_prep_recipe(rec)
  
  overlap_df(
    rec, steps = steps, type = type, target = target,
    contrast_id = contrast_id
  ) %>%
    heatmaply::heatmaply_cor(
      x = .,
      point_size_mat = .,
      colors = c(
        "#006837", "#1A9850", "#66BD63", "#A6D96A", "#D9EF8B",
        "#FFFFBF", "#FEE08B", "#FDAE61", "#F46D43", "#D73027", "#A50026"
      ),
      limits = c(0, 1),
      node_type = "scatter",
      trace = "none",
      dist_method = "canberra", 
      fontsize_col = font_size,
      fontsize_row = font_size,
      heatmap_layers = theme(
        axis.text = element_text(colour = "black", family = 'Arial')
      ), 
      hclust_method = "complete"
    ) 
}

## Intersection plot----

#' Plot results using UpSet plot
#'
#' @param rec A prepped `Recipe` object.
#' @param steps Character vector with step_ids to take in account.
#' @param ordered_by How the intersections in the matrix should be ordered by.
#'   Options include frequency (entered as "freq"), degree, or both in any
#'   order.
#' @param font_size Size of the font. 
#' @param target Optional modeled target to include.
#' @param contrast_id Optional modeled contrast identifier to include.
#'
#' @aliases intersection_plt
#' @return UpSet plot
#' @export
#' @autoglobal
#' @importFrom UpSetR upset
#' @tests
#' data(test_rec)
#' data(test_prep_rec)
#' expect_error(intersection_plt(test_rec))
#' expect_s3_class(intersection_plt(test_prep_rec), "upset")
#' expect_s3_class(
#'   intersection_plt(test_prep_rec, steps = steps_ids(test_prep_rec, "da")), 
#'   "upset"
#' )
#' expect_s3_class(intersection_plt(test_prep_rec, ordered_by = "degree"), "upset")
#' @examples
#' data(test_prep_rec)
#'
#' ## Running the function returns a UpSet plot ordered by frequency.
#' intersection_plt(test_prep_rec)
#'
#' ## Alternatively, you can order the plot by degree
#' intersection_plt(test_prep_rec, ordered_by = "degree")
#'
#' ## intersection_plt function needs a PrepRecipe. If you pass a a non-prep
#' ## Recipe the output is an error.
#' data(test_rec)
#' err <- testthat::expect_error(intersection_plt(test_rec))
#' err
intersection_plt <- function(rec,
                             steps = steps_ids(rec, "da"),
                             ordered_by = c("freq", "degree"), 
                             font_size = 2,
                             target = NULL,
                             contrast_id = NULL) {
  
  check_prep_recipe(rec)
  UpSetR::upset(
    data = intersection_df(
      rec, steps, target = target, contrast_id = contrast_id
    ),
    sets = steps,
    sets.bar.color = "#56B4E9",
    order.by = ordered_by, 
    text.scale = font_size
  )
}

## Exclusion plot----

#' Plot the number of shared DA OTUs between methods.
#'
#' @param rec A prepped `Recipe` object.
#' @param steps Character vector with step_ids to take in account.
#' @param target Optional modeled target to include.
#' @param contrast_id Optional modeled contrast identifier to include.
#'
#' @aliases exclusion_plt
#' @return ggplot2-class object
#' @export
#' @autoglobal
#' @import ggplot2
#' @tests
#' data(test_rec)
#' data(test_prep_rec)
#' 
#' expect_error(exclusion_plt(test_rec))
#' expect_s3_class(exclusion_plt(test_prep_rec), "ggplot")
#' expect_s3_class(
#'   exclusion_plt(test_prep_rec, steps = steps_ids(test_prep_rec, "da")[-1]),
#'   "ggplot"
#' )
#' @examples
#' data(test_prep_rec)
#'
#' ## Running the function returns a barplot plot,
#' exclusion_plt(test_prep_rec)
#'
#' ## If you want to exclude a method for the plot, you can remove it with the
#' ## step parameter. In the following example we eliminate from the graph the
#' ## results of maaslin
#' exclusion_plt(test_prep_rec, steps = steps_ids(test_prep_rec, "da")[-1])
#'
#' ## exclusion_plt function needs a PrepRecipe. If you pass a a non-prep
#' ## Recipe the output is an error.
#' data(test_rec)
#' err <- testthat::expect_error(exclusion_plt(test_rec))
#' err
exclusion_plt <- function(rec, steps = steps_ids(rec, "da"), target = NULL,
                          contrast_id = NULL) {
  
  check_prep_recipe(rec)

  intersections <- intersection_df(
    rec, steps = steps, target = target, contrast_id = contrast_id
  )
  key_columns <- setdiff(names(intersections), steps)
  intersections <- intersections %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(steps),
      names_to = "name",
      values_to = "value"
    )

  df <- steps %>%
    purrr::map_dfr(function(step_id) {
      retained_keys <- intersections %>%
        dplyr::filter(.data$name == .env$step_id, .data$value == 1L) %>%
        dplyr::distinct(dplyr::across(dplyr::all_of(key_columns)))

      intersections %>%
        dplyr::semi_join(retained_keys, by = key_columns) %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(key_columns))) %>%
        dplyr::summarise(sum = sum(.data$value), .groups = "drop") %>%
        dplyr::count(.data$sum) %>%
        dplyr::mutate(method = .env$step_id, total = sum(.data$n))
    })
  
  df %>%
    ggplot(aes(
      x = stats::reorder(method, total),
      y = n,
      fill = factor(sum)
    )) +
    geom_bar(stat = "identity", alpha = 0.9) +
    scale_fill_brewer(palette = "Spectral", direction = -1) +
    coord_flip() +
    theme_minimal() +
    labs(
      y = "Total number of differentially Abundant OTUs",
      x = "method identifier",
      fill = "Shared"
    )
}

## Abundance plot----

#' Abundance boxplot
#'
#' @param rec A prepped `Recipe` object.
#' @param taxa_ids Character vector with taxa_ids to plot. If taxa_ids is NULL
#'   the significant characteristics present in all of the executed methods will
#'   be plotted.
#' @param type Character vector indicating the type of the result. Options:
#'   c("boxplot", "heatmap").
#' @param transform Transformation to apply. The options include:
#'   'compositional' (ie relative abundance), 'Z', 'log10', 'log10p',
#'   'hellinger', 'identity', 'clr', 'alr', or any method from the
#'   vegan::decostand function. If the value is NULL, no normalization is
#'   applied and works with the raw counts.
#' @param scale Scaling constant for the abundance values when transform =
#'   "scale".
#' @param top_n Maximum number of taxa to represent. Default: 20.
#' @param target Modeled target used for annotations. Required when a modeled
#'   recipe contains more than one target.
#' @param contrast_id Modeled contrast used for automatic taxon selection.
#'   Required when `taxa_ids = NULL` and more than one contrast is available.
#'
#' @return ggplot2 or HeatmapList
#' @export
#' @autoglobal
#' @tests
#' data(test_rec)
#' data(test_prep_rec)
#' taxa_ids <- c("Otu_96", "Otu_78", "Otu_88", "Otu_35", "Otu_94", "Otu_34")
#' 
#' pdf(file = tempfile())
#' withr::defer(dev.off())
#' 
#' expect_error(abundance_plt(test_rec))
#' expect_s3_class(suppressMessages(
#'   abundance_plt(test_prep_rec)), "ggplot"
#' )
#' expect_s4_class(suppressMessages(
#'   abundance_plt(test_prep_rec, type = "heatmap")), "HeatmapList"
#' )
#' expect_s3_class(suppressMessages(
#'   abundance_plt(test_prep_rec, taxa_ids = taxa_ids)), "ggplot"
#' )
#' expect_s3_class(suppressMessages(
#'   abundance_plt(test_prep_rec, top_n = 10)), "ggplot"
#' )
#' expect_s3_class(suppressMessages(
#'   abundance_plt(test_prep_rec, transform = "clr")), "ggplot"
#' )
#' expect_s4_class(
#'   abundance_plt(test_prep_rec, taxa_ids = taxa_ids, type = "heatmap"), 
#'   "HeatmapList"
#' )
#' expect_s3_class(suppressMessages(
#'   abundance_plt(test_prep_rec, transform = "scale", scale = 10)), "ggplot"
#' )
#' 
#' withr::defer(dev.off())
#' @examples
#' data(test_prep_rec)
#'
#' ## Running the function returns a boxplot,
#' abundance_plt(test_prep_rec)
#'
#' ## Giving the value "heatmap" to the type parameter, the resulting graph
#' ## a heatmap.
#' # abundance_plt(test_prep_rec, type = "heatmap")
#'
#' ## By default, those taxa significant in all methods are plotted. If you want
#' ## to graph some determined features, you can pass them as vector through the
#' ## taxa_ids parameter.
#' # taxa_ids <- c("Otu_96", "Otu_78", "Otu_88", "Otu_35", "Otu_94", "Otu_34")
#' # abundance_plt(test_prep_rec, taxa_ids = taxa_ids)
#' # abundance_plt(test_prep_rec, taxa_ids = taxa_ids, type = "heatmap")
#'
#' ## abundance_plt function needs a PrepRecipe. If you pass a a non-prep
#' ## Recipe the output is an error.
#' data(test_rec)
#' err <- testthat::expect_error(abundance_plt(test_rec))
#' err
abundance_plt <- function(rec,
                          taxa_ids = NULL,
                          type = "boxplot",
                          transform = "compositional",
                          scale = 1, 
                          top_n = 20,
                          target = NULL,
                          contrast_id = NULL)  {
  
  check_prep_recipe(rec)
  
  if (!type %in% c("boxplot", "heatmap")) {
    cli::cli_abort(
      c("x" = "{.arg type} must be one of: {.val {c('boxplot', 'heatmap')}}."),
      class = "dar_error_invalid_plot_type"
    )
  }

  target <- resolve_plot_target(rec, target)
  contrast_id <- resolve_abundance_contrast(rec, taxa_ids, target, contrast_id)

  if (type == "boxplot") { 
    plt <- .abundance_boxplot(
      rec, taxa_ids, transform, scale, top_n, target, contrast_id
    )
  } else if (type == "heatmap") { 
    plt <- .abundance_heatmap(
      rec, taxa_ids, transform, scale, top_n, target, contrast_id
    )
  }
  
  plot(plt)
}

## Mutual finding plot----

#' Mutual finding plot
#'
#' Plots number of differentially abundant features mutually found by defined
#' number of methods, colored by the differential abundance direction and
#' separated by comparison.
#'
#' @param rec A prepped `Recipe` object.
#' @param count_cutoff Indicates the minimum number of methods in which an OTU
#'   must be present (Default: NULL). If count_cutoff is NULL count_cutoff is
#'   equal to `length(steps_ids(rec, "da")) * 2 / 3`.
#' @param comparisons By default, this function plots all comparisons. However,
#'   if the user indicates the comparison or comparisons of interest, only the
#'   selected ones will be plotted.
#' @param steps Character vector with step_ids to take in account. Default all
#'   "da" methods.
#' @param top_n Maximum number of taxa to represent. Default: 20.
#' @param target Optional modeled target to include.
#' @param contrast_id Optional modeled contrast identifier to include.
#'
#' @return ggplot2
#' @export
#' @autoglobal
#' @tests
#' data(test_rec)
#' data(test_prep_rec)
#' 
#' expect_error(mutual_plt(test_rec))
#' expect_s3_class(suppressMessages(mutual_plt(test_prep_rec)), "ggplot")
#' expect_s3_class(suppressMessages(
#'   mutual_plt(test_prep_rec, count_cutoff = 2)), "ggplot"
#' )
#' expect_error(suppressMessages(mutual_plt(test_prep_rec, count_cutoff = 0)))
#' expect_error(suppressMessages(mutual_plt(test_prep_rec, count_cutoff = 10)))
#' expect_s3_class(suppressMessages(
#'   mutual_plt(test_prep_rec, comparisons = "hts_msm")), "ggplot"
#' )
#' expect_error(suppressMessages(
#'   mutual_plt(test_prep_rec, comparisons = "hts_sm"))
#' )
#' expect_s3_class(suppressMessages(
#'   mutual_plt(test_prep_rec, steps = steps_ids(test_prep_rec, "da")[-1])), 
#'   "ggplot"
#' )
#' expect_error(suppressMessages(mutual_plt(test_prep_rec, steps = "da")))
#' expect_s3_class(suppressMessages(
#'   mutual_plt(test_prep_rec, top_n = 10)), "ggplot"
#' )
#' expect_error(suppressMessages(mutual_plt(test_prep_rec, top_n = 0)))
#' @examples
#' data(test_prep_rec)
#'
#' ## Running the function returns a tile plot,
#' mutual_plt(test_prep_rec)
#'
#' ## The count_cutoff indicates the minimum number of methods in which an OTU
#' ## must be present. By default the value is equal to
#' ## length(steps_ids(rec, "da")) * 2 / 3 but it is customizable.
#' mutual_plt(
#'   test_prep_rec, 
#'   count_cutoff = length(steps_ids(test_prep_rec, "da"))
#' )
#'
#' ## A single comparisons can be plotted through the comparison parameter.
#' mutual_plt(test_prep_rec, comparisons = c("hts_msm"))
#'
#' ## If you want to exclude a method for the plot, you can remove it with the
#' ## step parameter. In the following example we eliminate from the graph the
#' ## results of maaslin.
#' mutual_plt(test_prep_rec, steps = steps_ids(test_prep_rec, "da")[-1])
#'
#' ## mutual_plt function needs a PrepRecipe. If you pass a a non-PrepRecipe
#' ## the output is an error.
#' data(test_rec)
#' err <- testthat::expect_error(mutual_plt(test_rec))
#' err
mutual_plt <- function(rec,
                       count_cutoff = NULL,
                       comparisons = NULL,
                       steps = steps_ids(rec, type = "da"),
                       top_n = 20,
                       target = NULL,
                       contrast_id = NULL) {
  
  check_prep_recipe(rec)
  
  validate_consensus_steps(rec, steps)
  if (!is.numeric(top_n) || length(top_n) != 1L || !is.finite(top_n) ||
      top_n <= 0) {
    cli::cli_abort(
      c("x" = "{.arg top_n} must be greater than {.val {0}}."),
      class = "dar_error_invalid_top_n"
    )
  }

  if (is.null(count_cutoff)) {
    count_cutoff <- steps %>% length() * 2 / 3 
    count_cutoff <- round(count_cutoff, 0)
    cli::cli_inform(
      c("i" = "{.arg count_cutoff} set to {.val {count_cutoff}}.")
    )
  }

  if (!is.numeric(count_cutoff) || length(count_cutoff) != 1L ||
      !is.finite(count_cutoff) || count_cutoff > length(steps)) {
    cli::cli_abort(
      c("x" = "{.arg count_cutoff} must be less than or equal to the number of methods."),
      class = "dar_error_invalid_count_cutoff"
    )
  }

  if (count_cutoff <= 0) {
    cli::cli_abort(
      c("x" = "{.arg count_cutoff} must be greater than {.val {0}}."),
      class = "dar_error_invalid_count_cutoff"
    )
  }

  all_significant <- .all_significant(rec, steps = steps)
  if (!is.null(get_model(rec))) {
    all_significant <- filter_model_hypotheses(
      all_significant, target = target, contrast_id = contrast_id
    )
  }
  df <- dplyr::filter(
    all_significant, .data$method_count >= .env$count_cutoff
  )
  hypothesis_keys <- if (is.null(get_model(rec))) {
    "taxa_id"
  } else {
    c("taxa_id", "contrast_id", "effect")
  }
  ranked <- all_significant %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(hypothesis_keys))) %>%
    dplyr::summarise(method_count = max(.data$method_count), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(.data$method_count))
  
  if (nrow(dplyr::distinct(
    df, dplyr::across(dplyr::all_of(hypothesis_keys))
  )) > top_n) {
    cli::cli_inform(c(
      "!" = "Taxa present in selected methods are greater than the cutoff {.arg top_n} = {.val {top_n}}.",
      "i" = "The top {.val {top_n}} significant taxa with the greatest overlap between methods will be used."
    ))
    
    selected <- ranked %>%
      dplyr::filter(.data$method_count >= .env$count_cutoff) %>%
      dplyr::slice_head(n = top_n) %>%
      dplyr::select(dplyr::all_of(hypothesis_keys))
    df <- dplyr::semi_join(df, selected, by = hypothesis_keys)
  }
  
  if (nrow(df) == 0) {
    cli::cli_inform(c(
      "!" = "0 taxa are present with {.arg count_cutoff} = {.val {count_cutoff}}.",
      "i" = "The top {.val {top_n}} significant taxa with the greatest overlap between methods will be used."
    ))
    
    selected <- ranked %>%
      dplyr::slice_head(n = top_n) %>%
      dplyr::select(dplyr::all_of(hypothesis_keys))
    df <- dplyr::semi_join(all_significant, selected, by = hypothesis_keys)
  }
  
  if (!is.null(comparisons)) {
    if (!all(comparisons %in% df$comparison)) {
      cli::cli_abort(c(
        "x" = "Invalid value in {.arg comparisons}.",
        "i" = "Valid comparisons are: {.val {unique(df$comparison)}}."
      ),
      class = "dar_error_invalid_comparison"
      )
    }
    df <- df %>% dplyr::filter(comparison %in% comparisons)
  }

  facet <- if (is.null(get_model(rec))) ~ comparison else ~ var + comparison

  df %>% 
    dplyr::left_join(
      .all_stats(rec),
      by = c("taxa_id", "contrast_id", "var", "comparison", "method")
    ) %>% 
    dplyr::mutate(
      method = stringr::str_remove_all(method, "[:alpha:]_[:alpha:].*")
    ) %>% 
    tidyr::unite("taxa", c(taxa_id, taxa), sep = "|") %>% 
    dplyr::group_by(method) %>% 
    dplyr::mutate(zscore = scales::rescale(effect_v)) %>% 
    ggplot(aes(taxa, method, fill = effect)) +
    geom_tile(width = 0.7, height = 0.8, alpha = 1) +
    facet_wrap(facet, ncol = 1, strip.position = "right") +
    theme_light() +
    scale_fill_manual(values = c(down = "#74ADD1", neutral = "grey70", up = "#F46D43")) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 0.9)) +
    labs(x = NULL, y = NULL, fill = "DA")
}

# --- INTERNAL HELPERS ---

#' @noRd
#' @keywords internal
#' @autoglobal
.annotated_counts <- function(rec) {
  otu_table(rec) %>% 
    tidyr::pivot_longer(-1, names_to = "sample_id") %>%
    dplyr::left_join(analysis_sample_data(rec), by = "sample_id") %>%
    dplyr::left_join(analysis_tax_table(rec), by = "taxa_id")
}

#' @noRd
resolve_plot_target <- function(rec, target = NULL) {
  available <- recipe_targets(rec)
  if (is.null(get_model(rec))) return(target %||% available[[1]])
  if (is.null(target) && length(available) != 1L) {
    cli::cli_abort(
      "{.arg target} is required when the model contains multiple targets.",
      class = "dar_error_ambiguous_plot_target"
    )
  }
  if (is.null(target)) target <- available
  if (!is.character(target) || length(target) != 1L || is.na(target) ||
      !nzchar(target) || !target %in% available) {
    cli::cli_abort(
      "{.arg target} must identify one modeled target: {.val {available}}.",
      class = "dar_error_invalid_contrast_selector"
    )
  }
  target
}

#' @noRd
resolve_abundance_contrast <- function(rec, taxa_ids, target,
                                       contrast_id = NULL) {
  if (is.null(get_model(rec))) return(NULL)
  plan <- resolve_model(rec)$contrast_plan
  available <- plan$contrast_id[plan$var %in% target]
  if (!is.null(taxa_ids) && is.null(contrast_id)) return(NULL)
  if (is.null(contrast_id) && length(available) != 1L) {
    cli::cli_abort(
      paste0(
        "{.arg contrast_id} is required for automatic taxon selection when ",
        "multiple modeled contrasts are available."
      ),
      class = "dar_error_ambiguous_plot_contrast"
    )
  }
  if (is.null(contrast_id)) contrast_id <- available
  if (!is.character(contrast_id) || length(contrast_id) != 1L ||
      is.na(contrast_id) || !nzchar(contrast_id) ||
      !contrast_id %in% available) {
    cli::cli_abort(
      "{.arg contrast_id} must identify one planned model contrast.",
      class = "dar_error_invalid_contrast_selector"
    )
  }
  contrast_id
}

#' @noRd
selected_significant <- function(rec, target = NULL, contrast_id = NULL,
                                 steps = steps_ids(rec, "da")) {
  data <- .all_significant(rec, steps = steps)
  if (!is.null(get_model(rec))) {
    data <- filter_model_hypotheses(
      data, target = target, contrast_id = contrast_id
    )
  }
  data
}

#' @noRd
#' @keywords internal
#' @autoglobal
.abundance_boxplot <- function(rec, taxa_ids, transform, scale, top_n = 20,
                               target, contrast_id) {
  if (is.null(taxa_ids)) {
    taxa_ids <- 
      selected_significant(rec, target, contrast_id) %>%
      dplyr::filter(method_count == length(steps_ids(rec, "da"))) %>% 
      dplyr::pull(taxa_id) %>% 
      unique()
    
    if (length(taxa_ids) > top_n) {
      cli::cli_inform(c(
        "!" = "Taxa present in all methods are greater than the cutoff {.arg top_n} = {.val {top_n}}.",
        "i" = "The top {.val {top_n}} significant taxa with the greatest overlap between methods will be used."
      ))
     
      taxa_ids <- 
        selected_significant(rec, target, contrast_id) %>%
        dplyr::group_by(taxa_id) %>% 
        dplyr::summarise(method_count = max(method_count)) %>% 
        dplyr::arrange(-method_count) %>% 
        dplyr::pull(taxa_id) %>% 
        .[seq_len(top_n)]
    }
    
    if (length(taxa_ids) == 0) {
      cli::cli_inform(c(
        "!" = "0 taxa are present in all tested methods.",
        "i" = "The top {.val {top_n}} significant taxa with the greatest overlap between methods will be used."
      ))
      
      taxa_ids <- 
        selected_significant(rec, target, contrast_id) %>%
        dplyr::group_by(taxa_id) %>% 
        dplyr::summarise(method_count = max(method_count)) %>% 
        dplyr::arrange(-method_count) %>% 
        dplyr::pull(taxa_id) %>% 
        .[seq_len(top_n)]
    }
  }
  
  t_rec <- rec
  if (!is.null(transform)) {
    t_rec@phyloseq <- microbiome::transform(
      t_rec@phyloseq, 
      transform = transform, 
      scale = scale
    )
  } else {
    transform <- "raw counts"
  }
  
  .annotated_counts(t_rec) %>% 
    dplyr::filter(taxa_id %in% taxa_ids) %>% 
    tidyr::unite("taxa", c(taxa_id, taxa), sep = "|") %>% 
    ggplot(aes(taxa, value, fill = !!dplyr::sym(target))) +
    geom_boxplot(alpha = 0.7) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 0.9)) +
    labs(x = NULL, y = paste0("Abundance (", transform, ")"))
}

#' @noRd
#' @keywords internal
#' @autoglobal
.abundance_heatmap <- function(rec, taxa_ids, transform, scale, top_n, target,
                               contrast_id) {
  ComplexHeatmap::ht_opt(
    message = FALSE, 
    COLUMN_ANNO_PADDING = unit(0.5, "cm")
  )
  
  if (is.null(taxa_ids)) {
    taxa_ids <- 
      selected_significant(rec, target, contrast_id) %>%
      dplyr::filter(method_count == length(steps_ids(rec, "da"))) %>% 
      dplyr::pull(taxa_id) %>% 
      unique()
    
    if (length(taxa_ids) > top_n) {
      cli::cli_inform(c(
        "!" = "Taxa present in all methods are greater than the cutoff {.arg top_n} = {.val {top_n}}.",
        "i" = "The top {.val {top_n}} significant taxa with the greatest overlap between methods will be used."
      ))
      
      taxa_ids <- 
        selected_significant(rec, target, contrast_id) %>%
        dplyr::group_by(taxa_id) %>% 
        dplyr::summarise(method_count = max(method_count)) %>% 
        dplyr::arrange(-method_count) %>% 
        dplyr::pull(taxa_id) %>% 
        .[seq_len(top_n)]
    }
    
    if (length(taxa_ids) == 0) {
      cli::cli_inform(c(
        "!" = "{.val 0} taxa are present in all tested methods.",
        "i" = "The top {.val {top_n}} significant taxa with the greatest overlap between methods will be used."
      ))
      
      taxa_ids <- 
        selected_significant(rec, target, contrast_id) %>%
        dplyr::group_by(taxa_id) %>% 
        dplyr::summarise(method_count = max(method_count)) %>% 
        dplyr::arrange(-method_count) %>% 
        dplyr::pull(taxa_id) %>% 
        .[seq_len(top_n)]
    }
  }
  
  t_rec <- rec
  if (!is.null(transform)) {
    t_rec@phyloseq <- 
      microbiome::transform(t_rec@phyloseq, transform, scale = scale)
  } else {
    transform <- "raw counts"
  }
    
  df <- 
    .annotated_counts(t_rec) %>% 
    dplyr::filter(taxa_id %in% taxa_ids) %>% 
    tidyr::unite("taxa", c(taxa_id, taxa), sep = "|")
  
  mtx <- 
    dplyr::select(df, taxa, value, sample_id) %>% 
    tidyr::pivot_wider(names_from = sample_id, values_from = value) %>% 
    data.frame(row.names = 1) %>% 
    as.matrix()
  
  annot <- 
    dplyr::select(df, sample_id, dplyr::all_of(target)) %>%
    dplyr::distinct() %>% 
    data.frame(row.names = 1) %>% 
    ComplexHeatmap::HeatmapAnnotation(
      df = ., 
      simple_anno_size = grid::unit(0.4, "cm"),
      annotation_name_gp = grid::gpar(fontsize = 8, fontface = 'bold'),
      annotation_name_side = 'right'
    )
  
  ComplexHeatmap::Heatmap(
    matrix = mtx,
    clustering_distance_columns = "euclidean",
    clustering_method_columns = 'ward.D2',
    cluster_rows = TRUE,
    show_row_dend = FALSE,
    col = gplots::colorpanel(1e4, "white", "steelblue", "darkred"),
    top_annotation = annot,
    row_names_gp = grid::gpar(fontsize = 8),
    show_column_names = FALSE,
    name = paste0("Abundance (", transform, ")"),
    show_heatmap_legend = ifelse(length(annot@anno_list) > 5, FALSE, TRUE)
  )
}

#' @noRd
#' @keywords internal
#' @autoglobal
.otu_method_count <- function(rec, steps = steps_ids(rec, "da")) {
  intersection_df(rec, steps = steps, tidy = TRUE) %>%
    dplyr::group_by(taxa_id) %>%
    dplyr::summarise(method_count = sum(value))
}

#' @noRd
effect_direction <- function(effect) {
  dplyr::case_when(
    effect > 0 ~ "up",
    effect < 0 ~ "down",
    TRUE ~ "neutral"
  )
}

#' @noRd
#' @keywords internal
#' @autoglobal
.otu_effect_direction <- function(rec, steps = steps_ids(rec, "da")) {
  tidy_results(rec, steps = steps, significant_only = TRUE) %>%
    dplyr::transmute(
      taxa_id = .data$taxa_id,
      taxa = .data$taxa,
      comparison = .data$comparison,
      effect = effect_direction(.data$effect_size),
      method = .data$step_id,
      contrast_id = .data$contrast_id,
      contrast_type = .data$contrast_type,
      var = .data$var
    )
}

#' @noRd
#' @keywords internal
#' @autoglobal
.all_significant <- function(rec, steps = steps_ids(rec, "da")) {
  validate_consensus_steps(rec, steps)
  if (!is.null(get_model(rec))) {
    return(
      .otu_effect_direction(rec, steps = steps) %>%
        dplyr::group_by(
          .data$taxa_id, .data$taxa, .data$contrast_id,
          .data$comparison, .data$contrast_type, .data$var, .data$effect
        ) %>%
        dplyr::mutate(method_count = dplyr::n_distinct(.data$method)) %>%
        dplyr::ungroup()
    )
  }
  dplyr::left_join(
    .otu_effect_direction(rec, steps = steps),
    .otu_method_count(rec, steps = steps),
    by = "taxa_id"
  )
}

#' @noRd
#' @keywords internal
#' @autoglobal
.all_stats <- function(rec) {
  tidy_results(rec) %>%
    dplyr::transmute(
      taxa_id = .data$taxa_id,
      contrast_id = .data$contrast_id,
      var = .data$var,
      comparison = .data$comparison,
      effect_v = .data$effect_size,
      padj = .data$adj_p_value,
      method = .data$step_id
    )
}
