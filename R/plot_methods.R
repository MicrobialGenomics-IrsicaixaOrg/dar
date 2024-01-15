## correlation plot----

#' Plot otuput of the `overlap_df` function as a heatmap. 
#'
#' @param rec A `Recipe` object.
#' @param steps Character vector with step_ids to take in account.
#' @param font_size Size of the axis font.
#' @param type Indicates whether to use all taxa ("all") or only those that are
#'   differentially abundant in at least one method ("da"). Default as "all". 
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
#'   corr_heatmap(test_prep_rec, steps = steps_ids(test_prep_rec, "da")[-1]), 
#'   "plotly"
#' )
#' expect_s3_class(corr_heatmap(test_prep_rec, type = "da"), "plotly")
#' @examples
#' data(test_prep_rec)
#'
#' ## Running the function returns a UpSet plot ordered by frequency.
#' corr_heatmap(test_prep_rec)
#'
#' ## If you want to exclude a method for the plot, you can remove it with the
#' ## step parameter. In the following example we eliminate from the graph the
#' ## results of maaslin
#' corr_heatmap(test_prep_rec, steps = steps_ids(test_prep_rec, "da")[-1])
#'
#' ## corr_heatmap function needs a PrepRecipe. If you pass a a non-prep
#' ## Recipe the output is an error.
#' data(test_rec)
#' err <- testthat::expect_error(corr_heatmap(test_rec))
#' err
methods::setGeneric(
  name = "corr_heatmap",
  def = function(rec,
                 steps = steps_ids(rec, "da"),
                 font_size = 15,
                 type = "all") {
    standardGeneric("corr_heatmap")
  }
)

#' @rdname corr_heatmap
#' @export
methods::setMethod(
  f = "corr_heatmap",
  signature = "Recipe",
  definition = function(rec, steps, font_size, type) {
    rlang::abort(c(
      "This function needs a PrepRecipe!",
      glue::glue(
        "Run {crayon::bgMagenta('prep(rec)')} and then try with ", 
        "{crayon::bgMagenta('corr_heatmap()')}"
      )
    ))
  }
)

#' @rdname corr_heatmap
#' @export
#' @autoglobal
#' @importFrom UpSetR upset
methods::setMethod(
  f = "corr_heatmap",
  signature = "PrepRecipe",
  definition = function(rec, steps, font_size, type) {
    plt <- overlap_df(rec, steps = steps, type = type) %>%
      heatmaply::heatmaply_cor(
        x = .,
        point_size_mat = .,
        colors = c(
          "#006837",
          "#1A9850",
          "#66BD63",
          "#A6D96A",
          "#D9EF8B",
          "#FFFFBF",
          "#FEE08B",
          "#FDAE61",
          "#F46D43",
          "#D73027",
          "#A50026"
        ),
        limits = c(0, 1),
        node_type = "scatter",
        trace = "none",
        dist_method = "canberra", 
        fontsize_col = font_size,
        fontsize_row = font_size,
        heatmap_layers = theme(axis.text = element_text(
          colour = "black", family = 'Arial')
        ), 
        hclust_method = "complete"
      ) 
    plt
  }
)

## Intersection plot----

#' Plot results using UpSet plot
#'
#' @param rec A `Recipe` object.
#' @param steps Character vector with step_ids to take in account.
#' @param ordered_by How the intersections in the matrix should be ordered by.
#'   Options include frequency (entered as "freq"), degree, or both in any
#'   order.
#' @param font_size Size of the font. 
#'
#' @aliases intersection_plt
#' @return UpSet plot
#' @export
#' @autoglobal
#' @tests
#' data(test_rec)
#' data(test_prep_rec)
#' expect_error(intersection_plt(test_rec))
#' expect_s3_class(intersection_plt(test_prep_rec), "upset")
#' expect_s3_class(
#'   intersection_plt(test_prep_rec, steps = steps_ids(test_prep_rec, "da")[-1]), 
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
#' ## If you want to exclude a method for the plot, you can remove it with the
#' ## step parameter. In the following example we eliminate from the graph the
#' ## results of maaslin
#' intersection_plt(test_prep_rec, steps = steps_ids(test_prep_rec, "da")[-1])
#'
#' ## intersection_plt function needs a PrepRecipe. If you pass a a non-prep
#' ## Recipe the output is an error.
#' data(test_rec)
#' err <- testthat::expect_error(intersection_plt(test_rec))
#' err
methods::setGeneric(
  name = "intersection_plt",
  def = function(rec,
                 steps = steps_ids(rec, "da"),
                 ordered_by = c("freq", "degree"), 
                 font_size = 2) {
    standardGeneric("intersection_plt")
  }
)

#' @rdname intersection_plt
#' @export
methods::setMethod(
  f = "intersection_plt",
  signature = "Recipe",
  definition = function(rec, steps, font_size) {
    rlang::abort(c(
      "This function needs a PrepRecipe!",
      glue::glue(
        "Run {crayon::bgMagenta('prep(rec)')} and then try with ", 
        "{crayon::bgMagenta('intersection_plt()')}"
      )
    ))
  }
)

#' @rdname intersection_plt
#' @export
#' @autoglobal
#' @importFrom UpSetR upset
methods::setMethod(
  f = "intersection_plt",
  signature = "PrepRecipe",
  definition = function(rec, steps, ordered_by, font_size) {
    UpSetR::upset(
      data = intersection_df(rec, steps),
      sets = steps,
      sets.bar.color = "#56B4E9",
      order.by = ordered_by, 
      text.scale = font_size
    )
  }
)

## Exclusion plot----

#' Plot the number of shared DA OTUs between methods.
#'
#' @param rec A `Recipe` object.
#' @param steps Character vector with step_ids to take in account.
#'
#' @aliases exclusion_plt
#' @return ggplot2-class object
#' @export
#' @autoglobal
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
#' ## intersection_plt function needs a PrepRecipe. If you pass a a non-prep
#' ## Recipe the output is an error.
#' data(test_rec)
#' err <- testthat::expect_error(exclusion_plt(test_rec))
#' err
methods::setGeneric(
  name = "exclusion_plt",
  def = function(rec, steps = steps_ids(rec, "da")) {
    standardGeneric("exclusion_plt")
  }
)

#' @rdname exclusion_plt
#' @export
#' @autoglobal
methods::setMethod(
  f = "exclusion_plt",
  signature = "Recipe",
  definition = function(rec, steps) {
    rlang::abort(c(
      "This function needs a PrepRecipe!",
      glue::glue(
        "Run {crayon::bgMagenta('prep(rec)')} and then try with ", 
        "{crayon::bgMagenta('exclusion_plt()')}"
      )
    ))
  }
)

#' @rdname exclusion_plt
#' @export
#' @autoglobal
#' @import ggplot2
methods::setMethod(
  f = "exclusion_plt",
  signature = "PrepRecipe",
  definition = function(rec, steps) {
    df <-
      steps_ids(rec, "da") %>%
      purrr::map_dfr(~ {
        df <-
          intersection_df(rec, steps = steps) %>%
          tidyr::pivot_longer(cols = -1)
        
        to_retain <-
          df %>%
          dplyr::filter(name == .x & value == 1) %>%
          dplyr::pull(taxa_id)
        
        df %>%
          dplyr::filter(taxa_id %in% to_retain) %>%
          dplyr::group_by(taxa_id) %>%
          dplyr::summarise(sum = sum(value)) %>%
          dplyr::count(sum) %>%
          dplyr::mutate(method = .x, total = sum(n))
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
      labs(y = "Total number of differentially Abundant OTUs",
           x = "method identifier",
           fill = "Shared")
  }
)

## Abundance plot----

#' Abundance boxplot
#'
#' @param rec A Recipe or Recipe step.
#' @param taxa_ids Character vector with taxa_ids to plot. If taxa_ids is NULL
#'   the significant characteristics present in all of the executed methods will
#'   be plotted.
#' @param type Character vector indicating the type of the result. Options:
#'   c("boxoplot", "heatmap").
#' @param transform Transformation to apply. The options include:
#'   'compositional' (ie relative abundance), 'Z', 'log10', 'log10p',
#'   'hellinger', 'identity', 'clr', 'alr', or any method from the
#'   vegan::decostand function. If the value is NULL, no normalization is
#'   applied and works with the raw counts.
#' @param scale Scaling constant for the abundance values when transform =
#'   "scale".
#' @param top_n Maximum number of taxa to represent. Default: 20.
#'
#' @return ggplot2
#' @export
#' @autoglobal
#' @tests
#' data(test_rec)
#' data(test_prep_rec)
#' taxa_ids <- c("Otu_96", "Otu_78", "Otu_88", "Otu_35", "Otu_94", "Otu_34")
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
methods::setGeneric(
  name = "abundance_plt",
  def = function(rec,
                 taxa_ids = NULL,
                 type = "boxplot",
                 transform = "compositional",
                 scale = 1, 
                 top_n = 20)  {
    standardGeneric("abundance_plt")
  }
)

#' @rdname abundance_plt
#' @export
#' @autoglobal
methods::setMethod(
  f = "abundance_plt",
  signature = "Recipe",
  definition = function(rec, taxa_ids, type, transform, scale, top_n) {
    rlang::abort(c(
      "This function needs a PrepRecipe!",
      glue::glue(
        "Run {crayon::bgMagenta('prep(rec)')} and then try with ", 
        "{crayon::bgMagenta('abundance_plt()')}"
      )
    ))
  }
)

#' @rdname abundance_plt
#' @export
#' @autoglobal
methods::setMethod(
  f = "abundance_plt",
  signature = "PrepRecipe",
  definition = function(rec, taxa_ids, type, transform, scale, top_n) {
    if (type == "boxplot") { 
      plt <- .abundance_boxplot(rec, taxa_ids, transform, scale, top_n) 
    }
    if (type == "heatmap") { 
      plt <- .abundance_heatmap(rec, taxa_ids, transform, scale, top_n) 
    }
    if (!type %in% c("boxplot", "heatmap")) {
      rlang::abort("type must be boxplot or heatmap")
    }
    plot(plt)
  }
)

#' @noRd
#' @keywords internal
#' @autoglobal
.annotated_counts <- function(rec) {
  otu_table(rec) %>% 
    tidyr::pivot_longer(-1, names_to = "sample_id") %>%
    dplyr::left_join(sample_data(rec), by = "sample_id") %>% 
    dplyr::left_join(tax_table(rec), by = "taxa_id")
}

#' @noRd
#' @keywords internal
#' @autoglobal
.abundance_boxplot <- function(rec, taxa_ids, transform, scale, top_n = 20) {
  if (is.null(taxa_ids)) {
    taxa_ids <- 
      .all_significant(rec) %>% 
      dplyr::filter(method_count == length(steps_ids(rec, "da"))) %>% 
      dplyr::pull(taxa_id) %>% 
      unique()
    
    if (length(taxa_ids) > top_n) {
      rlang::inform(c(
        "!" = glue::glue(
          "Taxa present in all methods are greater than the cutoff top_n = ", 
          "{top_n}"
        ),
        "i" = glue::glue(
          "The top {top_n} significant taxa with the greatest overlap between", 
          " methods will be used"
        )
      ))
     
      taxa_ids <- 
        .all_significant(rec) %>% 
        dplyr::group_by(taxa_id) %>% 
        dplyr::summarise(method_count = max(method_count)) %>% 
        dplyr::arrange(-method_count) %>% 
        dplyr::pull(taxa_id) %>% 
        .[seq_len(top_n)]
    }
    
    if (length(taxa_ids) == 0) {
      rlang::inform(c(
        "!" = "0 taxa are present in all tested methods",
        "i" = glue::glue(
          "The top {top_n} significant taxa with the greatest overlap ", 
          "between methods will be used"
        )
      ))
      
      taxa_ids <- 
        .all_significant(rec) %>% 
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
    ggplot(aes(taxa, value, fill = !!dplyr::sym(get_var(rec)[[1]]))) +
    geom_boxplot(alpha = 0.7) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 0.9)) +
    labs(x = NULL, y = paste0("Abundance (", transform, ")"))
}

#' @noRd
#' @keywords internal
#' @autoglobal
.abundance_heatmap <- function(rec, taxa_ids, transform, scale, top_n) {
  ComplexHeatmap::ht_opt(
    message = FALSE, 
    COLUMN_ANNO_PADDING = unit(0.5, "cm")
  )
  
  if (is.null(taxa_ids)) {
    taxa_ids <- 
      .all_significant(rec) %>% 
      dplyr::filter(method_count == length(steps_ids(rec, "da"))) %>% 
      dplyr::pull(taxa_id) %>% 
      unique()
    
    if (length(taxa_ids) > top_n) {
      rlang::inform(c(
        "!" = glue::glue(
          "Taxa present in all methods are greater than the cutoff top_n = ", 
          "{top_n}"
        ),
        "i" = glue::glue(
          "The top {top_n} significant taxa with the greatest overlap between", 
          " methods will be used"
        )
      ))
      
      taxa_ids <- 
        .all_significant(rec) %>% 
        dplyr::group_by(taxa_id) %>% 
        dplyr::summarise(method_count = max(method_count)) %>% 
        dplyr::arrange(-method_count) %>% 
        dplyr::pull(taxa_id) %>% 
        .[seq_len(top_n)]
    }
    
    if (length(taxa_ids) == 0) {
      rlang::inform(c(
        "!" = "0 taxa are present in all tested methods",
        "i" = glue::glue(
          "The top {top_n} significant taxa with the greatest overlap between", 
          " methods will be used"
        )
      ))
      
      taxa_ids <- 
        .all_significant(rec) %>% 
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
    dplyr::select(df, sample_id, get_var(rec)[[1]]) %>% 
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

## Mutual finding plot----

#' Mutual finding plot
#'
#' Plots number of differentially abundant features mutually found by defined
#' number of methods, colored by the differential abundance direction and
#' separated by comparison.
#'
#' @param rec A Recipe or Recipe step.
#' @param count_cutoff Indicates the minimum number of methods in which an OTU
#'   must be present (Default: NULL). If count_cutoff is NULL count_cutoff is
#'   equal to `length(steps_ids(rec, "da")) * 2 / 3`.
#' @param comparisons By default, this function plots all comparisons. However,
#'   if the user indicates the comparison or comparisons of interest, only the
#'   selected ones will be plotted.
#' @param steps Character vector with step_ids to take in account. Default all
#'   "da" methods.
#' @param top_n Maximum number of taxa to represent. Default: 20.
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
methods::setGeneric(
  name = "mutual_plt",
  def = function(rec,
                 count_cutoff = NULL,
                 comparisons = NULL,
                 steps = steps_ids(rec, type = "da"),
                 top_n = 20) {
    standardGeneric("mutual_plt")
  }
)

#' @rdname mutual_plt
#' @export
methods::setMethod(
  f = "mutual_plt",
  signature = "Recipe",
  definition = function(rec, count_cutoff, comparisons, steps, top_n) {
    rlang::abort(c(
      "This function needs a PrepRecipe!",
      glue::glue(
        "Run {crayon::bgMagenta('prep(rec)')} and then try with ", 
        "{crayon::bgMagenta('mutual_plt()')}"
      )
    ))
  }
)

#' @rdname mutual_plt
#' @export
#' @autoglobal
methods::setMethod(
  f = "mutual_plt",
  signature = "PrepRecipe",
  definition = function(rec, count_cutoff, comparisons, steps, top_n) {
    
    if (top_n == 0) {
      rlang::abort("top_n must be greater than 0")
    }
    
    if (is.null(count_cutoff)) {
      count_cutoff <- steps %>% length() * 2 / 3 
      count_cutoff <- round(count_cutoff, 0)
      rlang::inform(c("i" = glue::glue("count_cutoff set to {count_cutoff}")))
    }
    
    if (count_cutoff > length(steps_ids(rec, "da"))) {
      rlang::abort("count_cutoff must be less or equal than the number of methods")
    }
    
    if (count_cutoff == 0) {
      rlang::abort("count_cutoff must be greater than 0")
    }
    
    if (!all(steps %in% steps_ids(rec, "da"))) {
      rlang::abort("steps must be a subset of step_ids(rec, 'da')")
    }
    
    df <- 
      .all_significant(rec) %>% 
      dplyr::filter(method_count >= count_cutoff & method %in% steps)
    
    if (nrow(df) > top_n) {
      rlang::inform(c(
        "!" = glue::glue(
          "Taxa present in selected methods are greater than the cutoff top_n = ", 
          "{top_n}"
        ),
        "i" = glue::glue(
          "The top {top_n} significant taxa with the greatest overlap between", 
          " methods will be used"
        )
      ))
      
      taxa_ids <- 
        .all_significant(rec) %>% 
        dplyr::filter(method %in% steps) %>% 
        dplyr::group_by(taxa_id) %>% 
        dplyr::summarise(method_count = max(method_count)) %>% 
        dplyr::arrange(-method_count) %>% 
        dplyr::pull(taxa_id) %>% 
        .[seq_len(top_n)]
      
      df <- .all_significant(rec) %>% dplyr::filter(taxa_id %in% taxa_ids)
    }
    
    if (nrow(df) == 0) {
      rlang::inform(c(
        "!" = glue::glue(
          "0 taxa are present with count_cutoff = {count_cutoff}"
        ),
        "i" = glue::glue(
          "The top {top_n} significant taxa with the greatest overlap between", 
          " methods will be used"
        )
      ))
      
      taxa_ids <- 
        .all_significant(rec) %>% 
        dplyr::filter(method %in% steps) %>% 
        dplyr::group_by(taxa_id) %>% 
        dplyr::summarise(method_count = max(method_count)) %>% 
        dplyr::arrange(-method_count) %>% 
        dplyr::pull(taxa_id) %>% 
        .[seq_len(top_n)]
      
      df <- .all_significant(rec) %>% dplyr::filter(taxa_id %in% taxa_ids)
    }
    
    if (!is.null(comparisons)) {
      if (!all(comparisons %in% df$comparison)) {
        rlang::abort("comparison must be a valid comparison")
      }
      df <- df %>% dplyr::filter(comparison %in% comparisons)
    }

    df %>% 
      dplyr::left_join(
        .all_stats(rec), by = c("taxa_id", "comparison", "method")
      ) %>% 
      dplyr::mutate(
        method = stringr::str_remove_all(method, "[:alpha:]_[:alpha:].*")
      ) %>% 
      tidyr::unite("taxa", c(taxa_id, taxa), sep = "|") %>% 
      dplyr::group_by(method) %>% 
      dplyr::mutate(zscore = scales::rescale(effect_v)) %>% 
      ggplot(aes(taxa, method, fill = effect)) +
      geom_tile(width = 0.7, height = 0.8, alpha = 1) +
      facet_wrap(~ comparison, ncol = 1, strip.position = "right") +
      theme_light() +
      scale_fill_manual(values = c("#74ADD1", "#F46D43")) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 0.9)) +
      labs(x = NULL, y = NULL, fill = "DA")
  }
)

#' @noRd
#' @keywords internal
#' @autoglobal
.otu_method_count <- function(rec) {
  intersection_df(rec, tidy = TRUE) %>%
    dplyr::group_by(taxa_id) %>%
    dplyr::summarise(method_count = sum(value))
}

#' @noRd
#' @keywords internal
#' @autoglobal
.otu_effect_direction <- function(rec) {
  names(rec@results) %>% 
    purrr::map_dfr( ~ {
      rec@results %>% 
        purrr::pluck(.x, 1) %>% 
        dplyr::filter(signif == TRUE) %>% 
        dplyr::mutate(
          effect = dplyr::if_else(effect > 0, "up", "down"), 
          method = .x
        ) %>% 
        dplyr::select(taxa_id, taxa, comparison, effect, method)
    }) 
}

#' @noRd
#' @keywords internal
#' @autoglobal
.all_significant <- function(rec) {
  dplyr::left_join(
    .otu_effect_direction(rec), 
    .otu_method_count(rec),
    by = "taxa_id"
  )
}

#' @noRd
#' @keywords internal
#' @autoglobal
.all_stats <- function(rec) {
  rec@results %>% 
    names() %>% 
    purrr::map_dfr(~ {
      rec@results[[.x]][[1]] %>%
        dplyr::select(taxa_id, comparison, effect_v = effect, dplyr::any_of(c(
          "padj" = "pajd",
          "padj" = "adjp",
          "padj" = "padj"
        ))) %>%
        dplyr::mutate(method = .x)
    })
}
