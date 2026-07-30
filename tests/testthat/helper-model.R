make_longitudinal_phy <- function() {
  metadata <- tidyr::crossing(
    condition = c("control", "treated"),
    subject_index = seq_len(2),
    time = c(0, 1, 2)
  ) |>
    dplyr::mutate(
      subject = paste0(substr(condition, 1, 1), subject_index),
      batch = factor(ifelse(subject_index == 1, "A", "B")),
      age = 30 + subject_index * 5 + ifelse(condition == "treated", 2, 0),
      sample_id = paste(condition, subject_index, time, sep = "_")
    ) |>
    dplyr::select(sample_id, condition, time, subject, batch, age)

  sample_number <- seq_len(nrow(metadata))
  counts <- rbind(
    taxon_1 = 30 + 12 * (metadata$condition == "treated") + 6 * metadata$time,
    taxon_2 = 45 - 8 * (metadata$condition == "treated") + 2 * metadata$time,
    taxon_3 = 20 + sample_number,
    taxon_4 = 18 + (sample_number %% 4),
    taxon_5 = 12 + 3 * metadata$time,
    taxon_6 = 25 + (sample_number %% 3)
  )
  colnames(counts) <- metadata$sample_id

  taxonomy <- matrix(
    paste0("Species_", seq_len(nrow(counts))),
    ncol = 1,
    dimnames = list(rownames(counts), "Species")
  )

  phyloseq::phyloseq(
    phyloseq::otu_table(counts, taxa_are_rows = TRUE),
    phyloseq::tax_table(taxonomy),
    phyloseq::sample_data(data.frame(
      dplyr::select(metadata, -sample_id),
      row.names = metadata$sample_id
    ))
  )
}

model_direction_truth <- c(taxon_1 = 1, taxon_2 = -1)

expect_model_directions <- function(
    result,
    expected = model_direction_truth,
    contrast_id = "condition[treated-control]"
) {
  result <- dar:::flatten_model_output(result)

  purrr::iwalk(expected, function(direction, taxa_id) {
    observed <- result |>
      dplyr::filter(
        .data$taxa_id == .env$taxa_id,
        .data$contrast_id == .env$contrast_id
      )
    context <- paste(taxa_id, contrast_id, sep = " / ")

    testthat::expect_equal(
      nrow(observed), 1L,
      info = paste("Expected exactly one modeled result for", context)
    )
    if (nrow(observed) != 1L) {
      return(invisible(NULL))
    }

    effect <- observed$effect[[1]]
    valid_effect <- is.numeric(effect) && length(effect) == 1L && is.finite(effect)
    testthat::expect_true(
      valid_effect,
      info = paste("Expected one finite numeric effect for", context)
    )
    if (!valid_effect) {
      return(invisible(NULL))
    }
    testthat::expect_true(
      abs(effect) > sqrt(.Machine$double.eps),
      info = paste("Expected a non-zero effect for", context)
    )
    testthat::expect_equal(
      sign(effect), unname(direction),
      info = paste("Unexpected numerator-minus-denominator direction for", context)
    )
  })

  invisible(result)
}
