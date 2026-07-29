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
