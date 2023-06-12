#!/usr/bin/env Rscript
library(getopt)
library(tidyverse)

option_matrix <- matrix(
  c(
    "infile", "i", 1, "character",
    "outfile", "o", 1, "character"
  ),
  byrow = TRUE, ncol = 4
)

options <- getopt(option_matrix)



read_gtcheck <- function(gtcheck_filename) {
  read_tsv(
    file = gtcheck_filename,
    skip = 20,
    col_names = c(
      "dc", "query_sample", "genotyped_sample", "discordance", "log_p_hwe",
      "number_of_sites"
    )
  ) %>%
  select(-dc)
}

create_distance_matrix <- function(gtcheck_long) {
  gtcheck_diagonal <-
    gtcheck_long %>%
    select(query_sample, genotyped_sample) %>%
    mutate(discordance = 0) %>%
    pivot_longer(
      query_sample:genotyped_sample,
      names_to = "type",
      values_to = "query_sample"
    ) %>%
    select(query_sample) %>%
    mutate(
      genotyped_sample = query_sample,
      discordance = 0.0
    ) %>%
    distinct()

  gtcheck_upper <-
    gtcheck_long %>%
    select(query_sample, genotyped_sample, discordance)

  gtcheck_lower <-
    gtcheck_long %>%
    select(
      genotyped_sample = query_sample,
      query_sample = genotyped_sample,
      discordance
    )

  gtcheck_discordances <-
    bind_rows(gtcheck_diagonal, gtcheck_upper, gtcheck_lower) %>%
    arrange(query_sample)

  gtcheck_discordances_matrix <-
    gtcheck_discordances %>%
    pivot_wider(
      names_from = genotyped_sample,
      values_from = discordance
    ) %>%
    column_to_rownames("query_sample") %>%
    as.matrix()

  return(gtcheck_discordances_matrix)
}


if (!interactive()) {
  gtcheck_long <- read_gtcheck(options$infile)
  gtcheck_discordances_matrix <- create_distance_matrix(gtcheck_long)
  pdf(file = options$outfile, paper = "a4")
  heatmap(gtcheck_discordances_matrix)
  dev.off()
}
