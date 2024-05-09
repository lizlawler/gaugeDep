library(readr)
library(tidyverse)

extract_metadata <- function(one_chain_csv_file) {
  meta <- read_csv(one_chain_csv_file)$metadata
  date <- meta[which(grepl("date", meta, ignore.case = TRUE))]
  inits <- meta[which(grepl("initial", meta, ignore.case = TRUE))]
  settings <- meta[which(grepl("settings", meta, ignore.case = TRUE))]
  elapsed_time <- meta[which(grepl("elapsed", meta, ignore.case = TRUE))]
  accept_rate <- meta[which(grepl("accept", meta, ignore.case = TRUE))]
  start_samples <- which(grepl("\\<iter\\>", meta))
  end_samples <- which(grepl("\\<metadata\\>", meta, ignore.case = TRUE)) - 1
  mcmc_samples <- meta[(start_samples+1):end_samples] |> as.data.frame()
  colnames(mcmc_samples) <- meta[start_samples]
  split_cols <- str_split(meta[start_samples], ",") |> unlist()
  mcmc_samples <- mcmc_samples |>
    as_tibble() |>
    separate_wider_delim(cols = 1, delim = ",", names = split_cols) |> 
    mutate(across(where(is.character), as.numeric))
  return(list(date = date,
              inits = inits,
              settings = settings,
              elapsed_time = elapsed_time,
              accept_rate = accept_rate,
              mcmc_samples = mcmc_samples))
}
