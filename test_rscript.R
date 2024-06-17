args <- commandArgs(trailingOnly = TRUE)
char1 <- args[1]
dig1 <- args[2]

library(cmdstanr)
library(posterior)
library(dplyr)
library(tidyr)

1+1
1 + as.numeric(dig1)
char1

getwd()
