# INITIAL SETTINGS --------------------------------------------------------

# clean environment
rm(list = ls())

# install packages with respective version
renv::restore()

# call packages
list.packages <-  c("sf", "tidyverse", "terra", "ggplot2", "stars", "ggmosaic")
vapply(list.packages, library, logical(1), character.only = TRUE, logical.return = TRUE, quietly = TRUE)
rm(list.packages)
