# INITIAL SETTINGS --------------------------------------------------------

# clean environment
rm(list = ls())

# install packages with respective version
renv::restore()

# load packages
list.packages <-  c("sf", "tidyverse", "terra", "stars", "ggmosaic", "mgcv", "pROC", "parallel", "sperrorest")
vapply(list.packages, library, logical(1), character.only = TRUE, logical.return = TRUE, quietly = TRUE)
rm(list.packages)

# custom plot theme
theme_plot <-  function() {
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        axis.text.y   = element_text(size=14),
        axis.text.x   = element_text(size=14),
        axis.title.y  = element_text(size=14, vjust= 2),
        axis.title.x  = element_text(size=14, vjust=-2),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5))
}

# function to generate object to iterate in fit
generate_myoptimal <- function(days) {
  preparatory <- names(d %>% dplyr::select(dplyr::matches("P\\d{1,2}")) %>% sf::st_drop_geometry())
  myoptimal <- as.data.frame(cbind(paste0("T", days), preparatory))
  names(myoptimal) <- c("temperature", "precipitation")
  return(myoptimal)
}

# function to accumulate precipitation backwards in time
accum_precip <- function(x) {
  for (i in 2:ncol(d_precipitation)){
    if (i == 2){
      x[i] = x[i]
    } else {
      x[i] <-  x[i] + x[i-1]
    }
  }
  return(x)
}

# function to plot monthly and yearly distribution of wildfires
month_freq <- function(d) {
  ggplot(d, aes(x=as.factor(month), y = after_stat(count), fill = forcats::fct_relevel(bin, "1", "0"))) + geom_bar(width = 0.7) +
    coord_cartesian(ylim= c(0, 600)) + ylab("") + xlab("Month") +
    stat_count(geom = "text", colour = "black", size = 4, aes(label = after_stat(count)),position=position_stack(vjust=0.5))+
    ggtitle("Wildfire presences-absences (Month)") + # less samples in 2012 since it starts from April
    scale_fill_discrete(name="", labels=c(paste("Presences:", paste(sum(d$bin==1))), paste("Absences:", paste(sum(d$bin==0))))) +
    theme_plot() + theme(legend.box.background = element_rect(color="black", linewidth=0.3),
                         legend.box.margin = margin(1, 1, 1, 1),
                         legend.text = element_text(size=11)) 
}

year_freq <- function(d) {
  ggplot(d, aes(x=as.factor(year), y = after_stat(count), fill = forcats::fct_relevel(bin, "1", "0"))) + geom_bar(width = 0.7) +
    coord_cartesian(ylim= c(0, 350)) + ylab("") + xlab("Year") +
    stat_count(geom = "text", colour = "black", size = 4, aes(label = after_stat(count)),position=position_stack(vjust=0.5))+
    ggtitle("Wildfire presences-absences (MYear)") + # less samples in 2012 since it starts from April
    scale_fill_discrete(name="", labels=c(paste("Presences:", paste(sum(d$bin==1))), paste("Absences:", paste(sum(d$bin==0))))) +
    theme_plot() + theme(legend.box.background = element_rect(color="black", linewidth=0.3),
                         legend.box.margin = margin(1, 1, 1, 1),
                         legend.text = element_text(size=11)) 
}


# LOAD INPUT DATA ---------------------------------------------------------

# aoi
aoi <- sf::st_read("./dat/raw/aoi.gpkg")

# covariates
grids <- list.files("./dat/raw/static", pattern = "50.tif$", full.names = T)
name_grids <- gsub("./dat/raw/static/", "", gsub("50.tif", "", grids))
name_grids[c(7, 10)] <- c("landcover", "elevation")
grids <- terra::rast(grids); names(grids) <- name_grids
rm(name_grids)

# space-time sampling
d <- readRDS("./dat/interim/spacetime_sampling.Rds") %>% dplyr::mutate(bin = forcats::fct_relevel(bin, "0", "1")) %>% tidyr::drop_na(P0, landcover)


# ACCUMULATE PRECIPITATION ----

# create additional predictors for modeling
d <-  dplyr::mutate(d, year = format(date_corr, format="%Y")) %>%
  dplyr::mutate(month = format(date_corr, format="%m")) %>%
  dplyr::mutate(day = format(date_corr, format="%d")) %>%
  dplyr::mutate(doy = lubridate::yday(date_corr)) %>%
  dplyr::mutate(across(year:doy, as.integer)) %>%
  dplyr::mutate(month_name = month.name[month]) %>% 
  dplyr::mutate(dow = weekdays.Date(date_corr)) %>% 
  dplyr::mutate(wday = dplyr::if_else(dow == "Saturday" | dow == "Sunday", "Weekend", "Working")) %>% 
  dplyr::relocate(c(month_name, dow, wday), .after = doy) %>% 
  dplyr::mutate_at(vars(year, month_name, dow, wday, hour, initiation_place,
                        vegetation_type, orography, ignition, cause, type), as.factor) %>% 
  dplyr::mutate(provinces = dplyr::if_else(provinces == 0, NA, provinces)) %>% 
  tidyr::drop_na(provinces) %>% 
  dplyr::mutate(provinces = droplevels(provinces)) 

# separating temperature and precipitation to produce cumulative values
d_precipitation <- dplyr::select(d, c(id_obs, dplyr::matches("P\\d{1,2}"))) %>% sf::st_drop_geometry()
d_rest <- dplyr::select(d, -c(dplyr::matches("P\\d{1,2}")))

# precipitation backwards in time
d_precipitation <- accum_precip(d_precipitation)

# combine outputs
d <- dplyr::left_join(d_rest, d_precipitation, by = "id_obs")

# clean up
rm(d_rest, d_precipitation)


# TRIVIAL TIMES -----------------------------------------------------------

# exclude rainy days during the wildfire occurrence
# 1.1 mm at the day of occurrence P0
d_filtered <- dplyr::filter(d, P0 < 1.1) # using shifted dates, it covers from 08:00 of the day after occurrence to 08:00 of the occurrence day
# saveRDS(d_filtered, "./dat/processed/spacetime_sampling_filtered.Rds")

# comparison plot
plot(gridExtra::arrangeGrob(
  year_freq(d),
  year_freq(d_filtered),
  ncol=2, nrow=1))

# keep filtered data
d <- d_filtered
rm(d_filtered)


# TEMPORAL CV MODEL -------------------------------------------------------

# time window selection dataframe
myoptimal <- do.call(rbind, lapply(0:10, generate_myoptimal))

# cross-validation with 10 folds and 10 repetitions 
resamp <- sperrorest::partition_cv(d,  nfold = 10, repetition = 10, seed1= 1) 

# object to iterate
grid <- data.frame(i = rep(1:341, each = 100),
                   j = rep(rep(1:10, each= 10), times=341),
                   k = rep(1:10, times = 3410))

# optimal time window based on predictive performance


# parallel settings
n_cores <- max(1, parallel::detectCores() - 6)
clust = parallel::makeCluster(n_cores)
result_auc = list()
parallel::clusterExport(cl = clust, varlist = c("myoptimal", "d", "resamp", "grid"))

# force parallel to use only 1 thread per process
Sys.setenv(OMP_NUM_THREADS = "1")
Sys.setenv(MKL_NUM_THREADS = "1")
Sys.setenv(NUMEXPR_NUM_THREADS = "1")
Sys.setenv(OPENBLAS_NUM_THREADS = "1")

# parallel loop
set.seed(1) 
a <- Sys.time(); a
result_auc <- parallel::parLapply(cl = clust, seq_len(nrow(grid)), function(idx){
  row <- grid[idx, ]
  i <- row$i
  j <- row$j
  k <- row$k
  formula <- as.formula(paste0("bin ~ ",
                               "s(", myoptimal[i,1], ", k=5) + ",
                               "s(", myoptimal[i,2], ", k=5) + ",
                               "s(doy, bs='cc', k=5) + ",
                               "s(year, bs='re')"))
  test_ids <- resamp[[j]][[k]]$test
  d_train <- d[-test_ids, ]
  d_test <- d[test_ids, ]
  fit <- mgcv::bam(formula, data=d_train, family=binomial, method="fREML", discrete = 100)
  d_test$fit <- predict(fit, newdata=d_test, type="response")
  auc <- as.numeric(pROC::auc(d_test$bin, d_test$fit))
  return(list(combination = i, repetition = j, fold = k, auc = auc))
})
b <- Sys.time()-a; b #2.7 min/30 cores/125G RAM
parallel::stopCluster(cl = clust)
rm(a, b, n_cores, clust, grid)

# arrange stored aucs
my_auc <- do.call(rbind, result_auc) %>% 
  tibble::as_tibble() %>% 
  dplyr::mutate_all(as.numeric)%>%
  dplyr::group_by(combination) %>% 
  summarise(
    auc_mean = mean(auc, na.rm =T),
    auc_median = median(auc, na.rm = TRUE),
    auc_sd = sd(auc, na.rm = TRUE),
    auc_iqr = IQR(auc, na.rm = TRUE)) %>%
  dplyr::select(-combination) %>% 
  dplyr::bind_cols(myoptimal, .) %>%
  dplyr::mutate(myoptimal, precipitation = gsub("P", "", precipitation)) %>%
  dplyr::mutate(temperature = gsub("T", "", temperature)) %>% 
  dplyr::mutate_all(as.numeric)

# rename an clean up
myoptimal <- my_auc
rm(my_auc)

# plot
title <- myoptimal %>% dplyr::filter(auc_mean == max(auc_mean))
ggplot(myoptimal, aes(y = temperature, x = precipitation, fill = auc_mean))+
  geom_tile(color = "white", lwd = 0.05, linetype = 1)+
  scale_fill_gradientn(colors = hcl.colors(50, "Reds", rev = 1))+
  ggtitle(paste0("AUC = ", round(title[3], 3), ", precipitation = ", title[2], ", temperature = ", title[1])) + 
  scale_y_continuous(expand = c(0, 0), breaks = 0:10) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 180, by = 10)) +
  coord_equal() +
  xlab("Day-window cumulative precipitation")+
  ylab("Day-window mean temperature") +
  guides(fill = guide_colourbar(barwidth = 30, barheight = 1, title = "AUROC")) +
  theme_plot() + theme(legend.position = "bottom")
# ggsave("./FIGURES/optwinCV.pdf", dpi = 500, height = 20, width = 30, units = "cm")





