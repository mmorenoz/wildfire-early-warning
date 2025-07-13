# INITIAL SETTINGS --------------------------------------------------------

# clean environment
rm(list = ls())

# install packages with respective version
renv::restore()

# load packages
list.packages <-  c("sf", "tidyverse", "sperrorest", "mgcv")
vapply(list.packages, library, logical(1), character.only = TRUE, logical.return = TRUE, quietly = TRUE)
rm(list.packages)

# custom plot theme
theme_plot <- function(){
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        legend.box.background = element_rect(color="black", linewidth=0.3),
        legend.box.margin = margin(1, 1, 1, 1),
        panel.grid.minor = element_blank(),
        axis.text.y   = element_text(size=14),
        axis.text.x   = element_text(size=14),
        axis.title.y  = element_text(size=14, vjust= 2),
        axis.title.x  = element_text(size=14, vjust=-2),
        legend.position = c(0.11,0.9),
        plot.margin = unit(c(1,2,1,1), "lines"),
        legend.text = element_text(size=11),
        plot.title = element_text(hjust = 0.5))
}

# load fit
load("./dat/processed/fit.Rda")


# VALIDATION --------------------------------------------------------------

# separate X and Y coordinates for plot
d <- d %>% 
  sf::st_centroid() %>% 
  dplyr::mutate(x = sf::st_coordinates(.)[,1]) %>% 
  dplyr::mutate(y = sf::st_coordinates(.)[,2]) %>% 
  st_drop_geometry()


## random cv ---------------------------------------------------------------

# setting up loop for random cross-validation
fold <- 10
repetition <- 10

# create random cv partition
set.seed(1)
partition <- sperrorest::partition_cv(d, nfold = fold, repetition = repetition, seed1 = 123) 
myroc_cv <- lapply(my.list<-vector(mode = 'list', 10), function(x) x<-vector(mode = 'list', 10))

# loop for validation
for (i in 1:repetition){
  id.hold <- partition[[i]]
  for (j in 1:fold){
    id.holdout <- id.hold[[j]]$test
    d.test <- d[id.holdout, ]
    d.train <- d[-id.holdout, ]
    fit <- mgcv::gam(formula, data = d.train, family = binomial, method = "REML")
    d.test$prediction <- predict(fit, d.test, type = "response", exclude = my_exclude)
    myroc <- unlist(pROC::roc(response = d.test$bin, predictor = d.test$prediction, auc=T))$auc
    myroc_cv[[i]][[j]] <- myroc
  }
}

# clean environment
rm(d.train, d.test, id.hold, id.holdout, myroc, i, j)

# store performances in data frame
myroc_cv = dplyr::as_tibble(do.call(rbind, lapply(myroc_cv, unlist))) %>% 
  tidyr::gather(key = "repetition", value = "auc") %>% 
  dplyr::mutate(repetition = rep(1:10, each = 10))


## spatial cv ---------------------------------------------------------------

# setting up loop for spatial cross-validation
fold <- 10
repetition <- 10

# create spatial scv partition
set.seed(1)
partition.s <- sperrorest::partition_kmeans(d, nfold = fold, repetition = repetition, seed1 = 123) 
myroc_scv <- lapply(my.list<-vector(mode = 'list',10), function(x) x<-vector(mode='list',10))

# loop for validation
for (i in 1:repetition){
  id.hold <- partition.s[[i]]
  for (j in 1:fold){
    id.holdout <- id.hold[[j]]$test
    d.test <- d[id.holdout, ]
    d.train <- d[-id.holdout, ]
    fit <- mgcv::gam(formula, data = d.train, family = binomial, method = "REML")
    d.test$prediction <- predict(fit, d.test, type = "response", exclude = my_exclude)
    myroc <- unlist(pROC::roc(response = d.test$bin, predictor = d.test$prediction, auc = T))$auc
    myroc_scv[[i]][[j]] <- myroc
  }
}

# clean environment
rm(d.train, d.test, id.hold, id.holdout, myroc, i, j)

# store performances in data frame
myroc_scv <- dplyr::as_tibble(do.call(rbind, lapply(myroc_scv, unlist))) %>%   
  tidyr::gather(key = "repetition", value = "auc") %>% 
  dplyr::mutate(repetition = rep(1:10, each = 10))


## factor cv ---------------------------------------------------------------

### land cover --------------------------------------------------------------

summary(d$landcover)

# create partitions
set.seed(1)
partition.land <- sperrorest::partition_factor(d, fac = "landcover") 

# setting up loop for cross-validation
fold <- length(partition.land[[1]])
myroc_lcv <- c()

# loop for validation
for (i in 1:fold){
  id.holdout <- partition.land[[1]][[i]]$test 
  d.test <- d[id.holdout, ]
  d.train <- d[-id.holdout, ]
  fit <- mgcv::gam(formula, data = d.train, family = binomial, method = "REML", drop.unused.levels = F)
  d.test$prediction <- predict(fit, d.test, type = "response", exclude = my_exclude)
  myroc <- unlist(pROC::roc(response = d.test$bin, predictor = d.test$prediction, auc = T))$auc
  myroc_lcv[i] <- myroc
}

# clean environment
rm(d.train, d.test, id.holdout, myroc, i)

# store performances in tibble
myroc_lcv <- dplyr::tibble(fold = c(1:5),
                           landcover = c(
                             "Croplands",
                             "Decidious forest",
                             "Coniferous forest",
                             "Herbaceous",
                             "Marshes"), 
                           auc = myroc_lcv) %>% 
  dplyr::mutate(auc = as.numeric(auc)) %>% 
  dplyr::mutate(fold = as.integer(fold)) %>%  
  dplyr::arrange(fold)


### month --------------------------------------------------------------

# convert to factor
d <- dplyr::mutate(d, month = as.factor(month))
summary(d$month)

# create partition
set.seed(1)
partition.month <- sperrorest::partition_factor(d, fac = "month") 

# setting up loop for cross-validation
fold <- length(partition.month[[1]])
myroc_mcv <- c()

# loop for validation
for (i in 1:fold){
  id.holdout <- partition.month[[1]][[i]]$test
  d.test <- d[id.holdout, ]
  d.train <- d[-id.holdout, ]
  fit <- mgcv::gam(formula, data = d.train, family = binomial, method = "REML", drop.unused.levels = F)
  d.test$prediction <- predict(fit, d.test, type = "response", exclude = my_exclude)
  myroc <- unlist(pROC::roc(response = d.test$bin, predictor = d.test$prediction, auc = T))$auc
  myroc_mcv[i] <- myroc
}

# clean environment
rm(d.train, d.test, id.holdout, myroc, i)

# store performances in data frame
myroc_mcv <- dplyr::tibble(fold = seq(1,12,1),
                           month = seq(1,12,1), 
                           auc = myroc_mcv) %>% 
  dplyr::mutate(auc = as.numeric(auc)) %>% 
  dplyr::mutate(fold = as.integer(fold)) %>%  
  dplyr::arrange(month)


### year --------------------------------------------------------------

# convert to factor
d <- dplyr::mutate(d, year = as.factor(year))
summary(d$year)

# create partition
set.seed(1)
partition.year <- sperrorest::partition_factor(d, fac = "year") 

# setting up loop for cross-validation
fold <- length(partition.year[[1]])
myroc_ycv <- c()

# loop for validation
for (i in 1:fold){
  id.holdout <- partition.year[[1]][[i]]$test
  d.test <- d[id.holdout, ]
  d.train <- d[-id.holdout, ]
  fit <- mgcv::gam(formula, data = d.train, family = binomial, method = "REML", drop.unused.levels = F)
  d.test$prediction <- predict(fit, d.test, type = "response", exclude = my_exclude)
  myroc <- unlist(pROC::roc(response = d.test$bin, predictor = d.test$prediction, auc = T))$auc
  myroc_ycv[i] <- myroc
}

# clean environment
rm(d.train, d.test, id.holdout, myroc, i)

# store performances in data frame
myroc_ycv = dplyr::tibble(fold = seq(1,24,1),
                          year = c(2000, 2001, 2002, 2003, 2004, 2005, 2006,
                                   2007, 2008, 2009, 2010, 2011, 2012, 2013,
                                   2014, 2015, 2016, 2017, 2018, 2019, 2020,
                                   2021, 2022, 2023), 
                          auc = myroc_ycv) %>%  
  dplyr::mutate(auc = as.numeric(auc)) %>% 
  dplyr::mutate(fold = as.integer(fold)) %>%  
  dplyr::arrange(year)


# store environment data
# save.image("./dat/processed/cross_validation.RData")


# PLOTS -------------------------------------------------------------------

# loading data
load("./dat/processed/cross_validation.RData")

## RANDOM AND SPATIAL ----

# combine cv and scv
cv <- dplyr::select(myroc_scv, auc) %>% 
  dplyr::rename(auc_scv=auc) %>% 
  dplyr::bind_cols(myroc_cv) %>%
  dplyr::rename(auc_cv=auc) %>% 
  dplyr::relocate(repetition, auc_cv, .before = auc_scv) %>% 
  tidyr::pivot_longer(cols = !repetition, names_to = "type", values_to = "auc") 

# plot
ggplot2::ggplot(cv, aes(auc, type)) + geom_boxplot(width = 0.5) + 
  ggtitle(paste0(paste0("IQR_cv= "), paste0(round(IQR(myroc_cv$auc),3)),
                 paste0(" median AUC_cv= "), paste0(round(median(myroc_cv$auc), 2)),"\n",
                 paste0(paste0("IQR_scv= "), paste0(round(IQR(myroc_scv$auc),3)),
                        paste0(" median AUC_scv= "), paste0(round(median(myroc_scv$auc), 2)))))+
  ylab("") + xlab ("AUC")+
  coord_flip()+
  theme_plot()
# ggsave("./plt/07_rcv_scv.pdf", width = 7, height = 8, dpi = 300, units = "in", device='pdf')


## LAND COVER ----

# order land cover levels
myroc_lcv <- myroc_lcv %>% 
  dplyr::mutate(landcover = ordered(landcover, levels = c("Croplands", "Decidious forest", "Coniferous forest", "Herbaceous", "Marshes")))

# plot
ggplot(myroc_lcv, aes(x=landcover, y=auc)) + geom_bar(stat="identity", width=0.85) +
  xlab("Land cover") + ylab("AUC") + ggtitle(paste("mean AUC:", round(mean(myroc_lcv$auc),3))) +
  coord_cartesian(ylim=c(0.60,1)) +
  theme_plot()
# ggsave("./plt/08_cv_landcover.pdf", width = 7, height = 8, dpi = 300, units = "in", device='pdf')


## MONTH ----

# months to factor
myroc_mcv <- dplyr::mutate(myroc_mcv, month = as.factor(month))

# plot
ggplot(myroc_mcv, aes(x=as.factor(month), y=auc)) + geom_bar(stat="identity", width=0.85) +
  xlab("Month") + ylab("AUC") + ggtitle(paste("mean AUC:", round(mean(myroc_mcv$auc),3))) +
  coord_cartesian(ylim=c(0.60,1)) +
  theme_plot()
# ggsave("./plt/09_cv_month.pdf", width = 7, height = 8, dpi = 300, units = "in", device='pdf')


## YEAR ----

# years to factor
myroc_ycv <- dplyr::mutate(myroc_ycv, year = as.factor(year))

# plot
ggplot(myroc_ycv, aes(x=as.factor(year), y=auc)) + geom_bar(stat="identity", width=0.85) +
  xlab("Year") + ylab("AUC") + ggtitle(paste("mean AUC:", round(mean(myroc_ycv$auc),3))) +
  coord_cartesian(ylim=c(0.60,1)) +
  theme_plot()
# ggsave("./plt/10_cv_year.pdf", width = 7, height = 8, dpi = 300, units = "in", device='pdf')

