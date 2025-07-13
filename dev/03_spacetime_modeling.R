# INITIAL SETTINGS --------------------------------------------------------

# clean environment
rm(list = ls())

# install packages with respective version
renv::restore()

# load packages
list.packages <-  c("sf", "tidyverse", "terra", "ggplot2", "stars", "ggmosaic", "mgcv", "pROC")
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

# function to generate threshold data
threshold_predictor <- function(doy_value, aspect_value, treecover_value) {
  dplyr::tibble(
    P28 = rep(seq(from = 0, to = 199), times = 90),
    T0 = rep(rev(seq(from = -19, to = 40)), each = 300),
    doy = doy_value,
    aspect = aspect_value,
    treecoverdensity = treecover_value,
    landcover = factor(3)  # Directly set as factor
  )
}

# plot threshold combinations
plot_threshold <- function(data, topt) {
  subtitle_text <- paste0(
    "Day of the year: ", data$doy,
    "\nMonth: ", month.name[as.integer(format(as.Date(data$doy, origin = "2000-01-01"), "%m"))],
    "\nAspect: ", unique(data$aspect),
    "\nTree density cover: ", unique(data$treecoverdensity)
  )
  ggplot(data = data) +  
    geom_raster(aes(x = P28, y = T0, fill = prob)) + 
    geom_contour(aes(x = P28, y = T0, z = prob, colour = after_stat(level)), breaks = topt, color = "red", lty = 1, lwd = 1) +
    scale_fill_gradientn(colors = hcl.colors(5, "YlOrRd", rev = TRUE), limits = c(0, 1), breaks = round(seq(0, 1, length.out = 5), 2)) +
    scale_x_continuous(expand = c(0, 0), limits = c(-0.5, 200), breaks = seq(0, 200, by = 25)) +
    scale_y_continuous(expand = c(0, -1), limits = c(-19, 40), breaks = seq(-15, 40, by = 5)) +
    theme_plot() +
    theme(aspect.ratio = 1) + ggtitle("", subtitle = subtitle_text) +
    xlab("P: 29-day cumulative precipitation (mm)") +
    ylab("T: 1-day mean temperature (°)")
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

# space-time filtered sampling
d <- readRDS("./dat/processed/spacetime_sampling_filtered.Rds") %>% dplyr::mutate(bin = forcats::fct_relevel(bin, "0", "1")) %>% tidyr::drop_na(P0, landcover)


# MODELED RELATIONSHIPS ---------------------------------------------------

# values above do not relate to wildfire occurrence
d <- d %>% 
  dplyr::mutate(distbuildings = dplyr::if_else(distbuildings > 4000, 4000, distbuildings)) %>% 
  dplyr::mutate(distroads = dplyr::if_else(distroads > 1000, 1000, distroads)) %>% 
  dplyr::mutate(year = as.factor(year)) %>%
  dplyr::mutate(P30 = dplyr::if_else(P30 > 300, 300, P30)) %>% 
  dplyr::mutate(P28 = dplyr::if_else(P28 > 300, 300, P28)) 

# formula using the selected time windows
formula <- bin ~
  s(T0) +
  s(P28) +
  s(doy, bs = "cc", k = 3) +
  s(aspect, bs = "cc") +
  s(treecoverdensity) +
  s(distbuildings) +
  landcover +
  s(mean_temperature, k = 3) +
  s(mean_precipitation, k = 3) +
  s(year, bs = "re") 

# fit
fit <- mgcv::gam(formula, data=d, family=binomial, select = TRUE, method = "REML");summary(fit)

# plot partial effects
# pdf("./plt/03_partial_effects.pdf", width = 11, height = 8, paper ="a4r")
par(mfrow = c(3,4));par(pty="s")
plot(fit, pages=1, all.terms=T, trans=plogis, ylim = c(0,1), shade = T, shade.col = "grey", cex.lab = 1.5, cex.axis = 1.5)
# dev.off()

# predictions and fitting performance
my_exclude <- c("s(year)")
d$prob <- predict(fit, d, "response", exclude = my_exclude); summary(d$prob)
roc_bin <- pROC::roc(d$bin, d$prob);roc_bin


# VARIABLE IMPORTANCE -----------------------------------------------------

# updated formula
formula <- formula(drop.terms(terms(formula), drop = 10, keep.response = TRUE))
                   
# fit full model
mod <- mgcv::gam(formula, data=d, family=binomial, select = TRUE, method = "REML")
summary(mod)

# iterative eval
formula <-mod$formula
formula.terms <- labels(terms(formula))
formula.vars <- all.vars(formula)
prop.dev <- rep(NA, length(formula.terms))

# null model
mod.null <- mgcv::gam(formula(paste0(formula.vars[1], '~1'), select=TRUE, method="REML"), family=binomial, data=d)
summary(mod.null)

# permutation variable importance
n.permut <- 10
prop.dev.perm <- list()

for(i in 1:length(formula.terms)){
  term <- formula.terms[i]
  dev.increase <- numeric(n.permut)
  
  for (j in 1:n.permut){
    d.perm <- d  
    term.name <- gsub("s\\(|\\)", "", term)          
    term.name <- strsplit(term.name, ",")[[1]][1]    
    term.name <- trimws(term.name) 
    
    d.perm[[term.name]] <- sample(d.perm[[term.name]])
    
    #fit model with permuted variable
    mod.perm <- mgcv::gam(formula, data = d.perm, family = binomial, select =TRUE, method = "REML")
    
    # deviance explained
    dev.increase[j] <- ((deviance(mod.perm)-deviance(mod))/deviance(mod.null))
  }
  prop.dev.perm[[i]] <- dplyr::tibble(
    term = rep(term, n.permut),
    permutation = 1:n.permut,
    prop.dev = dev.increase
  )
}

# arrange results
perm_results <- do.call(rbind, prop.dev.perm) %>% 
  group_by(term) %>% 
  summarise(
    mean_dev = mean(prop.dev, na.rm =TRUE),
    sd_dev = sd(prop.dev, na.rm = TRUE))

# clean up
rm(dev.increase, formula.terms, formula.vars, d.perm, mod.null, mod.perm,
   i, j, n.permut, term, term.name, prop.dev.perm, prop.dev, mod)

# plot of permutation var importance
ggplot(perm_results, aes(x = reorder(term, mean_dev), y =mean_dev)) +
  geom_point()+
  geom_errorbar(aes(ymin=mean_dev-sd_dev, ymax=mean_dev+sd_dev), width = 0.2)+
  labs(x = "", y = "Permutation variable importance (deviance explained)") +
  coord_flip() +
  theme_plot()
# ggsave("./plt/04_variable_importance.pdf", dpi = 300, height = 20, width = 30, units = "cm")


# THRESHOLDING ------------------------------------------------------------

# OPT threshold
topt <- pROC::coords(roc_bin, x = "best", input="threshold",  best.method = "youden") %>% 
  dplyr::mutate(type = "optimal");topt

# TPR95 threshold
tsen <- subset(pROC::coords(roc_bin), sensitivity > 0.95) %>%  # 95% TPR and max specificity (or min FPR)
  dplyr::filter(sensitivity == min(sensitivity)) %>% 
  dplyr::filter(specificity == max(specificity)) %>% 
  dplyr::mutate(type = "sen"); tsen

# FPR5 threshold
tspe <- subset(pROC::coords(roc_bin), specificity > 0.95) %>%  # 95% TNR and max sensitivity (or max TPR)
  dplyr::filter(specificity == min(specificity)) %>% 
  dplyr::filter(sensitivity == max(sensitivity)) %>%
  dplyr::mutate(type = "spe");tspe

# object for storing thresholds
t <- dplyr::bind_rows(topt, tsen, tspe);t
rm(topt, tsen, tspe)

# exclude terms in prediction
my_exclude = c("s(aspect)", "s(treecoverdensity)",
               "s(distbuildings)", "s(mean_precipitation)",
               "s(mean_temperature)", "landcover",
               # "s(doy)",
               "s(year)")

# set a hypothetical doy (doy = 60) for illustration
# doy can be excluded as well
d_threshold <- dplyr::tibble(P28 = rep(seq(from = 0, to = 199), times = 90),
                             T0 = rep(rev(seq(from = -19, to = 40)), each = 300),
                             doy = 60,
                             landcover = 3) %>% 
  dplyr::mutate_all(as.numeric) %>% 
  dplyr::mutate_at(vars(landcover), as.factor)

# predictions
d_threshold$prob <- predict(fit, type = "response", newdata.guaranteed=TRUE, newdata = d_threshold, exclude = my_exclude);summary(d_threshold$prob)
d_bin1 <- dplyr::filter(d, bin ==1)
d_bin0 <- dplyr::filter(d, bin == 0)

# static threshold plot
ggplot(data=d_threshold) +  
  geom_raster (aes(x=P28, y = T0, fill=prob)) + 
  geom_contour(aes(x=P28, y = T0, z = prob, colour = after_stat(level)), breaks = t$threshold[1],  color="#C62527", lty=1, lwd=1) +
  geom_contour(aes(x=P28, y = T0, z = prob, colour = after_stat(level)), breaks = t$threshold[2],  color="#0D9B58", lty=2, lwd=1) +
  geom_contour(aes(x=P28, y = T0, z = prob, colour = after_stat(level)), breaks = t$threshold[3],  color="#5A5F94", lty=2, lwd=1) +
  scale_fill_gradientn(colors = hcl.colors(5, "YlOrRd", rev = 1), limits = c(0, 1), breaks=round(seq(0, 1, length.out = 5),2))+
  scale_x_continuous(expand = c(0, 0), limits = c(-0.5, 200), breaks = seq(0, 200, by = 25)) + 
  scale_y_continuous(expand = c(0, -1), limits = c(-19, 40), breaks = seq(-15, 40, by = 5)) +
  geom_point(data=d_bin1, aes(x=P30 , y=T1),shape=4, alpha=0.2,  size=2, col="black")+
  geom_point(data=d_bin0, aes(x=P30 , y=T1),shape=16,alpha=0.2, size=0.5, col="black") +
  theme_plot()+
  theme(aspect.ratio = 1)+
  ggtitle("", subtitle= paste0("Day of the year: ", unique(d_threshold$doy), "\nMonth: ",
                               month.name[as.integer(format(as.Date(unique(d_threshold$doy), origin = "2000-01-01"), format = "%m"))]))+
  xlab("P: 29-day cumulative precipitation (mm)") +
  ylab("T: 1-day mean temperature (°)") +
  annotate(geom="text", x=75, y=21.5, angle = 25, label=paste0("Optimal (TNR", round(t[1,2],2)*100, " TPR", round(t[1,3],2)*100,")"), col="#C62527") +
  annotate(geom="text", x=75, y=10.5, angle = 25, label=paste0("TNR", round(t[2,2],2)*100, " TPR", round(t[2,3],2)*100),   col="#0D9B58") +
  annotate(geom="text", x=75, y=32.5, angle = 25, label=paste0("TNR", round(t[3,2],2)*100, " TPR", round(t[3,3],2)*100),   col="#5A5F94")
# ggsave("./plt/05_threshold.pdf", dpi = 300, height = 15, width = 15, units = "cm")

# clean up
rm(d_bin0, d_bin1, d_threshold)

## examples thresholds -----------------------------------------------------

# wildfire probabilities under different tree cover densities, aspect and doy can be compared
my_exclude = c("s(distbuildings)", "s(mean_precipitation)",
  "s(mean_temperature)", "landcover",
  "s(year)")

# different combinations of factors
threshold1 <- threshold_predictor(doy_value = 30, aspect_value = 0, treecover_value = 30)
threshold2 <- threshold_predictor(doy_value = 240, aspect_value = 0, treecover_value = 30)
threshold3 <- threshold_predictor(doy_value = 30, aspect_value = 225, treecover_value = 100)
threshold4 <- threshold_predictor(doy_value = 240, aspect_value = 225, treecover_value = 100)

# predictions
threshold1$prob <- predict(fit, type="response", newdata.guaranteed=TRUE, newdata=threshold1, exclude = my_exclude)
threshold2$prob <- predict(fit, type="response", newdata.guaranteed=TRUE, newdata=threshold2, exclude = my_exclude)
threshold3$prob <- predict(fit, type="response", newdata.guaranteed=TRUE, newdata=threshold3, exclude = my_exclude)
threshold4$prob <- predict(fit, type="response", newdata.guaranteed=TRUE, newdata=threshold4, exclude = my_exclude)

# plots
threshold1 <- plot_threshold(threshold1, t$threshold[1])
threshold2 <- plot_threshold(threshold2, t$threshold[1])
threshold3 <- plot_threshold(threshold3, t$threshold[1])
threshold4 <- plot_threshold(threshold4, t$threshold[1])

# combined plots
pm <- gridExtra::grid.arrange(threshold1, threshold2, threshold3, threshold4, ncol=2, nrow = 2)
# ggsave("./plt/06_threshold_examples.pdf", pm, dpi = 500, height = 20, width = 20, units = "cm")

