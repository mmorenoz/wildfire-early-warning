# INITIAL SETTINGS --------------------------------------------------------

# clean environment
rm(list = ls())

# install packages with respective version
renv::restore()

# load packages
list.packages <-  c("sf", "tidyverse", "terra", "ggplot2", "stars", "ggmosaic")
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


# LOAD INPUT DATA ---------------------------------------------------------

# aoi
aoi <- sf::st_read("./dat/raw/aoi.gpkg")

# covariates
grids <- list.files("./dat/raw/static", pattern = "50.tif$", full.names = T)
name_grids <- gsub("./dat/raw/static/", "", gsub("50.tif", "", grids))
name_grids[c(7, 10)] <- c("landcover", "elevation")
grids <- terra::rast(grids); names(grids) <- name_grids
rm(name_grids)

# inventory
d <- sf::st_read("./dat/interim/wildfire_point.gpkg") %>% 
  dplyr::select(-c(comune, location)) %>% 
  dplyr::mutate(bin = 1) %>%
  dplyr::filter(area > 10) %>% 
  dplyr::filter(date >= as.Date("2000/01/01") & date <= as.Date("2023/12/31")) %>% 
  sf::st_centroid(.) %>% 
  sf::st_buffer(., 50) %>% 
  dplyr::rename(id_original = id) %>% 
  dplyr::mutate(id = 1:nrow(.)) %>% 
  dplyr::mutate(date_corr = date + 1) %>% 
  dplyr::relocate(c(id, bin, area, date, year, month, day, doy, date_corr), .after = id_original)


# ABSENCES IN PRESENCE LOCATIONS ------------------------------------------

# the inventory has data from Jan 2000 to Dec 2023
# keep the analyses for the same temporal domain
summary(d$date)
difftime(as.Date("2023/12/31"), as.Date("2000/01/01"), units="days") + 1

# create 50 replicates at each wildfire location
# keep only one replicate with bin = 1 (wildfire) and the rest as bin = 0 (absences)
set.seed(7)
d <- dplyr::slice(d, rep(1:n(), each = 50)) %>%
  dplyr::mutate(id_obs = 1:nrow(.)) %>%
  dplyr::relocate(id_obs, .after = id) %>% 
  dplyr::mutate(date_corr = dplyr::if_else(duplicated(id), NA, date_corr)) %>% 
  dplyr::mutate(bin = dplyr::if_else(duplicated(id), 0, bin)) %>% 
  dplyr::mutate(date_corr = dplyr::if_else(is.na(date_corr), 
                                           sample(seq(as.Date('2000/01/01'),
                                                      as.Date('2023/12/31'),
                                                      by="day"), replace = T, nrow(.)), date_corr)) 

# separate presences and their (absence) replicates
d_bin1 <- dplyr::filter(d, bin == 1)
d_bin0 <- dplyr::filter(d, bin == 0) 

# filter/exclude the (absences) replicates (at the same location as presences) that:
# have a date after the wildfire occurrence: set up a 5 year time window  
# have a date before the wildfire occurrence: set up a 30 day (1 month) interval 
d_bin0_filtered <- d_bin0[0, ]
d_bin0_filtered <- do.call(rbind, lapply(d_bin1$id, function(x) {
  datebin1 <- dplyr::filter(d_bin1, id == x)$date_corr
  locid <- dplyr::filter(d_bin0, id == x) %>%
    dplyr::mutate(date_corr = dplyr::if_else(date_corr <= (datebin1 - 30) | date_corr >= datebin1 + 1825, date_corr, NA)) %>%
    tidyr::drop_na(date_corr)
  return(locid)
}))

# restrict the sampling in presence locations: fix 2 absences per presence locations
absence_loc <- function(x, loc, fact) {
  x <- x %>%
    dplyr::filter(id == loc) %>%
    dplyr::sample_n(size = round(fact), replace = FALSE)
}

fact <- 2
d_bin0_filtered <- do.call(dplyr::bind_rows, lapply(min(d_bin0_filtered$id):max(d_bin0_filtered$id), function(loc) {
  absence_loc(d_bin0_filtered, loc, fact)
}))


# ABSENCES IN ABSENCE LOCATIONS -------------------------------------------

# read land cover classes and replace by trivial terrain
landcoverclass <- read.csv("./dat/raw/static/landcover_newlegend.csv", sep = ",", header = T) 
trivial <- terra::subst(grids$landcover, landcoverclass$new, landcoverclass$trivial)

# burned area polygons for absences
burned <- sf::st_read("./dat/interim/wildfire_polygon.gpkg") %>%
  dplyr::filter(date >= as.Date("2000/01/01") & date <= as.Date("2023/12/31")) %>% 
  dplyr::mutate(bin = 1) %>% 
  dplyr::filter(area > 10) %>% 
  dplyr::select(bin) %>% 
  stars::st_rasterize(st_as_stars(st_bbox(grids), nx = 3353, ny = 3215)) %>% 
  terra::rast() %>% 
  terra::mask(aoi)

# invert mask
burned <- terra::ifel(burned == 1, 0, 1) 

# definition of sampling area
sampling <- trivial + burned
sampling <- terra::ifel(sampling == 2, 1, NA)
rm(landcoverclass, trivial, burned)

# sampling absences in space
set.seed(7)
n <- 3 #times the number of presences
absences <-  terra::spatSample(sampling, 3*nrow(d_bin1), as.points = T, method = "random", values = F, cells = F, na.rm = T, exhaustive = T) %>%
  sf::st_as_sf() %>%
  dplyr::mutate(bin = 0) %>% 
  dplyr::rename(geom = geometry) %>% 
  sf::st_buffer(., 50) 

# replicate absences (50)
absences <- dplyr::mutate(absences, id = (nrow(d_bin1)+1):((n+1)*nrow(d_bin1))) %>%  
  dplyr::mutate(date_corr = sample(seq(as.Date('2000/01/01'),
                                       as.Date('2023/12/31'),
                                       by = 'day'), replace = T, nrow(.)), date_corr) %>% 
  dplyr::slice(rep(1:n(), each = 50)) %>% 
  dplyr::mutate(id_obs = (nrow(d)+1):(nrow(d)+nrow(.))) %>% 
  dplyr::relocate(c(id, id_obs, bin, date_corr), .before = geom) %>% 
  dplyr::mutate(date_corr = dplyr::if_else(duplicated(id), NA, date_corr)) %>% 
  dplyr::mutate(date_corr = dplyr::if_else(is.na(date_corr), 
                                           sample(seq(as.Date('2000/01/01'),
                                                      as.Date('2023/12/31'),
                                                      by="day"), replace = T, nrow(.)), date_corr)) %>% 
  dplyr::mutate(date = date_corr - 1) %>% 
  dplyr::mutate(year = format(date, format="%Y")) %>% 
  dplyr::mutate(month = format(date, format="%m")) %>% 
  dplyr::mutate(day = format(date, format="%d")) %>% 
  dplyr::mutate(doy = lubridate::yday(date)) %>% 
  dplyr::mutate(across(year:doy, as.integer)) %>% 
  dplyr::relocate(c(date, year, month, day, doy), .before = date_corr)

# combine the absences at presence locations and absence locations
d_bin0_filtered <- dplyr::bind_rows(d_bin0_filtered, absences)
rm(absences, n)

# avoid overlapping times between observations (absences)
# assume 30 days 
locid <- unique(d_bin0_filtered$id)
d_bin0_filtered_days <- d_bin0[0, ]
a <- Sys.time(); a
for (x in locid){
  d_bin0_temp <- dplyr::filter(d_bin0_filtered, id == x) 
  datebin0 <- unique(d_bin0_temp$date_corr)
  for (y in 1:length(datebin0)){
    date_ite <- as.Date(d_bin0_temp$date_corr[y])
    if (!is.na(date_ite)){ 
      d_bin0_temp <- dplyr::mutate(d_bin0_temp, date_corr = dplyr::if_else(date_corr == date_ite, date_corr,
                                                                           dplyr::if_else(dplyr::between(date_corr, date_ite - 30, date_ite + 30),
                                                                                          NA, date_corr))) %>% 
        tidyr::drop_na(date_corr)
    }
  }
  d_bin0_filtered_days <- rbind(d_bin0_filtered_days, d_bin0_temp)
}
b <- Sys.time()-a; b #~30 min

# clean up
rm(x, y, date_ite, datebin0, locid, d_bin0_temp, a, b)


# BALANCING ACROSS TIME ---------------------------------------------------

## YEAR -----

# sample 2000 observations 
d_bin0 <- d_bin0_filtered_days
n <- 20000

# 20,000 absences are to be distributed within 24 years (2000 - 2023)
y365 <- 365/(365*24); y365*100 
24*round(y365*n) 

# yearly sampling
absence_yearly <- function(x, year_s, fact) {
  x <- x %>%
    dplyr::filter(year == year_s) %>%
    dplyr::sample_n(size = round(fact), replace =F)
}

# yearly balancing
d_bin0 <- do.call(dplyr::bind_rows, lapply(2000:2023, function(year) {
  absence_yearly(d_bin0, year, y365*n)
}))
table(d_bin0$year)


## MONTHS -----

# attempt a ratio of 1:5
n <-  nrow(d_bin1)*5 

# 4,990 absences are to be distributed within 12 months
m31 <-  31/365;m31*100 
m30 <-  30/365; m30*100 
m28 <-  28/365; m28*100  

# final amount of absences
7*round(m31*n) + 4*round(m30*n) + round(m28*n)

# monthly sampling
absence_monthly <- function(x, month_s, fact) {
  x <-  x %>%
    dplyr::filter(month == month_s & bin ==0) %>%
    dplyr::sample_n(size = round(fact), replace = F)
}

# sample for each group of months based on days
months_28 <- 2
months_30 <- c(4, 6, 9, 11)
months_31 <- c(1, 3, 5, 7, 8, 10, 12)

# generate absence samples by month length
d_bin0_28 <-  absence_monthly(d_bin0, 2, m28*n)  
d_bin0_30<- purrr::map_dfr(months_30, ~ absence_monthly(d_bin0, .x, m30 * n))
d_bin0_31 <- purrr::map_dfr(months_31, ~ absence_monthly(d_bin0, .x, m31 * n))

# combine absences 
d_bin0 <-  dplyr::bind_rows(d_bin0_28, d_bin0_30, d_bin0_31)

# combine presences and absences
d <- dplyr::bind_rows(d_bin0, d_bin1) %>%
  dplyr::mutate(bin = as.factor(bin)) %>% 
  dplyr::mutate(bin = forcats::fct_relevel(bin, "1", "0"))

# clean up
rm(d_bin0, d_bin1, d_bin0_28, d_bin0_30, d_bin0_31, d_bin0_filtered,
   d_bin0_filtered_days, m28, m30, m31, months_28, months_30,
   months_31, n, y365)


# PLOTS -------------------------------------------------------------------

# year
ggplot(d, aes(x=as.factor(year), y = after_stat(count), fill = forcats::fct_relevel(bin, "1", "0"))) + geom_bar(width = 0.7) +
  coord_cartesian(ylim= c(0, 400)) + ylab("") + xlab("Year") +
  stat_count(geom = "text", colour = "black", size = 4, aes(label = after_stat(count)),position=position_stack(vjust=0.5))+
  ggtitle("Wildfire presences-absences (Year)") + 
  scale_fill_discrete(name="", labels=c(paste("Presences:", paste(sum(d$bin==1))), paste("Absences:", paste(sum(d$bin==0))))) +
  theme_plot() + theme(legend.box.background = element_rect(color="black", linewidth=0.3),
                       legend.box.margin = margin(1, 1, 1, 1),
                       legend.position = c(0.10,0.93),
                       legend.text = element_text(size=11)) 

# month
ggplot(d, aes(x=as.factor(month), y = after_stat(count), fill = forcats::fct_relevel(bin, "1", "0"))) + geom_bar(width = 0.7) +
  coord_cartesian(ylim= c(0, 750)) + ylab("") + xlab("Month") +
  stat_count(geom = "text", colour = "black", size = 4, aes(label = after_stat(count)),position=position_stack(vjust=0.5))+
  ggtitle("Wildfire presences-absences (Month)") + # less samples in 2012 since it starts from April
  scale_fill_discrete(name="", labels=c(paste("Presences:", paste(sum(d$bin==1))), paste("Absences:", paste(sum(d$bin==0))))) +
  theme_plot() + theme(legend.box.background = element_rect(color="black", linewidth=0.3),
                       legend.box.margin = margin(1, 1, 1, 1),
                       legend.position = c(0.10,0.93),
                       legend.text = element_text(size=11)) 


# TEMPERATURE EXTRACTION  -------------------------------------------------

# path for temp
path <- "./dat/raw/dynamic/temp_grids/2022/DTMEAN_202204.nc"
read_nc <- function(x) {
  variable <- terra::rast(x)
  names(variable) <- gsub("temperature_DATE=", "", names(variable)) %>%
    gsub("tmean_DATE=", "", .) %>% 
    gsub("\\..*", "", .) %>%
    as.numeric() %>% 
    as.Date(origin = "1970-01-01")
  return(variable)
}
# plot(read_nc(path)[[1]])

# list dates and object to iterate
day <-  0 # initial day
days_back <-  30 # days back
d_list <- dplyr::filter(d, !is.na(month)) %>% 
  split(., seq(nrow(.)))

# adding names to dataframe
name <- unlist(lapply(0:days_back, function(x) paste0("T", x)))
d_variable <- as.data.frame(matrix(nrow = 1, ncol = (days_back + 2)))
names(d_variable) <- c("id_obs", name)
rm(name)

# paralleled loop to extract precipitation
# define cores
cl <-  parallel::makeCluster(parallel::detectCores() - 5)
doParallel::registerDoParallel(cl)

# to store the processing time
library(foreach)
a <-  Sys.time()
d_temperature <-  foreach(i = 1:length(d_list), .combine="rbind", .packages = c("terra", "sf", "exactextractr")) %dopar% {
  date <-  as.Date(d_list[[i]]$date)
  date_list <-  seq(date-days_back, date, by = "day")
  date_list <-  rev(date_list)
  variable_stack <-  terra::rast()
  for (j in 1:length(date_list)){
    year <-  format(date_list[j], "%Y")
    month <- format(date_list[j], "%m")
    day <-  format(date_list[j], "%d")
    variable <- read_nc(paste0(paste0("./dat/raw/dynamic/temp_grids/"), 
                               paste0(year), paste0("/DTMEAN_"), paste0(year), paste0(month), paste0(".nc"))) %>% 
      terra::subset(grep(paste0(paste0(year), paste0("."), paste0(month), paste0("."), paste0(day)),
                         names(.), value = T))
    terra::add(variable_stack) <- variable 
  }
  # change i for 1 is using .combine="rbind" inside the foreach loop
  d_variable[1, ] <-  suppressWarnings(cbind(d_list[[i]]$id_obs,
                                             data.frame(exactextractr::exact_extract(variable_stack,
                                                                                     d_list[[i]],
                                                                                     fun = "mean"))))
  return(d_variable)
}
print(Sys.time()-a) 
# 4.6 min for 5989 obs, 31 cores, 125G RAM
parallel::stopCluster(cl)


# PRECIPITATION EXTRACTION ------------------------------------------------

# path for precipitation
path <- "../../PROSLIDE/DATA/prec_grids/2023/DAILYPCP_202308.nc"
read_nc <- function(x) {
  variable <- terra::rast(x)
  names(variable) <- gsub("precipitation_DATE=", "", names(variable)) %>%
    gsub("prec_DATE=", "", .) %>% 
    gsub("\\..*", "", .) %>%
    as.numeric() %>% 
    as.Date(origin = "1970-01-01")
  return(variable)
}

# list dates and object to iterate
day <-  0 # initial day
days_back <-  30 # days back
d_list <- dplyr::filter(d, !is.na(month)) %>% 
  split(., seq(nrow(.)))

# adding names to dataframe
name <- unlist(lapply(0:days_back, function(x) paste0("P", x)))
d_variable <- as.data.frame(matrix(nrow = 1, ncol = (days_back + 2)))
names(d_variable) <- c("id_obs", name)
rm(name)

# paralleled loop to extract precipitation
# define cores
cl <-  parallel::makeCluster(parallel::detectCores() - 5)
doParallel::registerDoParallel(cl)

# to store the processing time
library(foreach)
a <-  Sys.time()
d_precipitation <-  foreach(i = 1:length(d_list), .combine="rbind", .packages = c("terra", "sf", "exactextractr")) %dopar% {
  date <-  as.Date(d_list[[i]]$date)
  date_list <-  seq(date-days_back, date, by = "day")
  date_list <-  rev(date_list)
  variable_stack <-  terra::rast()
  for (j in 1:length(date_list)){
    year <-  format(date_list[j], "%Y")
    month <- format(date_list[j], "%m")
    day <-  format(date_list[j], "%d")
    variable <- read_nc(paste0(paste0("./dat/raw/dynamic/prec_grids/"), 
                               paste0(year), paste0("/DAILYPCP_"), paste0(year), paste0(month), paste0(".nc"))) %>% 
      terra::subset(grep(paste0(paste0(year), paste0("."), paste0(month), paste0("."), paste0(day)),
                         names(.), value = T))
    terra::add(variable_stack) <- variable 
  }
  # change i for 1 is using .combine="rbind" inside the foreach loop
  d_variable[1, ] <-  suppressWarnings(cbind(d_list[[i]]$id_obs,
                                             data.frame(exactextractr::exact_extract(variable_stack,
                                                                                     d_list[[i]],
                                                                                     fun = "mean"))))
  return(d_variable)
}
print(Sys.time()-a) 
# 4.6 min for 5989 obs, 31 cores, 125G RAM
parallel::stopCluster(cl)


# COMBINE DATA ----

# join the precipitation and temperature data
d_meteo <- dplyr::left_join(d_temperature, d_precipitation, by = "id_obs")

# join meteo data with original spacetime dataset
d <- dplyr::left_join(d, d_meteo, by = "id_obs")

# clean up
rm(d_list, a, day, days_back, path, d_precipitation, d_temperature, d_variable, d_meteo, cl)


# STATIC EXTRACTION -------------------------------------------------------

# extraction of predictors
d <- d %>% 
  dplyr::mutate(aspect = exactextractr::exact_extract(grids$aspect, d, "mean")) %>% 
  dplyr::mutate(northness = exactextractr::exact_extract(grids$northness, d, "mean")) %>%
  dplyr::mutate(eastness = exactextractr::exact_extract(grids$eastness, d, "mean")) %>% 
  dplyr::mutate(elevation = exactextractr::exact_extract(grids$elevation, d, "mean")) %>% 
  dplyr::mutate(slope = exactextractr::exact_extract(grids$slope, d, "mean")) %>% 
  dplyr::mutate(tri = exactextractr::exact_extract(grids$tri, d, "mean"))  %>% 
  dplyr::mutate(landcover = as.factor(exactextractr::exact_extract(grids$landcover, d, "majority"))) %>% 
  dplyr::mutate(treecoverdensity = exactextractr::exact_extract(grids$treecoverdensity, d, "mean")) %>% 
  dplyr::mutate(buildings_fact = exactextractr::exact_extract(grids$buildings, d, "max")) %>% 
  dplyr::mutate(distbuildings = exactextractr::exact_extract(grids$distbuildings, d, "mean")) %>% 
  dplyr::mutate(roads_fact = exactextractr::exact_extract(grids$roads, d, "max")) %>% 
  dplyr::mutate(distroads = exactextractr::exact_extract(grids$distroads, d, "mean")) %>% 
  dplyr::mutate(mean_precipitation = exactextractr::exact_extract(grids$mean_precipitation, d, "mean")) %>% 
  dplyr::mutate(mean_temperature = exactextractr::exact_extract(grids$mean_temperature, d, "mean")) %>% 
  dplyr::mutate(provinces = exactextractr::exact_extract(grids$provinces, d, "majority")) %>% 
  dplyr::mutate(buildings_fact = dplyr::if_else(is.na(buildings_fact), 0, buildings_fact)) %>% 
  dplyr::mutate(roads_fact = dplyr::if_else(is.na(roads_fact), 0, roads_fact)) %>% 
  dplyr::mutate_at(vars(buildings_fact, roads_fact, bin, provinces), as.factor) %>% 
  dplyr::relocate(aspect:provinces, .before = T0)


# store
saveRDS(d, "./dat/interim/spacetime_sampling.Rds") # 2 absence dates at presence locations

