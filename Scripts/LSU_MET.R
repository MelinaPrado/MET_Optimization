########################################
# LSU MET analysis
# Roberto Fritsche-Neto and Melina Prado
# rfneto@agcenter.lsu.edu
# Last update: Jun 10 2024
########################################

# Cleaning environment
rm(list = ls()); ls(); gc() 

# Setting working directory
setwd("../MET_Optimization/Output/")

# Packages
require(foreach)      # Parallel execution of loops
require(doParallel)   # Parallel backend for foreach
require(doMC)         # Multicore functionality of the parallel package
library(SpATS)        # Spatial analysis of field trials with splines
library(car)          # Companion to Applied Regression
library(ggplot2)      # Data visualization
library(EnvRtype)     # Environmental typing for GxE analysis
library(SoilType)     # Soil typing for GxE analysis
library(caret)        # EC quality control and RFE algorithm
library(sommer)       # Linear mixed models
library(tidyverse)    # Data manipulation
library(cluster)      # Clustering algorithms
library(factoextra)   # Clustering algorithms & visualization
library(gridExtra)    # Side by side images
library(heatmaply)    # Heatmaps
library(ggrepel)      # Map plots with repelling text labels
library(maps)         # Map data
library(sf)           # Map plotting
library(ggfortify)    # Data visualization tools for statistical data using ggplot2
source('https://raw.githubusercontent.com/allogamous/EnvRtype/master/Supplementary%20Source%20and%20Data/get_weather_v2')  # Function to get weather data

# Importing data
data <- readRDS("../Data/phenotype.RData")

# Importing data from planting and harvesting date
DOP <- readRDS("../Data/DOP.RData")

# Importing metadata
metadata <- readRDS("../Data/metadata.RData")

# LSU MET Locations
trials <- unique(data$Location)

#############################################
# First step - fitting a model for each trial
#############################################

# Setting the number of cores that will be used
detectCores()
registerDoParallel(cores = detectCores()-1 ) # typing the number of cores
getDoParWorkers()

# Running all the single-trials in parallel 
grid1 <- expand.grid(trials, unique(data$Year)) #Just 2021 and 2022 

# Calculating BLUEs, BLUPs and standard errors for each Location x Year combination

system.time( #time to execute the code
  results.st <- foreach(i = 1:nrow(grid1), #iterate over each row of grid1
                        .packages = c("SpATS", "car"), 
                        .combine = "rbind",
                        .export = c("SpATS", "predict.SpATS", "getHeritability", "outlierTest"),
                        .multicombine = TRUE, 
                        .errorhandling = "remove",
                        .verbose = TRUE    
  ) %dopar% { #each iteration run on a separate core
    
    # Subsetting trials (By Location Name and Year)
    sample <- droplevels.data.frame(data[data$Location == grid1[i,1] & data$Year == grid1[i,2],])
    
    # Detecting and eliminating outlier (Yield)
    fit <- lm(Yield ~ replicate + rowNumber + colNumber + genotype + Gen_type, data = sample)
    outlier <- names(outlierTest(fit)$p)
    sample[outlier, "Yield"] <- NA
    
    # Storing number of rows and columns (PSANOVA)
    sample$row <- as.numeric((sample$rowNumber))
    sample$col <- as.numeric((sample$colNumber))
    nrow <- max(sample$row)
    ncol <- max(sample$col)
    nseg.row <- nrow
    nseg.col <- ncol
    
    # Fitting model with spatial correction and fixed genotype
    fitF <- SpATS(response = "Yield", 
                  fixed = ~ 1, 
                  random = ~ replicate + rowNumber + colNumber, 
                  spatial = ~ PSANOVA(col, row, nseg = c(nseg.col, nseg.row)), 
                  genotype = "genotype", 
                  genotype.as.random = FALSE, 
                  data = sample)
    
    # Estimate BLUEs
    blues <- predict.SpATS(fitF, which = "genotype")
    blues <- blues[,c("genotype", "predicted.values", "standard.errors")] 
    colnames(blues)[2:3] <- c("BLUE", "sep_BLUE")
    
    # Fitting model with spatial correction and random genotype
    fitR <- SpATS(response = "Yield", 
                  fixed = ~ 1, 
                  random = ~ replicate + rowNumber + colNumber, 
                  spatial = ~ PSANOVA(col, row, nseg = c(nseg.col, nseg.row)), 
                  genotype = "genotype", 
                  genotype.as.random = TRUE, 
                  data = sample)
    
    # Obtaining the heritability - package function
    h2g <- getHeritability(fitR)
    
    # Broad-sense heritability - Cullis method
    Vg <- fitR$var.comp["genotype"]
    ng <- length(unique(sample$genotype))
    C11_g <- fitR$vcov$C11_inv
    trC11_g <-sum(diag(C11_g))
    av2 <- 2/ng * (trC11_g - (sum(C11_g) - trC11_g) / ng-1) # mean var of a difference between genotypic BLUPS
    H2.Cullis <- 1 - av2 / (2 * Vg)
    
    # Estimate BLUPs for Grain Yield
    blups <- predict.SpATS(fitR, which = "genotype")
    blups <- blups[,c("predicted.values", "standard.errors")]
    colnames(blups)[1:2] <- c("BLUP", "sep_BLUP")
    
    # Reliability
    rel <- mean(1 - blups$sep_BLUP^2 / fitR$var.comp["genotype"])
    
    # weights for ID's - adjust residual for further analysis
    vcov.mme <- fitR$vcov$C11_inv
    w <- diag(vcov.mme)
    
    # Storing data
    output <- data.frame(blues,
                         w = w,
                         Location = as.character(unique(sample$Location)), 
                         Year = as.character(unique(sample$Year)),
                         H.cullis = H2.Cullis,
                         trait = fitR$model$response[1]
    )
  }
)

# merge metadata and remove fillers
results.st <- merge(results.st, metadata)
saveRDS(results.st, file = "results_st.RData")

# Multi environment model - Location BLUP
fitMET <- mmer(BLUE ~ Year,
               random= ~ Location + genotype + genotype:Location,
               weights = w,
               rcov = ~ units,
               data = results.st, 
               verbose = FALSE)

# Utilizing BLUP for each location
BLUPs.e <- predict.mmer(object = fitMET, D = "Location")$pvals

########################
# Second step - envRtype
########################

# matching positions
met <- readRDS("../Data/LSU_Rice_locations.RData")
met <- met[met$Location %in% unique(results.st$Location),]

# Considering historical weather data
period <- 2011:2020 
trial <- met$Location 

# Loop index
grid2 <- expand.grid(trial, period)
colnames(grid2) <- c("Location", "period")
grid2 <- merge(grid2, met)
grid2 <- merge(grid2, DOP)

# Rewriting dates using a single regular expression to remove "2022" and "2021"
grid2$plantingDate <- gsub("202[12]", "", grid2$plantingDate, fixed = FALSE)
grid2$harvestDate <- gsub("202[12]", "", grid2$harvestDate, fixed = FALSE)

# Environmental covariates using EnvRtype package
system.time(
  env.data <- foreach(i = 1:nrow(grid2), 
                      .packages = c("EnvRtype"), 
                      .combine = "rbind",
                      .export = c("get_weather"),
                      .multicombine = TRUE, 
                      .errorhandling = "remove",
                      .verbose = TRUE    
  ) %dopar% {
    
    # Sampling by grid2
    sample <- grid2[i,]  
    
    # Imports weather data from NASA-POWER GIS
    output <- get_weather(env.id = sample$Location,
                          country = NULL,
                          lat = sample$Lat,
                          lon = sample$Long,
                          start.day = paste0(sample$period, sample$plantingDate),
                          end.day = paste0(sample$period, sample$harvestDate),
                          parallel = F)
  }  
)  

# Tuning for cardinal of temperature - obtaining other traits
aux <- met[,c("Location","Altitude.m.")] #Altitude
colnames(aux) <- c("env", "Alt")
aux2 <- merge(env.data, aux)

# Calculating a series of parameters based on the get_weather() object
df.clim <- processWTH(env.data = aux2,Tbase1 = 12,Tbase2 = 24,Topt1 = 33,Topt2 = 37, Alt = aux2$Alt)

# Covariates that will be used in the analysis
(var.i = names(df.clim)[c(9:16,22:29)])

# Phenological tages that will be used in the analysis
stages = c('EM_MAX.TIL',"MAX.TIL_PAN.INIT", "PAN.INIT_PRE.FLW", "PRE.FLW_FLW", "FLW_POST.FLW", "POST.FLW_MAT")

# Intervals that will be used in the analysis
interval <- c(0, 45, 60, 75, 90, 105)

# Building the matrix of environmental covariates using reaction norms
E <- W_matrix(env.data = df.clim, #data.frame of environmental variables gotten from get_weather()
               env.id = 'env', #name of the columns to be used as id for environments
               var.id = var.i, #Covariates
               statistic = 'mean', # Indicates what statistic must be runned
               scale = F, 
               center = F,
               by.interval = T, #Indicates if temporal intervals must be computed insied of each environment.
               sd.tol = 5,
               time.window = interval, 
               names.window = stages)

E[E == "NaN"] <- 0
E <- scale(E) 

#######################
# Third step - SoilType
#######################

# load isric dataset
data(soil.data) 

# Running in parallel
system.time(
  SC <- foreach(i = 1:nrow(met), 
                .packages = c("caret", "stringr"), 
                .combine = "rbind",
                .export = c("predict", "train", "trainControl", "str_replace_all"),
                .multicombine = TRUE, 
                .errorhandling = "remove",
                .verbose = TRUE    
  ) %dopar% {
    
    # Subsetting data  
    trial.MET <- droplevels.data.frame(met[i,])
    
    # Retrieving data
    output <- get_soil(env.id = trial.MET$Location, 
                       lat = trial.MET$Lat, 
                       long = trial.MET$Long, 
                       max.depth = 20, 
                       isric.data = soil.data)
  }   
)

# Soil Covariates visualization
SCov <- reshape2::dcast(SC, env ~ Trait, value.var = "Predicted.Value", mean)
rownames(SCov) <- SCov[,1]
SCov <- SCov[,2:ncol(SCov)]

# Eliminating soil traits without any information
SCov <- SCov[, apply(SCov, 2, function(x){!all(is.na(x))})]
SCov[is.na(SCov)] <- 0
SCov <- scale(SCov)

# Joining E with Scov -- Forming W matrix
E <- E[match(rownames(SCov), rownames(E)),]
all(rownames(E) == rownames(SCov))
W <- cbind(E, SCov)

#####################################################
# Fourth step - Quality Control and feature selection 
#####################################################

# Eliminating those EC almost perfected correlated
highlyCorrelated <- findCorrelation(cor(W), cutoff = 0.95)
print(highlyCorrelated)
W.clean <- W[,-highlyCorrelated]
W.clean[W.clean == "NaN"] <- 0
W.clean[W.clean == "Na"] <- 0

# Identifying the most important EC via RFE to explain the GY across locations
data.rfe <- data.frame(W.clean[match(BLUPs.e$Location, rownames(W.clean)),], GY = BLUPs.e$predicted.value)
data.rfe <- data.rfe[!is.na(data.rfe[,2]),]
set.seed(29121983)

# Generating a control object that can be used to specify the details of the feature selection algorithms
control <- rfeControl(functions = rfFuncs, method = "repeatedcv", number = 5, repeats = 5)
set.seed(29121983)
opt.reg <- rfe(x = data.rfe[, 1:ncol(W.clean)], y = data.rfe[, ncol(data.rfe)], sizes = c(1:10), 
               metric = "Rsquared", maximize = TRUE, rfeControl = control)
opt.reg
predictors(opt.reg)

# Saving the most important EC to explain GY via ML
sink("themostimportantECtoexplainGYviaML.txt")
opt.reg
sink()

# Saving predictors
sink("EC.predictors.txt")
data.frame(EC = predictors(opt.reg))
sink()

# Cleaning W again using RFE
W.clean <- W.clean[, predictors(opt.reg)]

# Saving W amtrix
saveRDS(W.clean, "W_LSU.RData")

######################################
# Fifth step - Clustering and Plotting
######################################

# Plotting EC heatmap
heatmaply(W.clean, 
          fontsize_row = 6,
          fontsize_col = 6,
          file = "EC_heatmap.png")

# Calculating and plotting environmental kinship
kinship <- env_kernel(env.data = as.matrix(W.clean), is.scaled = T, gaussian = T, sd.tol = 5)[[2]]
saveRDS(kinship, "kinship.RData")

heatmaply(kinship, 
          fontsize_row = 6,
          fontsize_col = 6,
          file = "EC_kinship_heatmap.png")

# Clustering

# Determining optimal number of clusters
set.seed(123)
dev.off()
fviz_nbclust(scale(W.clean), kmeans, method = "wss")+ labs(title = NULL)

# K-Means Clustering
k <- kmeans(scale(W.clean), centers = 5, nstart = 25)

# Visualizing clustering results
p <- fviz_cluster(k, data = scale(W.clean), geom = "point", label = FALSE, main = "") +
  geom_label_repel(data = k$data, aes(label = name),
                   box.padding = unit(0.5, "lines"),  
                   fill = "grey",  
                   color = "white",  
                   fontface = "bold", 
                   size = 3,  
                   label.size = NA  
  ) +
  theme(text = element_text(size = 15))
p

ggsave(filename = "ClustersLSU.jpg",
       plot = p,
       #device = 'tiff',
       width = 400,
       height = 200,
       units = 'mm',
       dpi = 300)

# Creating data frame with locations and clusters
(clusters <- data.frame(Location = names(k$cluster), 
                        cluster = k$cluster))
write.table(clusters, "clusters.txt")

# Randomly choosing a location to show MET optimization
opt.trials <- matrix(NA, nrow = length(k$size), ncol = 1)
for (i in 1:length(k$size)) {
  opt.trials[i,] <- sort(
    sample(clusters$Location[clusters$cluster == i], 1, replace = F))
}  

# Creating new data.frame with clusters and MET information
trials2 <- merge(met, clusters)
trials2$cluster <- as.factor(trials2$cluster)
trials2$Trials <- gsub(" [tT]rials", "", trials2$Trials)

# Estimating the proportion of trials per cluster - LSU MET
trials2$Trials <- as.numeric(trials2$Trials)
pizza_LSU <- as.data.frame(tapply(trials2$Trials, trials2$cluster, sum) / sum(trials2$Trials) * 100)
pizza_LSU <- data.frame(Cluster = rownames(pizza_LSU), Trials = as.numeric(pizza_LSU[,1]))
capture.output(tapply(trials2$Trials, trials2$cluster, sum) / sum(trials2$Trials) * 100,
               file = "trials_per_cluster_LSU.txt")
opt.trials <- trials2[trials2$Location %in% unique(c(opt.trials)), ]

# Creating label with percentages
label <- paste(pizza_LSU$cluster, sprintf("(%.1f%%)", pizza_LSU$Trials))

# Plotting
pizza_LSU_graph <- ggplot(pizza_LSU, aes(x = "", y = Trials, fill = label)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  theme_void() +
  labs(fill = "Cluster (Trials %)") +  
  theme(legend.text = element_text(size = 12),  
        legend.key.size = unit(1.5, "lines"))  

# Saving Plot
ggsave(filename = 'pizza_LSU_graph.png',
       plot = pizza_LSU_graph,
       #device = 'tiff',
       width = 200,
       height = 200,
       units = 'mm',
       dpi = 300)

###########################
# Sixth step - Map plotting
###########################

# ----- LSU Closest trials - Trials optimization
opt.trials_closest <- c("Mamou","Winnsboro","Rice Research Station",
                        "Palmetto","Wintermann Rice Research Station")
opt.trials <- trials2[trials2$Location %in% unique(c(opt.trials_closest)), ]

# Getting USA states and counties
states_map <- map_data("state")
cities_map <- map_data("county")

# LSU Opt trials 
opt_trialsLSU <- ggplot() +
  geom_polygon(data = states_map, aes(x = long, y = lat, group = group), fill = "white", color = "#4A708B", size = 0.8) +
  geom_polygon(data = cities_map, aes(x = long, y = lat, group = group), fill = NA, color = "#ADD8E6", size = 0.5) +
  geom_point(data = opt.trials, aes(x = Long, y = Lat, color = cluster), size = 6, alpha = 0.6) +
  theme_minimal() +
  coord_cartesian(xlim = c(min(states_map$long), max(states_map$long)), 
                  ylim = c(min(states_map$lat), 37))
opt_trialsLSU

# LSU Trials 
trialsLSU <-ggplot() +
  geom_polygon(data = states_map, aes(x = long, y = lat, group = group), fill = "white", color = "#4A708B", size = 0.8) +
  geom_polygon(data = cities_map, aes(x = long, y = lat, group = group), fill = NA, color = "#ADD8E6", size = 0.5) +
  geom_point(data = trials2, aes(x = Long, y = Lat, color = cluster), size = 6, alpha = 0.6) +
  theme_minimal() +
  coord_cartesian(xlim = c(min(states_map$long), max(states_map$long)), 
                  ylim = c(min(states_map$lat), 37)) 
trialsLSU

# Combining plots
opt_LSU_trials <- grid.arrange(trialsLSU, opt_trialsLSU, ncol = 1)

ggsave(filename = "opt_LSU_trials.jpg",
       plot = opt_LSU_trials,
       #device = 'tiff',
       width = 400,
       height = 200,
       units = 'mm',
       dpi = 300)

###########################
# Seventh step - PCA Biplot
###########################

# Creating data.frame to plot environmental characterization
W.clean <- as.data.frame(W.clean)
W.clean$Location <- rownames(W.clean)
W.clean <- left_join(W.clean, clusters)
cluster_ECs <- matrix(NA, length(unique(W.clean$cluster)),8)

# Naming variables to plot
rownames(cluster_ECs) <- c("CLUSTER 1","CLUSTER 2","CLUSTER 3","CLUSTER 4","CLUSTER 5")
colnames(cluster_ECs) <- c(colnames(W.clean)[1:8])
for(i in 1:length(unique(W.clean$cluster))){
  data <- W.clean[W.clean$cluster==i,1:8]
  for(j in 1:8){
    cluster_ECs[i,j] <- mean(data[,j])
  }
}

#Performs a principal components analysis on the matrix 
pca_resultado <- prcomp(cluster_ECs, scale. = TRUE)

# Plotting
p <- autoplot(pca_resultado, data = cluster_ECs, loadings = TRUE, 
              loadings.label = TRUE, loadings.label.size = 5, loadings.colour = 'black',
              loadings.label.repel = TRUE) +
  geom_text_repel(aes(label = rownames(cluster_ECs), colour = rownames(cluster_ECs)),
                  size = 5,      
                  point.padding = unit(0.5, "lines")) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black") +
  theme(legend.position = "none")
p

# Saving
ggsave(filename = "LSU_cluster_ECs.jpg",
       plot = p,
       width = 250,
       height = 250,
       units = 'mm',
       dpi = 600)


