########################################
# USA Rice Belt TPE analysis
# Roberto Fritsche-Neto and Melina Prado
# rfneto@agcenter.lsu.edu
# Last update: jun 12 2024
########################################

# Cleaning environment
rm(list = ls()); ls(); gc() 

# Setting working directory
setwd("../MET_Optimization/Output/")

# Packages
require(foreach)      # Parallel execution of loops
require(doParallel)   # Parallel backend for foreach
require(doMC)         # Multicore functionality of the parallel package
library(ggplot2)      # Data visualization
require(EnvRtype)     # Environmental typing for GxE analysis
library(SoilType)     # Soil typing for GxE analysis
library(tidyverse)    # Data manipulation
library(cluster)      # Clustering algorithms
library(factoextra)   # Clustering algorithms & visualization
library(gridExtra)    # Side by side images
library(heatmaply)    # Heatmaps
library(ggrepel)      # Map plots with repelling text labels
library(ggplot2)      # Data visualization
library(maps)         # Map data
library(sf)           # Map plotting
library(data.table)   # Data.frame manipulation
library(ggfortify)    # Data visualization tools for statistical data using ggplot2
library(ggrepel)      # Map plots with repelling text labels
library(dplyr)        # Data manipulation
source('https://raw.githubusercontent.com/allogamous/EnvRtype/master/Supplementary%20Source%20and%20Data/get_weather_v2')

# Setting the number of cores that will be used
detectCores()
registerDoParallel(cores = detectCores()-1 ) # typing the number of cores
getDoParWorkers()

# Loading files
data <- readRDS("../Data/dataset_usa_rice.RData")

########################
# First step - envRtype
########################

# Considering historical weather data
period <- 2011:2020

# Considering all counties
trial <- unique(data$state_county)

# Loop index
grid <- expand.grid(trial, period)
colnames(grid) <- c("state_county", "year")
grid2 <- merge(grid, data)

# Combining the months of planting with the years 2011-2020, instead of 2023

# picking only the last 6 data digits
substrRight <- function(x, n = 6){
  substr(x, nchar(x)-n+1, nchar(x)) #substr(x, start, stop)
}

# Creating historical data
grid2$planting <- as.character(paste0(grid2$year, apply(data.frame(grid2$planting), 1, substrRight)))
grid2$harvest <- as.character(paste0(grid2$year, apply(data.frame(grid2$harvest), 1, substrRight)))

# # Environmental covariates using EnvRtype package
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
    output <- get_weather(env.id = sample$state_county,
                          country = NULL,
                          lat = sample$lat,
                          lon = sample$lon,
                          start.day = sample$planting,
                          end.day = sample$harvest,
                          parallel = F,
                          temporal.scale = "daily")
  }  
)  

# Tuning for cardinal of temperature, and obtaining other traits
aux <- data
colnames(aux)[4] <- "env"
aux2 <- merge(env.data, aux)

#  Calculating a series of parameters based on the get_weather() object
df.clim <- processWTH(env.data = aux2, Tbase1 = 12, Tbase2 = 24, Topt1 = 33, Topt2 = 37, Alt = aux2$elev)

# Covariates that will be used in the analysis
(var.i = names(df.clim)[c(9:16,35:42)]) 

# Phenological tages that will be used in the analysis
stages = c('EM_MAX.TIL',"MAX.TIL_PAN.INIT", "PAN.INIT_PRE.FLW", "PRE.FLW_FLW", "FLW_POST.FLW", "POST.FLW_MAT")

# Intervals that will be used in the analysis
interval <- c(0, 45, 60, 75, 90, 105) 

# Building the matrix of environmental covariates using reaction norms
E <- W_matrix(env.data = df.clim, #data.frame of environmental variables gotten from get_weather()
               env.id = 'env', #name of the columns to be used as id for environments
               var.id = var.i, #Covariates
               statistic = 'mean',# Indicates what statistic must be runned
               scale = F, 
               center = F,
               by.interval = T, #Indicates if temporal intervals must be computed insied of each environment.
               sd.tol = 5,
               time.window = interval, 
               names.window = stages)

E[E == "NaN"] <- 0
E <- scale(E) 

#######################
# Second step - SoilType
#######################

# load isric dataset
data(soil.data) 

# Running in parallel
system.time(
  SC <- foreach(i = 1:nrow(data), 
                .packages = c("caret", "stringr"), 
                .combine = "rbind",
                .export = c("predict", "train", "trainControl", "str_replace_all"),
                .multicombine = TRUE, 
                .errorhandling = "remove",
                .verbose = TRUE    
  ) %do% {
    
    # Subsetting data
    sample <- data[i,]
    
    # Retrieving data
    output <- get_soil(env.id = sample$state_county,
                       lat = sample$lat,
                       lon = sample$lon,
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
W_USA <- cbind(E, SCov)

# Pick only the predictors used in the LSU MET
EC.predictors <- c(read.table("EC.predictors.txt")[,1])
W_USA <- W_USA[, colnames(W_USA) %in% EC.predictors]
W_USA <- scale(W_USA)

#############################################
# Third step - Creating South and LA dataset
#############################################

# South dataset - USA dataset without California trials
W_South <- W_USA[!rownames(W_USA) %like% "california",]
W_South <- scale(W_South)

# LA dataset - USA dataset only with Louisiana trials
W_LA <- W_USA[rownames(W_USA) %like% "louisiana",]
W_LA <- scale(W_LA)

######################################
# Fourth step - Clustering and Plotting
######################################

# Loop - Different datasets
datasets <- c("USA", "South", "LA")

# Using the WSS method to aid in deciding the number of clusters 
n.clusters <- c(2, 4, 2)
x <- 2
for (x in 1:3) {
  
  # Loop dataset
  
  W_loop <- paste0("W_",datasets[x])
  W <- get(W_loop)
  
  # Plotting EC heatmap
  heatmaply(W, 
            fontsize_row = 6,
            fontsize_col = 6,
            file = paste0(W_loop,".png"))
  
  # Calculating and plotting environmental kinship
  kinship <- env_kernel(env.data = as.matrix(W), is.scaled = T, gaussian = T, sd.tol = 5)[[2]]
  
  # Plotting heatmap
  heatmaply(kinship, 
            fontsize_row = 6,
            fontsize_col = 6,
            file = paste0("EC_kinship_",datasets[x],".png"))
  
  # Clustering
  # Determining optimal number of clusters
  set.seed(123)
  dev.off()
  fviz_nbclust(scale(W), kmeans, method = "wss")+ labs(title = NULL)
  
  # K-Means Clustering
  k <- kmeans(scale(W), centers = n.clusters[x], nstart = 25)
  
  # Visualizing clustering results
  p <- fviz_cluster(k, data = scale(W), geom = "point", label = FALSE, main = "") +
    geom_label_repel(data = k$data, aes(label = name),
                     box.padding = unit(0.5, "lines"),  
                     fill = "grey",  
                     color = "white",  
                     fontface = "bold",  
                     max.overlaps = Inf,
                     size = 3,  
                     label.size = NA  
    ) +
    theme(text = element_text(size = 15))
  p
  
  # Saving
  ggsave(filename = paste0("Clusters_", datasets[x], ".png"),
         plot = p,
         #device = 'tiff',
         width = 400,
         height = 200,
         units = 'mm',
         dpi = 300)
  
  # Creating data frame with locations and clusters
  (clusters <- data.frame(state_county = names(k$cluster), 
                          cluster = k$cluster))
  
  # Creating data.frame - Trials per cluster
  trials <- merge(clusters, data)
  trials$cluster <- as.factor(trials$cluster)
  
  # Estimating the proportion of trials per cluster
  pizza <- data.frame(Cluster = rownames(pizza), Trials = as.numeric(pizza[,1]))
  capture.output(tapply(trials$acreage, trials$cluster, sum) / sum(trials$acreage) * 100,
                 file = paste0("trials_per_cluster_", datasets[x], ".txt"))

  # Creating label with percentages
  label <- paste(pizza$Cluster, sprintf("(%.1f%%)", pizza$Trials))
  
  # Plotting - Pizza graph
  pizza_graph <- ggplot(pizza, aes(x = "", y = Trials, fill = label)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar(theta = "y") +
    theme_void() +
    labs(fill = "Cluster (Trials %)") +  
    theme(legend.text = element_text(size = 12),  
          legend.key.size = unit(1.5, "lines")) 
  pizza_graph
  
  # Saving plot
  ggsave(filename = paste0("pizza_",datasets[x],"_graph.png"),
         plot = pizza_graph,
         #device = 'tiff',
         width = 200,
         height = 200,
         units = 'mm',
         dpi = 300)
  
  
  ###########################
  # Fifth step - Map plotting
  ###########################
  
  # Getting USA states and counties
  states_map <- map_data("state")
  cities_map <- map_data("county")
  
  # Creating data.frame
  trials2 <- trials
  colnames(trials2)[c(1,10,11)] <- c("Location","Long", "Lat")
  
  # USA Opt trials 
  map <- ggplot() +
    geom_polygon(data = states_map, aes(x = long, y = lat, group = group), fill = "white", color = "#838B83", size = 0.8) +
    geom_polygon(data = cities_map, aes(x = long, y = lat, group = group), fill = NA, color = "#C1CDC1", size = 0.5) +
    geom_point(data = trials2, aes(x = Long, y = Lat, color = cluster), size = 5, alpha = 0.6) +
    theme_minimal() +
    coord_cartesian(xlim = c(min(states_map$long), max(states_map$long)), 
                    ylim = c(min(states_map$lat), max(states_map$lat))) +
    theme(legend.position = "none") 
  map
  
  # Saving
  ggsave(filename = paste0("trials_", datasets[x], ".png"),
         plot = map,
         #device = 'tiff',
         width = 400,
         height = 200,
         units = 'mm',
         dpi = 300)
  
  ###########################
  # Sixth step - PCA Biplot
  ###########################
  
  # Creating data.frame to plot environmental characterization 
  cluster_ECs <- matrix(NA, length(unique(trials2$cluster)),8)
  
  # Naming variables to plot
  rownames(cluster_ECs) <- paste("CLUSTER", seq_len(n.clusters[x]))
  colnames(cluster_ECs) <- c(colnames(W)[1:8])
  
  # Calculates the clusters means
  for(i in 1:length(unique(trials2$cluster))){
    subset_trial <- trials2[trials2$cluster==i,]
    data <- W[rownames(W) %in% subset_trial$Name,1:8]
    for(j in 1:8){
      cluster_ECs[i,j] <- mean(data[,j])
    }
  }
  
  #Performs a principal components analysis on the matrix 
  pca_resultado <- prcomp(cluster_ECs, scale. = TRUE)
  
  # Plotting
  p <- autoplot(pca_resultado, data = cluster_ECs, loadings = TRUE, 
                loadings.label = TRUE, loadings.label.size = 5, loadings.colour = 'black',
                loadings.label.repel = TRUE)+
    geom_text_repel(aes(label = rownames(cluster_ECs), colour = rownames(cluster_ECs)),
                    size = 5,      
                    point.padding = unit(0.5, "lines")) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "black") +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "black") +
    theme(legend.position = "none")
  p
  
  # Saving
  ggsave(filename = paste0(datasets[x],"_cluster_ECs.png"),
         plot = p,
         #device = 'tiff',
         width = 250,
         height = 250,
         units = 'mm',
         dpi = 600)
  
}

########################################
# Seventh step - Comparing to LSU trials
########################################

# Calculating centroids
centroids <- k$centers

# Scaling new locations
W_LSU <- readRDS("../Output/W.clean.RData")
new_data_scaled <- scale(W_LSU)

# Distance - New points to centroids
distances <- as.matrix(dist(rbind(centroids, new_data_scaled)))[1:nrow(centroids), (nrow(centroids) + 1):(nrow(centroids) + nrow(new_data_scaled))]

# Atributte clusters
cluster_assignments <- apply(distances, 2, which.min)

# Count number os locations per cluster
cluster_counts <- table(cluster_assignments)

# Percentage
total_locations <- length(cluster_assignments)
cluster_percentages <- (cluster_counts / total_locations) * 100

# Creating data.frames
cluster_df <- data.frame(
  Cluster = names(cluster_counts),
  Trials = as.numeric(cluster_percentages)
)

# Exibir o data.frame
print(cluster_df)

#write.csv(cluster_df, file = "cluster_LSU_TPE.csv")

# Pizza plot
label <- paste(cluster_df$Cluster, sprintf("(%.1f%%)", cluster_df$Trials))

pizza_LSU_TPE_graph <- ggplot(cluster_df, aes(x = "", y = Trials, fill = label)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  theme_void() +
  labs(fill = "Percentage of LSU trials by TPE.") +  
  theme(legend.text = element_text(size = 12),  
        legend.title = element_text(size = 16),
        legend.key.size = unit(1.5, "lines"))  
pizza_LSU_TPE_graph

ggsave(filename = 'pizza_LSU_TPE_graph.jpg',
       plot = pizza_LSU_TPE_graph,
       #device = 'tiff',
       width = 200,
       height = 200,
       units = 'mm',
       dpi = 600)

cluster_assignments <- as.data.frame(cluster_assignments)
#write.csv(cluster_assignments, file = "cluster_assignments.csv")