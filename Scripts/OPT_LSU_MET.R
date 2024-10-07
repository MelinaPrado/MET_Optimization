###########################################################
# LSU MET Optmization and environmental matrices comparison
# Roberto Fritsche-Neto and Melina Prado
# rfneto@agcenter.lsu.edu
# Last update: jun 12 2024
###########################################################

# Cleaning environment
rm(list = ls()); ls(); gc() 

# Setting the number of cores that will be used
detectCores()
registerDoParallel(cores = detectCores()-1 ) # typing the number of cores
getDoParWorkers()

# Setting working directory
setwd("../MET_Optimization/Output/")

# Packages
library(sommer)       # Linear mixed models
library(ggplot2)      # Data visualization
require(foreach)      # Parallel execution of loops
require(doParallel)   # Parallel backend for foreach
require(doMC)         # Multicore functionality of the parallel package
library(SpATS)        # Spatial analysis of field trials with splines
library(car)          # Companion to Applied Regression
library(plyr)         # Split-apply-combine pattern in R
library(gridExtra)    # Side by side images
library(heatmaply)    # Heatmaps
library(tidyverse)    # Data manipulation
library(cluster)      # Clustering algorithms
library(factoextra)   # Clustering algorithms & visualization
library(ggrepel)      # Map plots with repelling text labels
library(maps)         # Map data
library(sf)           # Map plotting

# Loading files
data <- readRDS("results_st.RData")
clusters <- read.table("clusters.txt")
data <- merge(data, clusters)
data$cluster <- as.factor(data$cluster)

#############################
# First step - Cluster effect
#############################

# Fitting a model to estimate the cluster effect

fit.cluster <- mmer(BLUE ~ cluster,
                    random= ~ genotype + genotype:cluster,
                    weights = w,
                    rcov= ~ units,
                    data = data, 
                    verbose = FALSE)

# predicting the BLUP - main effect
BLUEs.cluster <- predict.mmer(object = fit.cluster, D = "cluster")$pvals
str(BLUEs.cluster)

# Defining limits
limits <- aes(ymax = predicted.value + std.error*1.96,
              ymin = predicted.value - std.error*1.96)

# Plotting
p <- ggplot(data = BLUEs.cluster, aes(x = factor(cluster), y = predicted.value, col = cluster)) + 
  geom_point(size = 3, shape = 19) +  
  geom_errorbar(limits, position = "identity") +
  labs(x = "Cluster", y = "Grain Yield") +
  theme_minimal() +  
  theme(axis.text.x = element_text(angle = 0, vjust = 0, size = 10))

print(p)

# Saving
ggsave(filename = "BLUES_clusters.jpg",
       plot = p,
       width = 140,
       height = 140,
       units = 'mm',
       dpi = 300
)

#########################
# Second step - MET Model
#########################

fitMET <- mmer(BLUE ~ Year + Year:Location,
               random= ~ genotype + genotype:Location,
               weights = w,
               rcov = ~ units,
               data = data, 
               verbose = FALSE)

# H2 Cullis
Vg <- fitMET$sigma$genotype
ng <- length(unique(data$genotype))
C22_g <- fitMET$PevU$genotype$BLUE
trC22_g <-sum(diag(C22_g))
av2 <- 2/ng * (trC22_g - (sum(C22_g) - trC22_g) / ng-1) 
H2.Cullis <- as.numeric(1 - av2 / (2 * Vg))
MET <- data.frame(Rep = 1,
                  Trials = 829,
                  Vg = as.numeric(fitMET$sigma$genotype),
                  Vge = as.numeric(fitMET$sigma$`genotype:Location`),
                  Vu = as.numeric(fitMET$sigma$units),
                  VgScaled = as.numeric(fitMET$sigma_scaled$genotype),
                  VgeScaled = as.numeric(fitMET$sigma_scaled$`genotype:Location`),
                  VuScaled = as.numeric(fitMET$sigma_scaled$units),
                  H.cullis = H2.Cullis, 
                  Scenario = "MET"
)

###########################
# Third step - MET_EC Model
###########################

# Model with EC
Ekinship <- readRDS("kinship.RData")

# Kronecker between Ig and Ekinship
data$GxE <- paste0(data$Location, data$genotype)

# Creating Ig
germplasms <- unique(data$genotype)
Gi <- diag(length(germplasms))
rownames(Gi) <- germplasms
colnames(Gi) <- germplasms

# Creating Kgxe
Kgxe <- kronecker(Gi, Ekinship)
location_Kgxe <- rep(rownames(Ekinship),times = length(unique(data$genotype)))
germplasm_Kgxe <- rep(rownames(Gi),each = length(unique(data$Location)))
GxE <- paste0(location_Kgxe, germplasm_Kgxe)
rownames(Kgxe) <- GxE
colnames(Kgxe) <- GxE

# Testing if it matches
all(unique(data$GxE) %in% GxE)

fitMET_EC <- mmer(BLUE ~ Year + Year:Location,
                  random= ~ genotype + vsr(GxE, Gu = Kgxe),
                  weights = w,
                  rcov = ~ units,
                  data = data, 
                  verbose = FALSE)

# H2 Cullis
Vg <- fitMET_EC$sigma$genotype
ng <- length(unique(data$genotype))
C22_g <- fitMET_EC$PevU$genotype$BLUE
trC22_g <-sum(diag(C22_g))
av2 <- 2/ng * (trC22_g - (sum(C22_g) - trC22_g) / ng-1) 
H2.Cullis <- as.numeric(1 - av2 / (2 * Vg))
MET_EC <- data.frame(Rep = 1,
                     Trials = 829,
                     Vg = as.numeric(fitMET_EC$sigma$genotype),
                     Vge = as.numeric(fitMET_EC$sigma$`u:GxE`),
                     Vu = as.numeric(fitMET_EC$sigma$units),
                     VgScaled = as.numeric(fitMET_EC$sigma_scaled$genotype),
                     VgeScaled = as.numeric(fitMET_EC$sigma_scaled$`u:GxE`),
                     VuScaled = as.numeric(fitMET_EC$sigma_scaled$units),
                     H.cullis = H2.Cullis, 
                     Scenario = "MET_EC"
)

###########################
# Fourth step - WC_MET Model
###########################

# Running all the clusters in parallel 
grid <- unique(data$cluster)

system.time(
  WC_MET <- foreach(i = 1:length(grid), 
                    .packages = c("sommer"), 
                    .combine = "rbind",
                    .export = c("mmer"),
                    .multicombine = TRUE, 
                    .errorhandling = "remove",
                    .verbose = TRUE    
  ) %dopar% {
    
    # Subsetting the data  
    sample <- droplevels.data.frame(data[data$cluster == grid[i],])
    
    # Fitting model
    fit <- mmer(BLUE ~ Year + Year:Location,
                random= ~ genotype + genotype:Location,
                weights = w,
                rcov = ~ units,
                data = sample, 
                verbose = FALSE)
    
    # H2 Cullis
    Vg <- fit$sigma$genotype
    ng <- length(unique(sample$genotype))
    C22_g <- fit$PevU$genotype$BLUE
    trC22_g <-sum(diag(C22_g))
    av2 <- 2/ng * (trC22_g - (sum(C22_g) - trC22_g) / ng-1) 
    H2 <- as.numeric(1 - av2 / (2 * Vg))      
    
    output <- data.frame(Rep = as.character(unique(sample$cluster)),
                         Trials = 829,
                         Vg = as.numeric(fit$sigma$genotype),
                         Vge = as.numeric(fit$sigma$`genotype:Location`),
                         Vu = as.numeric(fit$sigma$units),
                         VgScaled = as.numeric(fit$sigma_scaled$genotype),
                         VgeScaled = as.numeric(fit$sigma_scaled$`genotype:Location`),
                         VuScaled = as.numeric(fit$sigma_scaled$units),
                         H.cullis = H2, 
                         Scenario = "WC_MET"
    )
  }
)

WC_MET

#############################
# Fifth step - OPT_MET Model
#############################

# LSU number of trials per location - Whole network (10 years)
N_trials <- readRDS("../Data/wholeNetwork_trials.RData")

replicates <- 10

system.time(
  OPT_MET <- foreach(i = 1:replicates, 
                     .packages = c("sommer", "plyr"), 
                     .combine = "rbind",
                     .export = c("mmer", "ddply"),
                     .multicombine = TRUE, 
                     .errorhandling = "remove",
                     .verbose = TRUE    
  ) %dopar% {
    
    # Subsetting the data  
    subdata <- ddply(data, ~cluster, function(x){
      out <- droplevels.data.frame(x[x$Location %in% sample(unique(x$Location), 1, replace = F),])})
    
    fit <- mmer(BLUE ~ Year + Year:Location,
                random= ~ genotype + genotype:Location,
                weights = w,
                rcov = ~ units,
                data = subdata, 
                verbose = FALSE)
    
    # H2 Cullis
    Vg <- fit$sigma$genotype
    ng <- length(unique(subdata$genotype))
    C22_g <- fit$PevU$genotype$BLUE
    trC22_g <-sum(diag(C22_g))
    av2 <- 2/ng * (trC22_g - (sum(C22_g) - trC22_g) / ng-1) 
    H2 <- as.numeric(1 - av2 / (2 * Vg))      
    
    # Number of trials per cluster - LSU whole network
    New_N_trials <- N_trials[N_trials$Name %in% unique(subdata$Location),]
    
    output <- data.frame(Rep = i,
                         Trials = sum(New_N_trials$Number_trials),
                         Vg = as.numeric(fit$sigma$genotype),
                         Vge = as.numeric(fit$sigma$`genotype:Location`),
                         Vu = as.numeric(fit$sigma$units),
                         VgScaled = as.numeric(fit$sigma_scaled$genotype),
                         VgeScaled = as.numeric(fit$sigma_scaled$`genotype:Location`),
                         VuScaled = as.numeric(fit$sigma_scaled$units),
                         H.cullis = H2, 
                         Scenario = "OPT_MET"
    )
  }
)

OPT_MET

#################################
# Sixth step - Plotting scenarios
#################################

# Joining results
final <- rbind(MET, MET_EC, WC_MET, OPT_MET)
scen <- c("MET", "MET_EC", "WC_MET", "OPT_MET")
headers <- c("Scenario", "Mean_H.cullis", "SD_H.cullis")

# Making a table with mean and sd values
final_tab <- matrix(NA,length(scen),length(headers))
colnames(final_tab) <- headers

# Calculating mean and sd
for (i in 1:length(scen)){
  scenario <- scen[i]
  final_tab[i,1] <- scenario
  final_tab[i,2] <- mean(final[final$Scenario==scenario,]$H.cullis)
  if (is.na(sd(final[final$Scenario==scenario,]$H.cullis))){
    final_tab[i,3] <- 0
  }else{
    final_tab[i,3] <- sd(final[final$Scenario==scenario,]$H.cullis)
  }
}

# Organizing table
final_tab <- as.data.frame(final_tab)
final_tab$Mean_H.cullis <- as.numeric(final_tab$Mean_H.cullis)
final_tab$SD_H.cullis <- as.numeric(final_tab$SD_H.cullis)
final_tab$Scenario <- as.factor(final_tab$Scenario)

# Cost per scenario
final_tab$Cost <- c(67.5,67.5,67.5,15.0) 

# Heritability by cost
final_tab$Herd_Cost <- final_tab$Mean_H.cullis / final_tab$Cost

# Heritability - Plot
p <- ggplot(final_tab, aes(x = Scenario, y = Mean_H.cullis, fill = Scenario)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = Mean_H.cullis - SD_H.cullis, ymax = Mean_H.cullis + SD_H.cullis), 
                width = 0.2, position = position_dodge(0.9)) +
  ylim(0, 1) +
  labs(x = "Scenario",
       y = "Heritability")
p

# Heritability by cost - Plot
p2 <- ggplot(data = final_tab, aes(x = factor(Scenario), y = (Mean_H.cullis / Cost) *10,
                                   fill = factor(Scenario))) + 
  geom_bar(stat = "identity", position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin = Mean_H.cullis - SD_H.cullis, ymax = Mean_H.cullis + SD_H.cullis), 
                width = 0.2, position = position_dodge(0.9)) +
  ylim(0, 0.6) +
  labs(x = "Scenario", y = "Heritability / total cost") +
  scale_fill_discrete(name = "Dataset")
p2

# Joining plots
Opt_Cost <- grid.arrange(p, p2, ncol = 2)

# Saving
ggsave(filename = 'Opt_Cost.jpg',
       plot = Opt_Cost,
       #device = 'tiff',
       width = 300,
       height = 200,
       units = 'mm',
       dpi = 300)

#################################################
# Seventh step - Yield-based environmental matrix
#################################################

# Using yield by location matrix
BLUEs.MET.table <- reshape2::acast(data, genotype ~  Location, value.var = "BLUE", mean)

# Calculating correlation between environments
(cor.MET <- cor(BLUEs.MET.table, method = "pearson", use = "pairwise.complete.obs"))

# Plotting the yield based environmental matrix
heatmaply(t(BLUEs.MET.table), 
          fontsize_row = 10,
          fontsize_col = 10,
          file = "BLUEs_MET_table.png")

# Plotting the yield based environmental kinship matrix
heatmaply(cor.MET, 
          fontsize_row = 6,
          fontsize_col = 6,
          file = "cor_MET.png")

########################################
# Eighth step - Clustering and Plotting
########################################

set.seed(123)
dev.off()

# Determining optimal number of clusters
fviz_nbclust(scale(cor.MET), kmeans, method = "wss")

# K-Means Clustering
k <- kmeans(scale(cor.MET), centers = 4, nstart = 25)
p <- fviz_cluster(k, data = scale(cor.MET), geom = "point", label = FALSE, main = "") +
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

# Saving
ggsave(filename = 'GY_cluster_cor.png',
       plot = p,
       width = 400,
       height = 200,
       units = 'mm',
       dpi = 300)

# Creating data frame with locations and clusters
(clusters <- data.frame(Location = names(k$cluster), 
                        cluster = k$cluster))

# Renaming columns
colnames(clusters) <- c("Location", "cluster")
colnames(data)[16] <- "EC_clusters"

# Creating data.frame - Trials per cluster
new_data <- left_join(data,clusters)

# Matching positions
met <- read.csv("LSU_Rice_locations.csv")
met <- met[met$Location %in% unique(data$Location),]
trials2 <- merge(met, clusters)
trials2$cluster <- as.factor(trials2$cluster)
trials2$Trials <- gsub(" [tT]rials", "", trials2$Trials)
trials2$Trials <- as.numeric(trials2$Trials)

# Estimating the proportion of trials per cluster
pizza_LSU <- data.frame(Cluster = rownames(pizza_LSU), Trials = as.numeric(pizza_LSU[,1]))
capture.output(tapply(trials2$Trials, trials2$cluster, sum) / sum(trials2$Trials) * 100,
               file = "trials_per_cluster_LSU.txt")

# Creating label with percentages
label <- paste(pizza_LSU$Cluster, sprintf("(%.1f%%)", pizza_LSU$Trials))

# Plotting
pizza_LSU_graph <- ggplot(pizza_LSU, aes(x = "", y = Trials, fill = label)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  theme_void() +
  labs(fill = "Cluster (Trials %)") + 
  theme(legend.text = element_text(size = 12),  
        legend.key.size = unit(1.5, "lines"))

pizza_LSU_graph

# Saving
ggsave(filename = 'pizza_LSU_cor_graph.jpg',
       plot = pizza_LSU_graph,
       width = 200,
       height = 200,
       units = 'mm',
       dpi = 300)

############################
# Ninth step - Map plotting
############################

states_map <- map_data("state")
cities_map <- map_data("county")

# LSU Trials 
opt_trialsLSU <-ggplot() +
  geom_polygon(data = states_map, aes(x = long, y = lat, group = group), fill = "white", color = "#4A708B", size = 0.8) +
  geom_polygon(data = cities_map, aes(x = long, y = lat, group = group), fill = NA, color = "#ADD8E6", size = 0.5) +
  geom_point(data = trials2, aes(x = Long, y = Lat, color = cluster), size = 5, alpha = 0.5) +
  theme_minimal() +
  coord_cartesian(xlim = c(min(states_map$long), max(states_map$long)), 
                  ylim = c(min(states_map$lat), max(states_map$lat)))# +

opt_trialsLSU

# Saving
ggsave(filename = "opt_trialsLSU.jpg",
       plot = opt_trialsLSU,
       #device = 'tiff',
       width = 400,
       height = 200,
       units = 'mm',
       dpi = 300)

# correlation between kernels
sink("cor_between_EC_and_GxE_kernels.txt")
cor(c(Ekinship), c(cor.MET))
sink()

###############################
# Tenth step - Yield by cluster
###############################

# Preparing data
data <- left_join(data, clusters)
data$cluster <- as.factor(data$cluster)

# Cluster BLUPs
fit.cluster <- mmer(BLUE ~ cluster,
                    random= ~ genotype + genotype:cluster,
                    weights = w,
                    rcov= ~ units,
                    data = data, 
                    verbose = FALSE)

# predicting the BLUP - main effect
BLUEs.cluster <- predict.mmer(object = fit.cluster, D = "cluster")$pvals
str(BLUEs.cluster)

# Defining limits
limits <- aes(ymax = predicted.value + std.error*1.96,
              ymin = predicted.value - std.error*1.96)

p <- ggplot(data = BLUEs.cluster, aes(x = factor(cluster), y = predicted.value, col = cluster)) + 
  geom_point(size = 3, shape = 19) +  
  geom_errorbar(limits, position = "identity") +
  labs(x = "Cluster", y = "Grain Yield") +
  theme_minimal() +  
  theme(axis.text.x = element_text(angle = 0, vjust = 0, size = 10))

print(p)

# Plotting
ggsave(filename = "BLUES_cor_clusters.jpg",
       plot = p,
       #device = 'pdf',
       width = 200,
       height = 200,
       units = 'mm',
       dpi = 300
)

################################################
# Eleventh step - Comparing environment matrices
################################################

# Environment matrices
mat1 <- readRDS("W_LSU.RData")
mat2 <- cor.MET

# Transforming the matrix into a vector for plotting
value_mat1 <- as.vector(mat1)
value_mat2 <- as.vector(mat2)

# Mean, median and sd for matrix 1
mean_mat1 <- mean(value_mat1)
median_mat1 <- median(value_mat1)
sd_mat1 <- sd(value_mat1)

# Mean, median and sd for matrix 2
mean_mat2 <- mean(value_mat2)
median_mat2 <- median(value_mat2)
sd_mat2 <- sd(value_mat2)

# Calculating the density to estimate the peak position (approximate mode)
density_mat1 <- density(value_mat1)
density_mat2 <- density(value_mat2)

# Finding the peak position (x) for each density
peak_mat1 <- density_mat1$x[which.max(density_mat1$y)]
peak_mat2 <- density_mat2$x[which.max(density_mat2$y)]

# Preparing data.frame
df <- data.frame(Values = c(value_mat1, value_mat2),
                 Matrix = rep(c("Enviromic-based kernel", "Yield-based GxE matrix"), each=length(value_mat1)))

# Plotting
p <- ggplot(df, aes(x = Values, fill = Matrix)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("blue", "red")) +
  theme_minimal() +
  labs(title = "", x = "Environmental relationship/ Environmental correlation", y = "Density")

# Adding annotations for Matrix 1
p <- p + annotate("text", x = 0.25, y = 1.8, label = paste("Mean:", round(mean_mat1, 2), "\nSD:", round(sd_mat1, 2), "\nPeak:", round(peak_mat1, 2)), color = "blue", size = 4.5, hjust = 1, vjust = -0.5)

# Adding annotations for Matrix 2
p <- p + annotate("text", x = 0.80, y = 1.8, label = paste("Mean:", round(mean_mat2, 2), "\nSD:", round(sd_mat2, 2), "\nPeak:", round(peak_mat2, 2)), color = "red", size = 4.5, hjust = 0, vjust = -0.5)

# Adjusting plot
p <- p + theme(
  axis.title.x = element_text(size = 14),
  axis.title.y = element_text(size = 14),
  axis.text.x = element_text(size = 12),
  axis.text.y = element_text(size = 12),
  legend.title = element_text(size = 12),
  legend.text = element_text(size = 12)
)

print(p)

# Saving
ggsave(filename = "Density_corr_matrices.jpg",
       plot = p,
       #device = 'pdf',
       width = 300,
       height = 200,
       units = 'mm',
       dpi = 600)
