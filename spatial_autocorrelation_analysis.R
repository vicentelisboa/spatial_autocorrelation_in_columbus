# Final project 
# Name: Vicente Lisboa
# Project: Analysis of Spatial Autocorrelation in neighborhoods of Columbus, OH 


# Import libraries -------------

pacman::p_load(sf,spData,tidyverse,dplyr,reshape2,haven,stringr,lubridate,openxlsx,data.table,
               readxl,gdata,haven,miceadds,ezR,janitor,forcats,scales,ggplot2,cowplot)

# Import data ---------

#Import columbus data from spData package

columbus_df <- sf::st_read(system.file("shapes/columbus.shp", package="spData")[1])
columbus <- as(columbus_df, "Spatial")

#Import weights
col.gal.nb <- read.gal(system.file("weights/columbus.gal", package="spData")[1])


# 2 Data analysis ----------------

# Resume table 

data_analysis <- data.frame (variable  = c("area", "hoval","inc","open","crime"),
                             min_value  = c(min(columbus_df$AREA), min(columbus_df$HOVAL),min(columbus_df$INC),min(columbus_df$OPEN),min(columbus_df$CRIME)),
                             mean_value  = c(mean(columbus_df$AREA), mean(columbus_df$HOVAL),mean(columbus_df$INC),mean(columbus_df$OPEN),mean(columbus_df$CRIME)),
                             max_value  = c(max(columbus_df$AREA), max(columbus_df$HOVAL),max(columbus_df$INC),max(columbus_df$OPEN),max(columbus_df$CRIME)))

data_analysis

#In the report I will copy it in excel just for the space and format

#Relation between HOVAL and CRIME

cri_hov <- ggplot(columbus_df, aes(x=CRIME, y=HOVAL)) +
  geom_line()+
  theme_minimal() +
  labs(
    x = "Crime rate (per thousand households)" ,
    y = "Housing value (in 1,000 USD)") 

cri_inc <- ggplot(columbus_df, aes(x=CRIME, y=INC)) +
  geom_line()+
  theme_minimal()+
  labs(
    x = "Crime rate (per thousand households)" ,
    y = "Household Income (in 1,000 USD)") 


plot_row <- plot_grid(cri_hov, cri_inc)

# Now add the title
title <- ggdraw() + 
  draw_label(
    "Housing value and Household income against crime rate",
    fontface = 'bold',
    x = 0,
    hjust = -0.3
  ) +
  theme(
    plot.margin = margin(0, 0, 0, 7)
  )

#Plot graph
plot_grid(
  title, plot_row,
  ncol = 1,
  rel_heights = c(0.1, 1)
)


#Crime rates graph 
#Map with the crime rates by neighborhood

cr_graph <- tm_shape(columbus) + tm_polygons(col = "CRIME", title = "Crimen Rates", style = "quantile") +
  tm_layout(legend.outside = TRUE)
cr_graph



# 3.1 Global Moran's I -------------

#nb2listw function supplements a neighbours list with spatial weights for the chosen coding scheme 

nbw <- spdep::nb2listw(col.gal.nb, style = "W")

#Calculating gmoran statistic
(gmoran <- moran.test(columbus$CRIME, nbw))

#Saving obtained values
moran_I_statistic = gmoran[["estimate"]][["Moran I statistic"]] 
moran_I_statistic_sd = gmoran[["statistic"]] 
moran_I_statistic_pvalue = gmoran[["p.value"]]

#Table with main values 

table_moran_I <- data.frame ('Moran I statistic' = moran_I_statistic,
                  "Moran I statistic standard deviate" = moran_I_statistic_sd,
                  'p-value'= moran_I_statistic_pvalue)

table_moran_I
#For the final report an excel table will be shown with all the parameters obtained


# 3.1.1 Montecarlo approach ---------------

#Monte Carlo approach that creates random patterns by reassigning the values among the fixed areas and calculates the Moran’s I for each of these patterns
gmoranMC <- moran.mc(columbus_df$CRIME, nbw, nsim = 999)
gmoranMC

#Creating histogram
hist(gmoranMC$res, main = "Montecarlo simulation", xlab = 'Value')
abline(v = gmoranMC$statistic, col = "red")


#Scatter plot
#Moran’s I scatterplot to visualize the spatial autocorrelation in the data
moran.plot(columbus_df$CRIME, nbw,xlab = 'Crime rate', ylab = 'Spatially lagged values - crime rate')


# 3.2 Local Moran’s I ------------

#Computing the  Local Moran's

lmoran <- localmoran(columbus_df$CRIME, nbw, alternative = "greater")
head(lmoran)


columbus$lmI <- lmoran[, "Ii"] # local Moran's I
columbus$lmZ <- lmoran[, "Z.Ii"] # z-scores
columbus$lmp <- lmoran[, "Pr(z > E(Ii))"] # p-values corresponding to alternative greater

#Local Moran's I plot 

local_moran_graph <- tm_shape(columbus) + tm_polygons(col = "lmI", title = "Local Moran's I", style = "quantile",midpoint = NA) +
  tm_layout(legend.title.size = 1,legend.text.size = 0.5)

local_moran_graph

# Maps with z-scores and p-values as follows.
z_score_graph <- tm_shape(columbus) + tm_polygons(col = "lmZ", title = "Z-score", breaks = c(-Inf, 1.65, Inf)) +
  tm_layout(legend.title.size = 1,legend.text.size = 0.6)

p_value_graph <- tm_shape(columbus) + tm_polygons(col = "lmp", title = "p-value", breaks = c(-Inf, 0.05, Inf)) +
  tm_layout(legend.title.size = 1,legend.text.size = 0.6)

tmap_arrange(z_score_graph, p_value_graph)


#Two sided hypothesis test 

#z-score values lower than -1.96 would indicate negative spatial autocorrelation, 
#and z-score values greater than 1.96 would indicate positive spatial autocorrelation.

hypo_test <- tm_shape(columbus) + tm_polygons(col = "lmZ",
title = "Local Moran's I", style = "fixed",
breaks = c(-Inf, -1.96, 1.96, Inf),
labels = c("Negative SAC", "Random SAC", "Positive SAC"),
palette = "RdBu", n = 3, midpoint = NA, border.alpha = 0.3) +
  tm_layout(legend.title.size = 1,legend.text.size = 0.6)

hypo_test


# 3.3 Clustering -----------

#Clustering classification 

lmoran <- localmoran(columbus_df$CRIME, nbw, alternative = "two.sided")
#Showing principal values 
head(lmoran)

# Saving p-values from column 5
columbus$lmp <- lmoran[, 5] 

#Plot with lagged and scale data (removed from report because larga of the document)
(mp <- moran.plot(as.vector(scale(columbus$CRIME)), nbw))

# Find out quadrant of (observation, lagged_observation) and significance
columbus$quadrant <- NA
columbus@data[(mp$x >= 0 & mp$wx >= 0) & (columbus$lmp <= 0.05), "quadrant"] <- 1
columbus@data[(mp$x <= 0 & mp$wx <= 0) & (columbus$lmp <= 0.05), "quadrant"] <- 2
columbus@data[(mp$x >= 0 & mp$wx <= 0) & (columbus$lmp <= 0.05), "quadrant"] <- 3
columbus@data[(mp$x <= 0 & mp$wx >= 0) & (columbus$lmp <= 0.05), "quadrant"] <- 4
columbus@data[(columbus$lmp > 0.05), "quadrant"] <- 5

#Clustering plot 

clusters <- tm_shape(columbus) + tm_fill(col = "quadrant", breaks = c(1, 2, 3, 4, 5, 6),
                             palette =  c("red", "blue", "lightpink", "skyblue2", "white"),
                             labels = c("High-High", "Low-Low", "High-Low", "Low-High", "Not significant"), title = "") +
  tm_legend(text.size = 1)  + tm_borders(alpha = .5) +
  tm_layout(frame = FALSE,  title = "Clusters",legend.title.size = 1.2, legend.text.size = 0.8 ,legend.outside = TRUE) 

clusters






rm(list = ls())