######################
######################
### FANOVA Figures ###
######################
######################


#clear environment and load libraries
rm(list = ls()); gc()
library(tidyverse)
library(INLA)
library(abind)
library(sf)
library(RColorBrewer)
library(scales)
library(maps)
library(plotly)
library(ggpubr)

#set colors
colormap <- colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32)


#helper functions
round_nearest_0.5 <- function(x) {
  round(x * 2) / 2
}

################################
### One-Way: Indy T2 and TKE ###
################################


################################################################################ --- T2

#load data
load('FANOVA_Data_Indy_T2_D02_City_Detrended_2.0.RData')
replicates = 2

# fanova_dat$citylocations[which(fanova_dat$citylocations == 0)]  = -1 #CHANGE THIS POSSIBLY

#load coordinates - d01 --- WE USE THIS FOR THE MESH
lat_d01 = readRDS('Data/d01lat_CROP')
lon_d01 = readRDS('Data/d01lon_CROP')
coords_d01 = expand.grid(lon_d01, lat_d01)
coord_sam = coords_d01
mesh1<-inla.mesh.create(loc=as.matrix(coord_sam))
SPDE = inla.spde2.matern(mesh1, alpha = 2)

#load FANOVA results
load('FANOVA_T2_Results_Indy_IID_D02Only_D01Mesh.RData')


#compile results
timepoints = 61
replicates = 2
actual_range = 1:dim(coords_d01)[1]
all_betas = matrix(NaN,
                   nrow = dim(coords_d01)[1],
                   ncol = timepoints)
all_lower = matrix(NaN,
                   nrow = dim(coords_d01)[1],
                   ncol = timepoints)
all_upper = matrix(NaN,
                   nrow = dim(coords_d01)[1],
                   ncol = timepoints)
for(t in 1:timepoints)
{
  betas = inla_output[[t]]$b1.field$mean
  lb = inla_output[[t]]$b1.field$`0.025quant`
  ub = inla_output[[t]]$b1.field$`0.975quant`
  estimate<-c()
  upper = c()
  lower = c()
  for(l in 1:(length(betas)/replicates)){
    estimate[l] <- mean(betas[seq(l, length(betas),by=length(betas)/replicates)])
    lower[l] = mean(lb[seq(l, length(ub),by=length(ub)/replicates)])
    upper[l] = mean(ub[seq(l, length(lb),by=length(lb)/replicates)])
  }
  
  all_betas[,t] = estimate[actual_range]
  all_lower[,t] = lower[actual_range]
  all_upper[,t] = upper[actual_range]
  
  print(t)
}


#get significant locations
wanted_time = 1:61
avg_beta = apply(all_betas[,wanted_time], 1, mean)
avg_lower = apply(all_lower[,wanted_time], 1, mean)
avg_upper = apply(all_upper[,wanted_time], 1, mean)
sig_locs = !c(avg_lower <= 0 & avg_upper >= 0)



#plot map of beta coeff
range = 1:dim(coords_d01)[1]
df <- cbind.data.frame(mesh1$loc[range,c(1,2)], beta = avg_beta)
colnames(df) = c('lon', 'lat', 'beta')

map <- df
states <- st_as_sf(map("state", plot = FALSE, fill = TRUE))
cities <- data.frame(city = c("Indianapolis"),
                     lat = c(39.791),
                     lon = c(-86.148))
cities <- st_as_sf(cities, coords = c("lon", "lat"), remove = FALSE,
                   crs = 4326, agr = "constant")


#generate beta plot
# bounds = round_nearest_0.5(max(abs(range(map$beta))) * c(-1,1))
bounds = 0.652 * c(-1,1)
beta_plot = ggplot() +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_line(color = "gray", size = 0.5),
    panel.grid.minor = element_blank(),
    panel.ontop = TRUE,
    legend.key.width = unit(1.25, "cm"),
    legend.key.height = unit(0.75, "cm"),
    legend.position = 'bottom'
  ) +
  geom_tile(map, mapping = aes(x = lon, y = lat, fill = beta)) +
  labs(fill = expression(hat(beta[city])), title = 'T2') +
  scale_fill_gradient2(low = "blue", mid = "white", high = muted("red"), midpoint = 0, limits = bounds) +  # Specifica i break per i tre range
  geom_sf(data = states, fill = "transparent", color = 'black', linewidth = 0.5) +
  coord_sf(xlim = c(-88.36627, -83.95374), ylim = c(38.00989, 41.48617), expand = TRUE) +
  xlab("") + ylab("") +
  theme(
    plot.background = element_blank(),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'),
    legend.title = element_text(vjust = 1, hjust = 1, size = 16),
    legend.text = element_text(size = 14)
  ) +
  theme(plot.background = element_blank()) +
  scale_x_continuous(
    breaks = c(-88, -87, -86, -85, -84),
    labels = c("88°", "87°", "86°", "85°", "84°")
  ) +
  scale_y_continuous(
    breaks = c(38, 39, 40, 41),
    labels = c("38°", "39°", "40°", "41°")
  )

# beta_plot




#plot map of sig. locations
range = 1:dim(coords_d01)[1]
df <- cbind.data.frame(mesh1$loc[range,c(1,2)], sig = sig_locs)
colnames(df) = c('lon', 'lat', 'sig')

map_sig <- df
sig_plot = ggplot() +
  geom_tile(data = map_sig,
            mapping = aes(x = lon, y = lat, fill = sig)) +
  labs(fill = "Significant?", title = 'T2: Significance', x='', y='') + 
  geom_sf(data = states, fill = "transparent", color = 'black', linewidth = 0.5) +
  coord_sf(xlim = c(-88.36627, -83.95374), ylim = c(38.00989, 41.48617), expand = TRUE) +
  theme(
    plot.background = element_blank(),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    plot.title = element_text(hjust = 0.5, face = 'bold', size = 16)
  ) +
  theme(plot.background = element_blank()) +
  scale_x_continuous(
    breaks = c(-88, -87, -86, -85, -84),
    labels = c("88°", "87°", "86°", "85°", "84°")
  ) +
  scale_y_continuous(
    breaks = c(38, 39, 40, 41),
    labels = c("38°", "39°", "40°", "41°")
  ) +
  scale_fill_manual(
    labels = c("No", "Yes"), values = c("grey80", "red")
  ) +
  theme(
    panel.background = element_blank(),
    axis.text = element_text(size = 20),
    panel.grid.major = element_line(color = "gray", size = 0.5),
    panel.grid.minor = element_blank(),
    panel.ontop = TRUE,
    legend.key.width = unit(0.5, "cm"),
    legend.key.height = unit(0.5, "cm"),
    legend.position = 'bottom',
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.key = element_rect(fill = NA)
  ) +
  guides(color = guide_legend(override.aes = list(size = 10),
                              label.position = "bottom"))

# sig_plot






################################################################ --- TKE






#load data
load('FANOVA_Data_Indy_TKE_D02_City_Detrended.RData')
replicates = 2


#load FANOVA results
load('FANOVA_TKE_Results_Indy_IID_D02Only_D01Mesh.RData')

#compile results
timepoints = 61
replicates = 2
actual_range = 1:dim(coords_d01)[1]
all_betas = matrix(NaN,
                   nrow = dim(coords_d01)[1],
                   ncol = timepoints)
all_lower = matrix(NaN,
                   nrow = dim(coords_d01)[1],
                   ncol = timepoints)
all_upper = matrix(NaN,
                   nrow = dim(coords_d01)[1],
                   ncol = timepoints)
for(t in 1:timepoints)
{
  full_ID = inla_output[[t]]$b1.field$ID
  betas = inla_output[[t]]$b1.field$mean
  lb = inla_output[[t]]$b1.field$`0.025quant`
  ub = inla_output[[t]]$b1.field$`0.975quant`
  estimate<-c()
  upper = c()
  lower = c()
  for(l in 1:(length(betas)/replicates)){
    estimate[l] <- mean(betas[seq(l, length(betas),by=length(betas)/replicates)])
    lower[l] = mean(lb[seq(l, length(ub),by=length(ub)/replicates)])
    upper[l] = mean(ub[seq(l, length(lb),by=length(lb)/replicates)])
  }
  
  all_betas[,t] = estimate[actual_range]
  all_lower[,t] = lower[actual_range]
  all_upper[,t] = upper[actual_range]
  
  print(t)
}



#get significant locations
wanted_time = 1:61
avg_beta = apply(all_betas[,wanted_time], 1, mean)
avg_lower = apply(all_lower[,wanted_time], 1, mean)
avg_upper = apply(all_upper[,wanted_time], 1, mean)
sig_locs = !c(avg_lower <= 0 & avg_upper >= 0)


#plot map of beta coeff
range = 1:dim(coords_d01)[1]
df <- cbind.data.frame(mesh1$loc[range,c(1,2)], beta = avg_beta)
colnames(df) = c('lon', 'lat', 'beta')

map_tke <- df
states <- st_as_sf(map("state", plot = FALSE, fill = TRUE))
cities <- data.frame(city = c("Indianapolis"),
                     lat = c(39.791),
                     lon = c(-86.148))
cities <- st_as_sf(cities, coords = c("lon", "lat"), remove = FALSE,
                   crs = 4326, agr = "constant")

colormap <- colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32)

bounds = max(abs(range(map_tke$beta))) * c(-1,1)
beta_plot_tke = ggplot() +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_line(color = "gray", size = 0.5),
    panel.grid.minor = element_blank(),
    panel.ontop = TRUE,
    legend.key.width = unit(1.25, "cm"),
    legend.key.height = unit(0.75, "cm"),
    legend.position = 'bottom'
  ) +
  geom_tile(map_tke, mapping = aes(x = lon, y = lat, fill = beta)) +
  labs(fill = expression(hat(beta[city])), title = 'TKE') +
  scale_fill_gradient2(low = "blue", mid = "white", high = muted("red"), midpoint = 0, limits = bounds) +  # Specifica i break per i tre range
  geom_sf(data = states, fill = "transparent", color = 'black', linewidth = 0.5) +
  coord_sf(xlim = c(-88.36627, -83.95374), ylim = c(38.00989, 41.48617), expand = TRUE) +
  xlab("") + ylab("") +
  theme(
    plot.background = element_blank(),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'),
    legend.title = element_text(vjust = 1, hjust = 1, size = 16),
    legend.text = element_text(size = 14)
  ) +
  theme(plot.background = element_blank()) +
  scale_x_continuous(
    breaks = c(-88, -87, -86, -85, -84),
    labels = c("88°", "87°", "86°", "85°", "84°")
  ) +
  scale_y_continuous(
    breaks = c(38, 39, 40, 41),
    labels = c("38°", "39°", "40°", "41°")
  )

# beta_plot_tke


#plot map of sig. locations
range = 1:dim(coords_d01)[1]
df <- cbind.data.frame(mesh1$loc[range,c(1,2)], sig = sig_locs)
colnames(df) = c('lon', 'lat', 'sig')

map_sig_tke <- df
sig_plot_tke = ggplot() +
  geom_tile(data = map_sig_tke,
            mapping = aes(x = lon, y = lat, fill = sig)) +
  labs(fill = "Significant?", title = 'TKE: Significance', x='', y='') + 
  geom_sf(data = states, fill = "transparent", color = 'black', linewidth = 0.5) +
  coord_sf(xlim = c(-88.36627, -83.95374), ylim = c(38.00989, 41.48617), expand = TRUE) +
  theme(
    plot.background = element_blank(),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    plot.title = element_text(hjust = 0.5, face = 'bold', size = 16)
  ) +
  theme(plot.background = element_blank()) +
  scale_x_continuous(
    breaks = c(-88, -87, -86, -85, -84),
    labels = c("88°", "87°", "86°", "85°", "84°")
  ) +
  scale_y_continuous(
    breaks = c(38, 39, 40, 41),
    labels = c("38°", "39°", "40°", "41°")
  ) +
  scale_fill_manual(
    labels = c("No", "Yes"), values = c("grey80", "red")
  ) +
  theme(
    panel.background = element_blank(),
    axis.text = element_text(size = 20),
    panel.grid.major = element_line(color = "gray", size = 0.5),
    panel.grid.minor = element_blank(),
    panel.ontop = TRUE,
    legend.key.width = unit(0.5, "cm"),
    legend.key.height = unit(0.5, "cm"),
    legend.position = 'bottom',
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.key = element_rect(fill = NA)
  ) +
  guides(color = guide_legend(override.aes = list(size = 10),
                              label.position = "bottom"))

# sig_plot_tke





beta_plots = ggarrange(beta_plot, beta_plot_tke, nrow = 2, ncol = 1,
                       labels = c('A', 'C'), common.legend = TRUE,
                       legend = 'bottom')
# beta_plots


sig_plots = ggarrange(sig_plot, sig_plot_tke, nrow = 2, ncol = 1,
                      labels = c('B', 'D'), common.legend = TRUE,
                      legend = 'bottom')

# sig_plots


all_plots = ggarrange(beta_plots, sig_plots, nrow = 1, ncol = 2)
all_plots




png(file = "FIGURES/F2-FANOVA_Results_Indy.png", width = 7 * 300, height = 9 * 300, res = 300)
all_plots
dev.off(); dev.off()



################################
### One-Way: NYC T2 and TKE ###
################################


################################################################################ --- T2

#load D03 grids
lat_d03 = readRDS('NYC/Data/Lat_grid_d03_NY')
lon_d03 = readRDS('NYC/Data/Lon_grid_d03_NY')
min_lat = min(lat_d03); max_lat = max(lat_d03)
min_lon = min(lon_d03); max_lon = max(lon_d03)

#load D02 grids
lat_d02 = readRDS('NYC/Data/Lat_grid_d01_NY')
lon_d02 = readRDS('NYC/Data/Lon_grid_d01_NY')

#check which d02 points are in D03
filter_lat = (lat_d02 >= min_lat & lat_d02 <= max_lat)
filter_lon = (lon_d02 >= min_lon & lon_d02 <= max_lon)
full_filter = filter_lat & filter_lon

#restrict the D02 domain
lat_restrict = which(apply(full_filter, 2, sum) >= 1)
lon_restrict = which(apply(full_filter, 1, sum) >= 1)


lat_mat = lat_d02[lon_restrict, lat_restrict]
lon_mat = lon_d02[lon_restrict, lat_restrict]

#create mesh
coords_d01 = cbind(c(lon_mat), c(lat_mat))
coord_sam = coords_d01
mesh1<-inla.mesh.create(loc=as.matrix(coord_sam))


#load data
load('NYC/FANOVA_Data_NYC_T2_D02_City_Detrended.RData')
replicates = 2



#load FANOVA results
load('NYC/FANOVA_T2_Results_NYC_IID_D02Only_D01Mesh.RData')



#compile results
timepoints = 61
replicates = 2
actual_range = 1:dim(coords_d01)[1]
all_betas = matrix(NaN,
                   nrow = dim(coords_d01)[1],
                   ncol = timepoints)
all_lower = matrix(NaN,
                   nrow = dim(coords_d01)[1],
                   ncol = timepoints)
all_upper = matrix(NaN,
                   nrow = dim(coords_d01)[1],
                   ncol = timepoints)
for(t in 1:timepoints)
{
  betas = inla_output[[t]]$b1.field$mean
  lb = inla_output[[t]]$b1.field$`0.025quant`
  ub = inla_output[[t]]$b1.field$`0.975quant`
  estimate<-c()
  upper = c()
  lower = c()
  for(l in 1:(length(betas)/replicates)){
    estimate[l] <- mean(betas[seq(l, length(betas),by=length(betas)/replicates)])
    lower[l] = mean(lb[seq(l, length(ub),by=length(ub)/replicates)])
    upper[l] = mean(ub[seq(l, length(lb),by=length(lb)/replicates)])
  }
  
  all_betas[,t] = estimate[actual_range]
  all_lower[,t] = lower[actual_range]
  all_upper[,t] = upper[actual_range]
  
  print(t)
}


#get significant locations
wanted_time = 1:61
avg_beta = apply(all_betas[,wanted_time], 1, mean)
avg_lower = apply(all_lower[,wanted_time], 1, mean)
avg_upper = apply(all_upper[,wanted_time], 1, mean)
sig_locs = !c(avg_lower <= 0 & avg_upper >= 0)



#plot map of beta coeff
range = 1:dim(coords_d01)[1]
df <- cbind.data.frame(mesh1$loc[range,c(1,2)], beta = avg_beta)
colnames(df) = c('lon', 'lat', 'beta')

map <- df
states <- st_as_sf(map("state", plot = FALSE, fill = TRUE))
cities <- data.frame(city = c("New York City"),
                     lat = c(40.7128),
                     lon = c(-74.0060))
cities <- st_as_sf(cities, coords = c("lon", "lat"), remove = FALSE,
                   crs = 4326, agr = "constant")

#generate beta plot
# bounds = round_nearest_0.5(max(abs(range(map$beta))) * c(-1,1))
bounds = 0.828 * c(-1,1)
beta_plot = ggplot() +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_line(color = "gray", size = 0.5),
    panel.grid.minor = element_blank(),
    panel.ontop = TRUE,
    legend.key.width = unit(1.25, "cm"),
    legend.key.height = unit(0.75, "cm"),
    legend.position = 'bottom'
  ) +
  geom_point(map, mapping = aes(x = lon, y = lat, color = beta), shape = 15, size = 10) +
  labs(color = expression(hat(beta[city])), title = 'T2') +
  scale_color_gradient2(low = "blue", mid = "white", high = muted("red"), midpoint = 0, limits = bounds) +  # Specifica i break per i tre range
  geom_sf(data = states, fill = "transparent", color = 'black', linewidth = 0.5) +
  # coord_sf(xlim = c(-76.21, -71.28), ylim = c(39.05, 42.6), expand = TRUE) +
  coord_sf(xlim = c(-76, -71.8), ylim = c(39.25, 42.2), expand = TRUE) +
  xlab("") + ylab("") +
  theme(
    plot.background = element_blank(),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'),
    legend.title = element_text(vjust = 1, hjust = 1, size = 16),
    legend.text = element_text(size = 14)
  ) +
  theme(plot.background = element_blank()) +
  scale_x_continuous(
    breaks = c(-76, -75, -74, -73, -72),
    labels = c("76°", "75°", "74°", "73°", "72°")
  ) +
  scale_y_continuous(
    breaks = c(40, 41, 42),
    labels = c("40°", "41°", "42°")
  )

# beta_plot




#plot map of sig. locations
range = 1:dim(coords_d01)[1]
df <- cbind.data.frame(mesh1$loc[range,c(1,2)], sig = sig_locs)
colnames(df) = c('lon', 'lat', 'sig')

map_sig <- df
sig_plot = ggplot() +
  geom_point(data = map_sig,
            mapping = aes(x = lon, y = lat, color = sig),
            shape = 15, size = 10) +
  labs(color = "Significant?", title = 'T2: Significance', x='', y='') + 
  geom_sf(data = states, fill = "transparent", color = 'black', linewidth = 0.5) +
  # coord_sf(xlim = c(-76.21, -71.28), ylim = c(39.05, 42.6), expand = TRUE) +
  coord_sf(xlim = c(-76, -71.8), ylim = c(39.25, 42.2), expand = TRUE) +
  theme(
    plot.background = element_blank(),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    plot.title = element_text(hjust = 0.5, face = 'bold', size = 16)
  ) +
  theme(plot.background = element_blank()) +
  scale_x_continuous(
    breaks = c(-76, -75, -74, -73, -72),
    labels = c("76°", "75°", "74°", "73°", "72°")
  ) +
  scale_y_continuous(
    breaks = c(40, 41, 42),
    labels = c("40°", "41°", "42°")
  ) +
  scale_color_manual(
    labels = c("No", "Yes"), values = c("grey80", "red")
  ) +
  theme(
    panel.background = element_blank(),
    axis.text = element_text(size = 20),
    panel.grid.major = element_line(color = "gray", size = 0.5),
    panel.grid.minor = element_blank(),
    panel.ontop = TRUE,
    legend.key.width = unit(0.5, "cm"),
    legend.key.height = unit(0.5, "cm"),
    legend.position = 'bottom',
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.key = element_rect(fill = NA)
  )

# sig_plot






################################################################ --- TKE






#load data
load('NYC/FANOVA_Data_NYC_TKE_D02_City_Detrended.RData')
replicates = 2


#load FANOVA results
load('NYC/FANOVA_TKE_Results_NYC_IID_D02Only_D01Mesh.RData')

#compile results
timepoints = 61
replicates = 2
actual_range = 1:dim(coords_d01)[1]
all_betas = matrix(NaN,
                   nrow = dim(coords_d01)[1],
                   ncol = timepoints)
all_lower = matrix(NaN,
                   nrow = dim(coords_d01)[1],
                   ncol = timepoints)
all_upper = matrix(NaN,
                   nrow = dim(coords_d01)[1],
                   ncol = timepoints)
for(t in 1:timepoints)
{
  full_ID = inla_output[[t]]$b1.field$ID
  betas = inla_output[[t]]$b1.field$mean
  lb = inla_output[[t]]$b1.field$`0.025quant`
  ub = inla_output[[t]]$b1.field$`0.975quant`
  estimate<-c()
  upper = c()
  lower = c()
  for(l in 1:(length(betas)/replicates)){
    estimate[l] <- mean(betas[seq(l, length(betas),by=length(betas)/replicates)])
    lower[l] = mean(lb[seq(l, length(ub),by=length(ub)/replicates)])
    upper[l] = mean(ub[seq(l, length(lb),by=length(lb)/replicates)])
  }
  
  all_betas[,t] = estimate[actual_range]
  all_lower[,t] = lower[actual_range]
  all_upper[,t] = upper[actual_range]
  
  print(t)
}



#get significant locations
wanted_time = 1:61
avg_beta = apply(all_betas[,wanted_time], 1, mean)
avg_lower = apply(all_lower[,wanted_time], 1, mean)
avg_upper = apply(all_upper[,wanted_time], 1, mean)
sig_locs = !c(avg_lower <= 0 & avg_upper >= 0)


#plot map of beta coeff
range = 1:dim(coords_d01)[1]
df <- cbind.data.frame(mesh1$loc[range,c(1,2)], beta = avg_beta)
colnames(df) = c('lon', 'lat', 'beta')

map_tke <- df
states <- st_as_sf(map("state", plot = FALSE, fill = TRUE))
cities <- data.frame(city = c("New York City"),
                     lat = c(40.7128),
                     lon = c(-74.0060))
cities <- st_as_sf(cities, coords = c("lon", "lat"), remove = FALSE,
                   crs = 4326, agr = "constant")


# bounds = max(abs(range(map_tke$beta))) * c(-1,1)
bounds = 0.828 * c(-1,1)
bounds = 
beta_plot_tke = ggplot() +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_line(color = "gray", size = 0.5),
    panel.grid.minor = element_blank(),
    panel.ontop = TRUE,
    legend.key.width = unit(1.25, "cm"),
    legend.key.height = unit(0.75, "cm"),
    legend.position = 'bottom'
  ) +
  geom_point(map_tke, mapping = aes(x = lon, y = lat, color = beta), shape = 15, size = 10) +
  labs(color = expression(hat(beta[city])), title = 'TKE') +
  scale_color_gradient2(low = "blue", mid = "white", high = muted("red"), midpoint = 0, limits = bounds) +  # Specifica i break per i tre range
  geom_sf(data = states, fill = "transparent", color = 'black', linewidth = 0.5) +
  coord_sf(xlim = c(-76, -71.8), ylim = c(39.25, 42.2), expand = TRUE) +
  xlab("") + ylab("") +
  theme(
    plot.background = element_blank(),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'),
    legend.title = element_text(vjust = 1, hjust = 1, size = 16),
    legend.text = element_text(size = 14)
  ) +
  theme(plot.background = element_blank()) +
  scale_x_continuous(
    breaks = c(-76, -75, -74, -73, -72),
    labels = c("76°", "75°", "74°", "73°", "72°")
  ) +
  scale_y_continuous(
    breaks = c(40, 41, 42),
    labels = c("40°", "41°", "42°")
  )

# beta_plot_tke


#plot map of sig. locations
range = 1:dim(coords_d01)[1]
df <- cbind.data.frame(mesh1$loc[range,c(1,2)], sig = sig_locs)
colnames(df) = c('lon', 'lat', 'sig')

map_sig_tke <- df
sig_plot_tke = ggplot() +
  geom_point(data = map_sig_tke,
             mapping = aes(x = lon, y = lat, color = sig),
             shape = 15, size = 10) +
  labs(color = "Significant?", title = 'TKE: Significance', x='', y='') + 
  geom_sf(data = states, fill = "transparent", color = 'black', linewidth = 0.5) +
  coord_sf(xlim = c(-76, -71.8), ylim = c(39.25, 42.2), expand = TRUE) +
  theme(
    plot.background = element_blank(),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    plot.title = element_text(hjust = 0.5, face = 'bold', size = 16)
  ) +
  theme(plot.background = element_blank()) +
  scale_x_continuous(
    breaks = c(-76, -75, -74, -73, -72),
    labels = c("76°", "75°", "74°", "73°", "72°")
  ) +
  scale_y_continuous(
    breaks = c(40, 41, 42),
    labels = c("40°", "41°", "42°")
  )+
  scale_color_manual(
    labels = c("No", "Yes"), values = c("grey80", "red")
  ) +
  theme(
    panel.background = element_blank(),
    axis.text = element_text(size = 20),
    panel.grid.major = element_line(color = "gray", size = 0.5),
    panel.grid.minor = element_blank(),
    panel.ontop = TRUE,
    legend.key.width = unit(0.5, "cm"),
    legend.key.height = unit(0.5, "cm"),
    legend.position = 'bottom',
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.key = element_rect(fill = NA)
  )

# sig_plot_tke





beta_plots = ggarrange(beta_plot, beta_plot_tke, nrow = 2, ncol = 1,
                       labels = c('A', 'C'), common.legend = TRUE,
                       legend = 'bottom')
# beta_plots


sig_plots = ggarrange(sig_plot, sig_plot_tke, nrow = 2, ncol = 1,
                      labels = c('B', 'D'), common.legend = TRUE,
                      legend = 'bottom')

# sig_plots


all_plots = ggarrange(beta_plots, sig_plots, nrow = 1, ncol = 2)
all_plots




png(file = "FIGURES/F3-FANOVA_Results_NYC.png", width = 7 * 300, height = 9 * 300, res = 300)
all_plots
dev.off(); dev.off()







##################################
### Multi-Way FANOVA: PRCP DFW ###
##################################




#Get colors and dallar coordinates
states <- st_as_sf(map("state", plot = FALSE, fill = TRUE))
cities <- data.frame(city = c("Dallas"),
                     lat = c(32.7767),
                     lon = c(-96.7970))
cities <- st_as_sf(cities, coords = c("lon", "lat"), remove = FALSE,
                   crs = 4326, agr = "constant")



#load data
load('Dallas/PRCP_CF_City_Emissions/FANOVA_DFW_March_PostCityALLTIMEPRCP_D02.RData')
fanova_dat = fanova_dat %>%
  filter(lon < -96, lon > -97.5, time > 12)
fanova_dat$time = fanova_dat$time - 12
locations = unique(cbind(fanova_dat$lon, fanova_dat$lat))

# 
# temp = fanova_dat %>%
#   arrange(desc(precipitation)) %>%
#   filter(time == 6, precipitation >= 10)
# large_locs = unique(cbind(temp$lon, temp$lat))
# 
# plot(locations, pch = 15, cex = 3)
# points(large_locs, pch = 15, cex = 3, col = 'red')

max_order = rev(order(fanova_dat$precipitation))
fanova_dat[max_order,]

#load results
method = 'IID'
size = 'D01'
fanova_results = readRDS(paste0('Dallas/PRCP_CF_City_Emissions/FANOVA Results/FANOVA_Results_DFW_March_PostCityLASTHOURPRCP_IID_', size, 'Mesh.rds'))
wanted_locs = unique(fanova_results$b1.field$ID)

#compile betas and significance
beta1 = fanova_results$b1.field %>%
  group_by(ID) %>%
  summarize(avg_beta = mean(mean),
            lower = mean(`0.025quant`),
            upper = mean(`0.975quant`))
significant1 = !(beta1$lower <= 0 & beta1$upper >=0)

beta2 = fanova_results$b2.field %>%
  group_by(ID) %>%
  summarize(avg_beta = mean(mean),
            lower = mean(`0.025quant`),
            upper = mean(`0.975quant`))
significant2 = !(beta2$lower <= 0 & beta2$upper >=0)



#load coordinate grids and filter down to D03 space
lat_d02 = c(readRDS(paste0('Dallas/PRCP_CF_City_Emissions/lat_grid_', size)))
lon_d02 = c(readRDS(paste0('Dallas/PRCP_CF_City_Emissions/lon_grid_', size)))
d02_grid = cbind(lon_d02, lat_d02)
total_d02 = dim(d02_grid)[1]
filter1 = d02_grid[,1] >= -98.5 & d02_grid[,1] <= -95
filter2 = d02_grid[,2] >= 31.5 & d02_grid[,2] <= 34
d02_filter = filter1 & filter2
d02_grid = d02_grid[d02_filter,]

#plot locations of interest

filter = wanted_locs <= dim(d02_grid)[1]
# plot(d02_grid, pch = 19)
# points(d02_grid[wanted_locs[filter], ], col = 'red', pch = 19)
new_grid = d02_grid[wanted_locs[filter], ]





#combine to make beta1 dataframe
beta1_dat <- cbind.data.frame(new_grid, beta1$avg_beta[filter])
colnames(beta1_dat) = c('lon', 'lat', 'beta1')
bounds = max(c(abs(min(beta1$avg_beta, na.rm = TRUE)),
               abs(max(beta1$avg_beta, na.rm = TRUE))))



beta1_plot = ggplot() +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_line(color = "gray", size = 0.5),
    panel.grid.minor = element_blank(),
    panel.ontop = TRUE,
    legend.key.width = unit(1.25, "cm"),
    legend.key.height = unit(0.75, "cm"),
    legend.position = 'bottom'
  ) +
  geom_point(data = beta1_dat, aes(x = lon,
                                   y = lat,
                                   color = beta1),
             shape = 15,
             size = 8) +
  labs(color = expression(hat(beta)[City]), x = '', y = '', title = 'City') +
  scale_color_gradient2(low = "blue",
                        mid = 'white',
                        high = muted("red"),
                        midpoint = 0,
                        limits = c(-bounds, bounds)) +
  geom_sf(data = states,
          fill = "transparent",
          color = 'black',
          linewidth = 0.5) +
  geom_text(data = cities, aes(x = lon, y = lat, label = city), size = 5, col = "black", fontface = "bold") +
  theme(
    plot.background = element_blank(),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    plot.title = element_text(hjust = 0.5, size = 18, face = 'bold'),
    legend.title = element_text(vjust = 1, hjust = 1, size = 16),
    legend.text = element_text(size = 14)
  ) +
  theme(plot.background = element_blank()) +
  coord_sf(xlim = c(-98.5, -95), 
           ylim = c(31.5, 34), 
           expand = TRUE)+ 
  scale_x_continuous(
    breaks = c(-98, -97, -96, -95),
    labels = c("98°", "97°", "96°", "95°")
  ) +
  scale_y_continuous(
    breaks = c(31.5, 32, 32.5, 33, 33.5, 34),
    labels = c("31.5°", "32°", "32.5°", "33°", "33.5°", "34°")
  ) 


# beta1_plot




##############################plot results for BETA1 Significance (City)

# par(mfrow = c(1,2))
# plot(beta1$lower)
# plot(beta1$upper)



beta1_dat_sig <- cbind.data.frame(new_grid, significant1[filter])
colnames(beta1_dat_sig) = c('lon', 'lat', 'significant1')

beta1_plot_sig = ggplot() +
  theme(
    plot.background = element_blank(),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    plot.title = element_text(hjust = 0.5, face = 'bold', size = 18)
  ) +
  theme(plot.background = element_blank()) +
  geom_point(data = beta1_dat_sig, aes(x = lon,
                                       y = lat,
                                       color = significant1),
             shape = 15, size = 8) +
  labs(color = expression(hat(beta)[City] * ' Significant?'), x = '', y = '', title = 'City: Significance') +
  geom_sf(data = states,
          fill = "transparent",
          color = 'black',
          linewidth = 0.5) +
  geom_text(data = cities, aes(x = lon, y = lat, label = city), size = 5, col = "black", fontface = "bold") +
  scale_color_manual(
    labels = c("No", "Yes"), values = c("grey80", "red")
  ) +
  theme(
    panel.background = element_blank(),
    axis.text = element_text(size = 20),
    panel.grid.major = element_line(color = "gray", size = 0.5),
    panel.grid.minor = element_blank(),
    panel.ontop = TRUE,
    legend.key.width = unit(0.5, "cm"),
    legend.key.height = unit(0.5, "cm"),
    legend.position = 'bottom',
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.key = element_rect(fill = NA)
  ) +
  guides(color = guide_legend(override.aes = list(size = 10),
                              label.position = "right")) +
  coord_sf(xlim = c(-98.5, -95), 
           ylim = c(31.5, 34), 
           expand = TRUE) + 
  scale_x_continuous(
    breaks = c(-98, -97, -96, -95),
    labels = c("98°", "97°", "96°", "95°")
  ) +
  scale_y_continuous(
    breaks = c(31.5, 32, 32.5, 33, 33.5, 34),
    labels = c("31.5°", "32°", "32.5°", "33°", "33.5°", "34°")
  ) 



# beta1_plot_sig







###############################plot results for BETA2 (Emissions)



#combine to make beta1 dataframe
beta2_dat <- cbind.data.frame(new_grid, beta2$avg_beta[filter])
colnames(beta2_dat) = c('lon', 'lat', 'beta2')
bounds2 = max(c(abs(min(beta2$avg_beta, na.rm = TRUE)),
                abs(max(beta2$avg_beta, na.rm = TRUE))))



beta2_plot = ggplot() +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_line(color = "gray", size = 0.5),
    panel.grid.minor = element_blank(),
    panel.ontop = TRUE,
    legend.key.width = unit(1.25, "cm"),
    legend.key.height = unit(0.75, "cm"),
    legend.position = 'bottom'
  ) +
  geom_point(data = beta2_dat, aes(x = lon,
                                   y = lat,
                                   color = beta2),
             shape = 15,
             size = 8) +
  labs(color = expression(hat(beta)[Emiss]), x = '', y = '', title = 'Emissions') +
  scale_color_gradient2(low = "blue",
                        mid = 'white',
                        high = muted("red"),
                        midpoint = 0,
                        limits = c(-bounds2, bounds2)) +
  geom_sf(data = states,
          fill = "transparent",
          color = 'black',
          linewidth = 0.5) +
  geom_text(data = cities, aes(x = lon, y = lat, label = city), size = 5, col = "black", fontface = "bold") +
  theme(
    plot.background = element_blank(),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    plot.title = element_text(hjust = 0.5, size = 18, face = 'bold'),
    legend.title = element_text(vjust = 1, hjust = 1, size = 16),
    legend.text = element_text(size = 14)
  ) +
  theme(plot.background = element_blank()) +
  coord_sf(xlim = c(-98.5, -95), 
           ylim = c(31.5, 34), 
           expand = TRUE) + 
  scale_x_continuous(
    breaks = c(-98, -97, -96, -95),
    labels = c("98°", "97°", "96°", "95°")
  ) +
  scale_y_continuous(
    breaks = c(31.5, 32, 32.5, 33, 33.5, 34),
    labels = c("31.5°", "32°", "32.5°", "33°", "33.5°", "34°")
  ) 


# beta2_plot




###############################plot results for BETA2 Significance (Emissions)

# par(mfrow = c(1,2))
# plot(beta2$lower)
# plot(beta2$upper)



beta2_dat_sig <- cbind.data.frame(new_grid, significant2[filter])
colnames(beta2_dat_sig) = c('lon', 'lat', 'significant2')

beta2_plot_sig = ggplot() +
  theme(
    plot.background = element_blank(),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    plot.title = element_text(hjust = 0.5, face = 'bold', size = 18)
  ) +
  theme(plot.background = element_blank()) +
  geom_point(data = beta2_dat_sig, aes(x = lon,
                                       y = lat,
                                       color = significant2),
             shape = 15, size = 8) +
  labs(color = expression(hat(beta)[Emiss] * ' Significant?'), x = '', y = '', title = 'Emissions: Significance') +
  geom_sf(data = states,
          fill = "transparent",
          color = 'black',
          linewidth = 0.5) +
  geom_text(data = cities, aes(x = lon, y = lat, label = city), size = 5, col = "black", fontface = "bold") +
  scale_color_manual(
    labels = c("No", "Yes"), values = c("grey80", "red")
  ) +
  theme(
    panel.background = element_blank(),
    axis.text = element_text(size = 20),
    panel.grid.major = element_line(color = "gray", size = 0.5),
    panel.grid.minor = element_blank(),
    panel.ontop = TRUE,
    legend.key.width = unit(0.5, "cm"),
    legend.key.height = unit(0.5, "cm"),
    legend.position = 'bottom',
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.key = element_rect(fill = NA)
  ) +
  guides(color = guide_legend(override.aes = list(size = 10),
                              label.position = "right")) +
  coord_sf(xlim = c(-98.5, -95), 
           ylim = c(31.5, 34), 
           expand = TRUE) + 
  scale_x_continuous(
    breaks = c(-98, -97, -96, -95),
    labels = c("98°", "97°", "96°", "95°")
  ) +
  scale_y_continuous(
    breaks = c(31.5, 32, 32.5, 33, 33.5, 34),
    labels = c("31.5°", "32°", "32.5°", "33°", "33.5°", "34°")
  ) 



# beta2_plot_sig




beta_plots = ggarrange(beta1_plot, beta2_plot,
                       nrow = 2, ncol = 1, 
                       labels = c('A', 'C'))


sig_plots = ggarrange(beta1_plot_sig, beta2_plot_sig,
                      nrow = 2, ncol = 1, 
                      labels = c('B', 'D'))


all_plots = ggarrange(beta_plots, sig_plots, 
                      nrow = 1, ncol = 2)



png(file = "FIGURES/F3-FANOVA_Results_DFW.png", width = 10 * 300, height = 10 * 300, res = 300)
all_plots
dev.off(); dev.off()






##############################
### PRCP DFW - Simulations ###
##############################


color_range = c('snow1', 'lightblue1', 'cyan3', 'blue', 'navy', 'midnightblue')

#Get colors and dallar coordinates
states <- st_as_sf(map("state", plot = FALSE, fill = TRUE))
cities <- data.frame(city = c("Dallas"),
                     lat = c(32.7767),
                     lon = c(-96.7970))
cities <- st_as_sf(cities, coords = c("lon", "lat"), remove = FALSE,
                   crs = 4326, agr = "constant")



#load data
load('Dallas/PRCP_CF_City_Emissions/FANOVA_DFW_March_PostCityALLTIMEPRCP_D02.RData')
fanova_dat = fanova_dat %>%
  filter(lon < -96, lon > -97.5, time > 12)
fanova_dat$time = fanova_dat$time - 12
locations = unique(cbind(fanova_dat$lon, fanova_dat$lat))
colnames(locations) = c('lon', 'lat')


#load results
method = 'IID'
size = 'D01'
fanova_results = readRDS(paste0('Dallas/PRCP_CF_City_Emissions/FANOVA Results/FANOVA_Results_DFW_March_PostCityLASTHOURPRCP_IID_', size, 'Mesh.rds'))
wanted_locs = unique(fanova_results$b1.field$ID)




#load coordinate grids and filter down to D03 space
lat_d02 = c(readRDS(paste0('Dallas/PRCP_CF_City_Emissions/lat_grid_', size)))
lon_d02 = c(readRDS(paste0('Dallas/PRCP_CF_City_Emissions/lon_grid_', size)))
d02_grid = cbind(lon_d02, lat_d02)
total_d02 = dim(d02_grid)[1]
filter1 = d02_grid[,1] >= -98.5 & d02_grid[,1] <= -95
filter2 = d02_grid[,2] >= 31.5 & d02_grid[,2] <= 34
d02_filter = filter1 & filter2
d02_grid = d02_grid[d02_filter,]

#plot locations of interest

filter = wanted_locs <= dim(d02_grid)[1]
# plot(d02_grid, pch = 19)
# points(d02_grid[wanted_locs[filter], ], col = 'red', pch = 19)
new_grid = d02_grid[wanted_locs[filter], ]

#partition fanova data to extract the 4 scenarios
cityemix = fanova_dat %>%
  filter(repl == 4) %>%
  group_by(location) %>%
  summarize(tot_prcp = sum(precipitation)) %>%
  cbind(locations)
nocityemix = fanova_dat %>%
  filter(repl == 3) %>%
  group_by(location) %>%
  summarize(tot_prcp = sum(precipitation)) %>%
  cbind(locations)
citynoemix = fanova_dat %>%
  filter(repl == 2) %>%
  group_by(location) %>%
  summarize(tot_prcp = sum(precipitation)) %>%
  cbind(locations)
nocitynoemix = fanova_dat %>%
  filter(repl == 1) %>%
  group_by(location) %>%
  summarize(tot_prcp = sum(precipitation)) %>%
  cbind(locations)


all_dat = rbind(cityemix, citynoemix, nocityemix, nocitynoemix)
max_total = max(all_dat$tot_prcp)
min_total = min(all_dat$tot_prcp)
bounds = c(min_total, ceiling(max_total))

#combine to make beta1 dataframe
plot1_dat <- cityemix[,-1]
colnames(plot1_dat) = c('beta1','lon', 'lat')



plot1 = ggplot() +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_line(color = "gray", size = 0.5),
    panel.grid.minor = element_blank(),
    panel.ontop = TRUE,
    legend.key.width = unit(1.25, "cm"),
    legend.key.height = unit(0.75, "cm"),
    legend.position = 'bottom'
  ) +
  geom_point(data = plot1_dat, aes(x = lon,
                                   y = lat,
                                   color = beta1),
             shape = 15,
             size = 8) +
  labs(color = 'mm', x = '', y = '', title = 'City + Emissions') +
  scale_color_gradientn(colors = color_range,
                       limit = bounds) +
  geom_sf(data = states,
          fill = "transparent",
          color = 'black',
          linewidth = 0.5) +
  geom_text(data = cities, aes(x = lon, y = lat, label = city), size = 5, col = "black", fontface = "bold") +
  theme(
    plot.background = element_blank(),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    plot.title = element_text(hjust = 0.5, size = 18, face = 'bold'),
    legend.title = element_text(vjust = 1, hjust = 1, size = 16),
    legend.text = element_text(size = 14)
  ) +
  theme(plot.background = element_blank()) +
  coord_sf(xlim = c(-98.5, -95), 
           ylim = c(31.5, 34), 
           expand = TRUE)+ 
  scale_x_continuous(
    breaks = c(-98, -97, -96, -95),
    labels = c("98°", "97°", "96°", "95°")
  ) +
  scale_y_continuous(
    breaks = c(31.5, 32, 32.5, 33, 33.5, 34),
    labels = c("31.5°", "32°", "32.5°", "33°", "33.5°", "34°")
  ) 


plot1




##############################plot results for BETA1 Significance (City)

#combine to make beta1 dataframe
plot2_dat <- nocityemix[,-1]
colnames(plot2_dat) = c('beta1','lon', 'lat')



plot2 = ggplot() +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_line(color = "gray", size = 0.5),
    panel.grid.minor = element_blank(),
    panel.ontop = TRUE,
    legend.key.width = unit(1.25, "cm"),
    legend.key.height = unit(0.75, "cm"),
    legend.position = 'bottom'
  ) +
  geom_point(data = plot2_dat, aes(x = lon,
                                   y = lat,
                                   color = beta1),
             shape = 15,
             size = 8) +
  labs(color = 'mm', x = '', y = '', title = 'No City + Emissions') +
  scale_color_gradientn(colors = color_range,
                        limit = bounds) +
  geom_sf(data = states,
          fill = "transparent",
          color = 'black',
          linewidth = 0.5) +
  geom_text(data = cities, aes(x = lon, y = lat, label = city), size = 5, col = "black", fontface = "bold") +
  theme(
    plot.background = element_blank(),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    plot.title = element_text(hjust = 0.5, size = 18, face = 'bold'),
    legend.title = element_text(vjust = 1, hjust = 1, size = 16),
    legend.text = element_text(size = 14)
  ) +
  theme(plot.background = element_blank()) +
  coord_sf(xlim = c(-98.5, -95), 
           ylim = c(31.5, 34), 
           expand = TRUE)+ 
  scale_x_continuous(
    breaks = c(-98, -97, -96, -95),
    labels = c("98°", "97°", "96°", "95°")
  ) +
  scale_y_continuous(
    breaks = c(31.5, 32, 32.5, 33, 33.5, 34),
    labels = c("31.5°", "32°", "32.5°", "33°", "33.5°", "34°")
  ) 


# plot2







###############################plot results for BETA2 (Emissions)



#combine to make beta1 dataframe
plot3_dat <- citynoemix[,-1]
colnames(plot3_dat) = c('beta1','lon', 'lat')



plot3 = ggplot() +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_line(color = "gray", size = 0.5),
    panel.grid.minor = element_blank(),
    panel.ontop = TRUE,
    legend.key.width = unit(1.25, "cm"),
    legend.key.height = unit(0.75, "cm"),
    legend.position = 'bottom'
  ) +
  geom_point(data = plot3_dat, aes(x = lon,
                                   y = lat,
                                   color = beta1),
             shape = 15,
             size = 8) +
  labs(color = 'mm', x = '', y = '', title = 'City + No Emissions') +
  scale_color_gradientn(colors = color_range,
                        limit = bounds) +
  geom_sf(data = states,
          fill = "transparent",
          color = 'black',
          linewidth = 0.5) +
  geom_text(data = cities, aes(x = lon, y = lat, label = city), size = 5, col = "black", fontface = "bold") +
  theme(
    plot.background = element_blank(),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    plot.title = element_text(hjust = 0.5, size = 18, face = 'bold'),
    legend.title = element_text(vjust = 1, hjust = 1, size = 16),
    legend.text = element_text(size = 14)
  ) +
  theme(plot.background = element_blank()) +
  coord_sf(xlim = c(-98.5, -95), 
           ylim = c(31.5, 34), 
           expand = TRUE)+ 
  scale_x_continuous(
    breaks = c(-98, -97, -96, -95),
    labels = c("98°", "97°", "96°", "95°")
  ) +
  scale_y_continuous(
    breaks = c(31.5, 32, 32.5, 33, 33.5, 34),
    labels = c("31.5°", "32°", "32.5°", "33°", "33.5°", "34°")
  ) 


# plot3




###############################plot results for BETA2 Significance (Emissions)

#combine to make beta1 dataframe
plot4_dat <- nocitynoemix[,-1]
colnames(plot4_dat) = c('beta1','lon', 'lat')



plot4 = ggplot() +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_line(color = "gray", size = 0.5),
    panel.grid.minor = element_blank(),
    panel.ontop = TRUE,
    legend.key.width = unit(1.25, "cm"),
    legend.key.height = unit(0.75, "cm"),
    legend.position = 'bottom'
  ) +
  geom_point(data = plot4_dat, aes(x = lon,
                                   y = lat,
                                   color = beta1),
             shape = 15,
             size = 8) +
  labs(color = 'mm', x = '', y = '', title = 'No City + No Emissions') +
  scale_color_gradientn(colors = color_range,
                        limit = bounds) +
  geom_sf(data = states,
          fill = "transparent",
          color = 'black',
          linewidth = 0.5) +
  geom_text(data = cities, aes(x = lon, y = lat, label = city), size = 5, col = "black", fontface = "bold") +
  theme(
    plot.background = element_blank(),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    plot.title = element_text(hjust = 0.5, size = 18, face = 'bold'),
    legend.title = element_text(vjust = 1, hjust = 1, size = 16),
    legend.text = element_text(size = 14)
  ) +
  theme(plot.background = element_blank()) +
  coord_sf(xlim = c(-98.5, -95), 
           ylim = c(31.5, 34), 
           expand = TRUE)+ 
  scale_x_continuous(
    breaks = c(-98, -97, -96, -95),
    labels = c("98°", "97°", "96°", "95°")
  ) +
  scale_y_continuous(
    breaks = c(31.5, 32, 32.5, 33, 33.5, 34),
    labels = c("31.5°", "32°", "32.5°", "33°", "33.5°", "34°")
  ) 


# plot4






all_plots = ggarrange(plot1, plot2, plot3, plot4, 
                      nrow = 2, ncol = 2,
                      common.legend = TRUE, legend = 'bottom',
                      labels = c('A', 'B', 'C', 'D'))



png(file = "FIGURES/F4-DFW_Total_PRCP.png", width = 10 * 300, height = 10 * 300, res = 300)
all_plots
dev.off(); dev.off()












