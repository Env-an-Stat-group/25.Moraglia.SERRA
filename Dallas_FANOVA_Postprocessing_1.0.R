##########################
##########################
### Dallas Data FANOVA ###
##########################
##########################

#clear environment
rm(list = ls()); gc()

#load libraries 
library(tidyverse)
library(Matrix)
library(abind)
library(INLA)



#load data
#fanova_dat = readRDS('FANOVA_Data_Dallas_Detrended_D03Only.rds')
load('FANOVA_Data_Dallas_Detrended_D03Only.RData')
replicates = max(fanova_dat$repl)


#load coordinates - d01 --- WE USE THIS FOR THE MESH
d02_grid = readRDS('dallas_D03_grid.rds')
total_d02 = dim(d02_grid)[1]
filter1 = d02_grid[,1] >= -99.95 & d02_grid[,1] <= -93.97
filter2 = d02_grid[,2] >= 30.35 & d02_grid[,2] <= 35.26
d02_filter = filter1 & filter2
coords = d02_grid[d02_filter,]
coord_sam = coords
mesh1<-inla.mesh.create(loc=as.matrix(coord_sam))

#load FANOVA results
load('FANOVA_Results_Dallas1.0_IID_D03Only.RData')
#fanova_results = readRDS('FANOVA_Results_Dallas1.0_IID_D03Only.rds')

####################################
### Dallas FANOVA Postprocessing ###
####################################


#declare list for results
fanova_results = list()
fanova_results[['coords']] = coords


#get data ready for i.i.d runs
wanted_data = list()
total_time = max(fanova_dat$time)
for (i in 1:total_time)
{
  wanted_data[[i]] = fanova_dat[which(fanova_dat$time == i),]
}


#compile results - BETA1 (CITY)
timepoints = max(fanova_dat$time)
actual_range = 1:dim(coords)[1]
all_betas = matrix(NaN,
                   nrow = dim(coords)[1],
                   ncol = timepoints)
all_lower = matrix(NaN,
                   nrow = dim(coords)[1],
                   ncol = timepoints)
all_upper = matrix(NaN,
                   nrow = dim(coords)[1],
                   ncol = timepoints)
for(t in 1:timepoints)
{
  full_ID = inla_output[[t]]$b1.field$ID
  ID = unique(inla_output[[t]]$b1.field$ID)
  betas = inla_output[[t]]$b1.field$mean
  lb = inla_output[[t]]$b1.field$`0.025quant`
  ub = inla_output[[t]]$b1.field$`0.975quant`
  estimate<-c()
  upper = c()
  lower = c()
  #print(paste(t, length(betas)/replicates))
  
  for(l in 1:length(ID)){
    index = which(full_ID == ID[l])
    estimate[l] = mean(betas[index])
    lower[l] = mean(lb[index])
    upper[l] = mean(ub[index])
  }
  
  all_betas[ID,t] = estimate
  all_lower[ID,t] = lower
  all_upper[ID,t] = upper
  
  print(t)
}
dont_want = c(1,2)
avg_beta = apply(all_betas[,-dont_want], 1, function(x) mean(x, na.rm = TRUE))
avg_lower = apply(all_lower[,-dont_want], 1, function(x) mean(x, na.rm = TRUE))
avg_upper = apply(all_upper[,-dont_want], 1, function(x) mean(x, na.rm = TRUE))

city_results = cbind(avg_beta, avg_lower, avg_upper)
colnames(city_results) = c('mean', 'lower', 'upper')

fanova_results[['cityflag']] = city_results





# #compile results - BETA2 (RESO)
# timepoints = max(fanova_dat$time)
# actual_range = 1:dim(coords)[1]
# all_betas = matrix(NaN,
#                    nrow = dim(coords)[1],
#                    ncol = timepoints)
# all_lower = matrix(NaN,
#                    nrow = dim(coords)[1],
#                    ncol = timepoints)
# all_upper = matrix(NaN,
#                    nrow = dim(coords)[1],
#                    ncol = timepoints)
# for(t in 1:timepoints)
# {
#   full_ID = inla_output[[t]]$b2.field$ID
#   ID = unique(inla_output[[t]]$b2.field$ID)
#   betas = inla_output[[t]]$b2.field$mean
#   lb = inla_output[[t]]$b2.field$`0.025quant`
#   ub = inla_output[[t]]$b2.field$`0.975quant`
#   estimate<-c()
#   upper = c()
#   lower = c()
#   #print(paste(t, length(betas)/replicates))
#   
#   for(l in 1:length(ID)){
#     index = which(full_ID == ID[l])
#     estimate[l] = mean(betas[index])
#     lower[l] = mean(lb[index])
#     upper[l] = mean(ub[index])
#   }
#   
#   all_betas[ID,t] = estimate
#   all_lower[ID,t] = lower
#   all_upper[ID,t] = upper
#   
#   #print(t)
# }
# dont_want = c(1,2)
# avg_beta = apply(all_betas[,-dont_want], 1, function(x) mean(x, na.rm = TRUE))
# avg_lower = apply(all_lower[,-dont_want], 1, function(x) mean(x, na.rm = TRUE))
# avg_upper = apply(all_upper[,-dont_want], 1, function(x) mean(x, na.rm = TRUE))
# 
# boundary_results = cbind(avg_beta, avg_lower, avg_upper)
# colnames(boundary_results) = c('mean', 'lower', 'upper')
# 
# fanova_results[['reso']] = boundary_results





#compile results - BETA3 (PHYSICS)
timepoints = max(fanova_dat$time)
actual_range = 1:dim(coords)[1]
all_betas = matrix(NaN,
                   nrow = dim(coords)[1],
                   ncol = timepoints)
all_lower = matrix(NaN,
                   nrow = dim(coords)[1],
                   ncol = timepoints)
all_upper = matrix(NaN,
                   nrow = dim(coords)[1],
                   ncol = timepoints)
for(t in 1:timepoints)
{
  full_ID = inla_output[[t]]$b3.field$ID
  ID = unique(inla_output[[t]]$b3.field$ID)
  betas = inla_output[[t]]$b3.field$mean
  lb = inla_output[[t]]$b3.field$`0.025quant`
  ub = inla_output[[t]]$b3.field$`0.975quant`
  estimate<-c()
  upper = c()
  lower = c()
  #print(paste(t, length(betas)/replicates))
  
  for(l in 1:length(ID)){
    index = which(full_ID == ID[l])
    estimate[l] = mean(betas[index])
    lower[l] = mean(lb[index])
    upper[l] = mean(ub[index])
  }
  
  all_betas[ID,t] = estimate
  all_lower[ID,t] = lower
  all_upper[ID,t] = upper
  
  print(t)
}
dont_want = c(1,2)
avg_beta = apply(all_betas[,-dont_want], 1, function(x) mean(x, na.rm = TRUE))
avg_lower = apply(all_lower[,-dont_want], 1, function(x) mean(x, na.rm = TRUE))
avg_upper = apply(all_upper[,-dont_want], 1, function(x) mean(x, na.rm = TRUE))

physics_a_results = cbind(avg_beta, avg_lower, avg_upper)
colnames(physics_a_results) = c('mean', 'lower', 'upper')


fanova_results[['physics']] = physics_a_results





print(max(ID))
print(dim(coords))

saveRDS(fanova_results, file = 'Dallas1.0_FANOVA_Results_D03Only.rds')