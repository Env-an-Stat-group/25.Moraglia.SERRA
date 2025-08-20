#######################
#######################
### FANOVA CRC Code ###
#######################
#######################

#clear environment
rm(list = ls());gc()


#load libraries
library(tidyverse)
library(INLA)


#load data
load('FANOVA_Data_Indy_T2_D02_City_Detrended_2.0.RData')
#load('FANOVA_Data_Indy_TKE_D02_City_Detrended.RData')
#load('FANOVA_Data_NYC_TKE_D02_City_Detrended.RData')
#load('FANOVA_Data_NYC_T2_D02_City_Detrended.RData')
replicates = 2


fanova_dat$citylocations[which(fanova_dat$citylocations == 0)]  = -1 #CHANGE THIS POSSIBLY



#load coordinates - d01 --- WE USE THIS FOR THE MESH
lat_d01 = readRDS('d01lat_CROP')
lon_d01 = readRDS('d01lon_CROP')
coords_d01 = expand.grid(lon_d01, lat_d01)



#specify/print number of cores
options(cores = 24)
inla.getOption("num.threads")

#set inla mode
inla.setOption(inla.mode="experimental")



####################
### FANVOA Model ###
####################

#load D03 grids
#lat_d03 = readRDS('Lat_grid_d03_NY')
#lon_d03 = readRDS('Lon_grid_d03_NY')
#min_lat = min(lat_d03); max_lat = max(lat_d03)
#min_lon = min(lon_d03); max_lon = max(lon_d03)

#load D02 grids
#lat_d02 = readRDS('Lat_grid_d02_NY')
#lon_d02 = readRDS('Lon_grid_d02_NY')

#check which d02 points are in D03
#filter_lat = (lat_d02 >= min_lat & lat_d02 <= max_lat)
#filter_lon = (lon_d02 >= min_lon & lon_d02 <= max_lon)
#full_filter = filter_lat & filter_lon

#restrict the D02 domain
#lat_restrict = which(apply(full_filter, 2, sum) >= 1)
#lon_restrict = which(apply(full_filter, 1, sum) >= 1)


#lat_mat = lat_d02[lon_restrict, lat_restrict]
#lon_mat = lon_d02[lon_restrict, lat_restrict]

#coords_d01 = cbind(c(lon_mat), c(lat_mat))
  
#start inla and SPDE code
coord_sam = coords_d01
mesh1<-inla.mesh.create(loc=as.matrix(coord_sam))
SPDE = inla.spde2.matern(mesh1,
                         alpha = 2) #I can set priors here for range, etc. What should they be?
                         
#SPDE = inla.spde2.matern(mesh1,B.tau=matrix(cbind(0, 1, 0, sin(pi*mesh1$loc[,1])),ncol=4),
#                         B.kappa=matrix(c(0, 0, 1, 0),nrow=1,ncol=4),
#                         theta.prior.mean=c(0, 0, 0),
#                         theta.prior.prec=c(0.1, 0.1, 0.1),
#                         prior.range.nominal = 0.25) #nonstationary declaration of SPDE


#SPDE = inla.spde2.matern(mesh1, prior.range.nominal = 2) #I can set priors here for range, etc. What should they be?

#get data ready for i.i.d runs
wanted_data = list()
total_time = max(fanova_dat$time)
for (i in 1:total_time)
{
  wanted_data[[i]] = fanova_dat[which(fanova_dat$time == i),]
}

#declare locations
locations = cbind(wanted_data[[1]]$lon, wanted_data[[1]]$lat)




#Run time points i.i.d
inla_output = list()
start = proc.time()
for(t in 1:total_time)
{
  
  #declare variables
  x1 = abs(as.double(wanted_data[[t]]$cityflag)) #remove absolute value
  x2 = as.double(wanted_data[[t]]$citylocations)
  #x3 = as.double(wanted_data[[t]]$reso_b)
  #Y = as.numeric(wanted_data[[t]]$tke)
  Y = as.numeric(wanted_data[[t]]$detrended)
  
  
  
  #build INLA model components
  A1 = inla.spde.make.A(SPDE, loc = locations, repl = wanted_data[[t]]$repl, n.repl = replicates, weights = x1)
  A2 = inla.spde.make.A(SPDE, loc = locations, repl = wanted_data[[t]]$repl, n.repl = replicates, weights = x2)
  #A3 = inla.spde.make.A(SPDE, loc = locations, repl = wanted_data[[t]]$repl, n.repl = replicates, weights = x3)
  spatial.idx.1 = inla.spde.make.index("b1.field", n.spde=SPDE$n.spde, n.repl = replicates)
  spatial.idx.2 = inla.spde.make.index("b2.field", n.spde=SPDE$n.spde, n.repl = replicates)
  #spatial.idx.3 = inla.spde.make.index("b3.field", n.spde=SPDE$n.spde, n.repl = replicates)
  stk.b1 = inla.stack(data = list(y = Y),
                      A = list(A1,1),
                      effects = list(spatial.idx.1,
                                     list(b1.intercept=rep(1,length(Y)))),
                      tag="est1")
  stk.b2 = inla.stack(data = list(y = Y),
                      A = list(A2,1),
                      effects = list(spatial.idx.2,
                                     list(b2.intercept=rep(1,length(Y)))),
                      tag="est2")
  #stk.b3 = inla.stack(data = list(y = Y),
  #                    A = list(A3,1),
  #                    effects = list(spatial.idx.3,
  #                                   list(b3.intercept=rep(1,length(Y)))),
  #                    tag="est3")
  #stk.data=inla.stack(stk.b1,stk.b2,stk.b3)
  stk.data=inla.stack(stk.b1, stk.b2)
  
  #declare formula --- try both model = "iid" and model = SPDE
  #formula = y ~ -1 +
  #  b1.intercept +
  #  f(b1.field, model = "iid", replicate = b1.field.repl) +
  #  b2.intercept +
  #  f(b2.field, model = "iid", replicate = b2.field.repl) +
  #  b3.intercept +
  #  f(b3.field, model = "iid", replicate = b3.field.repl)
  
  
  formula = y ~ -1 +
    b1.intercept +
    f(b1.field, model = "iid", replicate = b1.field.repl) +
    b2.intercept +
    f(b2.field, model = "iid", replicate = b2.field.repl) #changed from model = "iid"
  
  #run inla
  result = inla(formula,
                family="normal",
                data=inla.stack.data(stk.data),
                control.predictor = list(A=inla.stack.A(stk.data),compute=TRUE),
                control.compute = list(dic=TRUE, openmp.strategy="huge"),
                verbose = FALSE,
                num.threads = "24:1")
  
  inla_output[[t]] = result$summary.random
  
}
proc.time() - start


#save output
save(inla_output, file = 'FANOVA_T2_Results_Indy_IID_D02Only_D01Mesh_CityLocs-11.RData')

