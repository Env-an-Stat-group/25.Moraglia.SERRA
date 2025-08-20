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


#Specify array
index <<- as.numeric(Sys.getenv('SGE_TASK_ID'))

#load data
fanova_dat = readRDS('FANOVA_DataMarch_Dallas_Detrended_HRRR.rds')
#load('FANOVA_Data2.0_Dallas_Detrended.RData')
replicates = max(fanova_dat$repl)

#fanova_dat$cityflag[which(fanova_dat$cityflag == -1)]  = 0

#load coordinates - d01 --- WE USE THIS FOR THE MESH
d02_grid = readRDS('dallas_D02_grid.rds')
total_d02 = dim(d02_grid)[1]
filter1 = d02_grid[,1] >= -99.95 & d02_grid[,1] <= -93.97
filter2 = d02_grid[,2] >= 30.35 & d02_grid[,2] <= 35.26
d02_filter = filter1 & filter2
coords = d02_grid[d02_filter,]



#specify/print number of cores
options(cores = 12)
inla.getOption("num.threads")
threads_list = c("12:1", "8:1", "4:1")
threads = threads_list[index]

#set inla mode
inla.setOption(inla.mode="experimental")



####################
### FANVOA Model ###
####################



#start inla and SPDE code
coord_sam = coords
mesh1<-inla.mesh.create(loc=as.matrix(coord_sam))
SPDE = inla.spde2.matern(mesh1,
                         alpha = 2) #I can set priors here for range, etc. What should they be?



#get data ready for i.i.d runs
wanted_data = list()
total_time = max(fanova_dat$time)
for (i in 1:total_time)
{
  wanted_data[[i]] = fanova_dat[which(fanova_dat$time == i),]
}


#Run time points i.i.d
inla_output = list()
start = proc.time()
for(t in 1:total_time)
{

  #declare locations
  locations = cbind(wanted_data[[t]]$lon, wanted_data[[t]]$lat)
  
  #declare variables
  x1 = as.double(wanted_data[[t]]$cityflag)
  #x2 = as.double(wanted_data[[t]]$boundary)
  x3 = as.double(wanted_data[[t]]$physics_a)
  x4 = as.double(wanted_data[[t]]$physics_b)
  #Y = as.numeric(wanted_data[[t]]$detrended)
  Y = as.numeric(wanted_data[[t]]$reflectivity) + 0.001
  
  
  
  #build INLA model components
  A1 = inla.spde.make.A(SPDE, loc = locations, repl = wanted_data[[t]]$repl, n.repl = replicates, weights = x1)
  #A2 = inla.spde.make.A(SPDE, loc = locations, repl = wanted_data[[t]]$repl, n.repl = replicates, weights = x2)
  A3 = inla.spde.make.A(SPDE, loc = locations, repl = wanted_data[[t]]$repl, n.repl = replicates, weights = x3)
  A4 = inla.spde.make.A(SPDE, loc = locations, repl = wanted_data[[t]]$repl, n.repl = replicates, weights = x4)
  spatial.idx.1 = inla.spde.make.index("b1.field", n.spde=SPDE$n.spde, n.repl = replicates)
  #spatial.idx.2 = inla.spde.make.index("b2.field", n.spde=SPDE$n.spde, n.repl = replicates)
  spatial.idx.3 = inla.spde.make.index("b3.field", n.spde=SPDE$n.spde, n.repl = replicates)
  spatial.idx.4 = inla.spde.make.index("b4.field", n.spde=SPDE$n.spde, n.repl = replicates)
  stk.b1 = inla.stack(data = list(y = Y),
                      A = list(A1,1),
                      effects = list(spatial.idx.1,
                                     list(b1.intercept=rep(1,length(Y)))),
                      tag="est1")
  #stk.b2 = inla.stack(data = list(y = Y),
  #                    A = list(A2,1),
  #                    effects = list(spatial.idx.2,
  #                                   list(b2.intercept=rep(1,length(Y)))),
  #                    tag="est2")
  stk.b3 = inla.stack(data = list(y = Y),
                      A = list(A3,1),
                      effects = list(spatial.idx.3,
                                     list(b3.intercept=rep(1,length(Y)))),
                      tag="est3")
  stk.b4 = inla.stack(data = list(y = Y),
                      A = list(A4,1),
                      effects = list(spatial.idx.4,
                                     list(b4.intercept=rep(1,length(Y)))),
                      tag="est4")
  #stk.data=inla.stack(stk.b1,stk.b2,stk.b3,stk.b4)
  stk.data=inla.stack(stk.b1,stk.b3,stk.b4)
  
  #declare formula --- try both model = "iid" and model = SPDE
  #formula = y ~ -1 +
  #  b1.intercept +
  #  f(b1.field, model = "iid", replicate = b1.field.repl) +
  #  b2.intercept +
  #  f(b2.field, model = "iid", replicate = b2.field.repl) +
  #  b3.intercept +
  #  f(b3.field, model = "iid", replicate = b3.field.repl) +
  #  b4.intercept +
  #  f(b4.field, model = "iid", replicate = b4.field.repl)
    
  formula = y ~ -1 +
    b1.intercept +
    f(b1.field, model = "iid", replicate = b1.field.repl) +
    b3.intercept +
    f(b3.field, model = "iid", replicate = b3.field.repl) +
    b4.intercept +
    f(b4.field, model = "iid", replicate = b4.field.repl)
  
  #run inla
  result = inla(formula,
                family="gamma",
                data=inla.stack.data(stk.data),
                control.predictor = list(A=inla.stack.A(stk.data),compute=TRUE),
                control.compute = list(dic=TRUE, openmp.strategy="huge"),
                verbose = FALSE,
                num.threads = threads)
  
  inla_output[[t]] = result$summary.random
  
}
proc.time() - start


#save output
save(inla_output, file = 'FANOVA_Results_DallasMarch_IID_HRRR_ReflGamma.RData')
saveRDS(inla_output, 'FANOVA_Results_DallasMarch_IID_HRRR_ReflGamma.rds')

