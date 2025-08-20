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
#fanova_dat = readRDS('FANOVA_DFW_March_PostCity3HOURPRCP1mm_D02_FINAL.rds')
load('FANOVA_DFW_March_PostCityALLTIMEPRCP_D02.RData')
fanova_dat = fanova_dat %>%
  filter(lon < -96, lon > -97.5, time > 12)
fanova_dat$time = fanova_dat$time - 12

#load('FANOVA_DataMarchPRCP_Dallas.RData')
#fanova_dat = fanova_dat %>% filter(repl != 1)
#fanova_dat$repl = fanova_dat$repl - 1
#fanova_dat$combined = fanova_dat$repl - 2
#fanova_dat$combined = factor(fanova_dat$combined, levels = c(-1,0,1), labels = c('citynoemiss', 'nocityemiss', 'cityemiss'))
replicates = length(unique(fanova_dat$repl))
head(fanova_dat)
summary(fanova_dat)

#load coordinate grids and filter down to D03 space
lat_d02 = c(readRDS('lat_grid_d01')) #can also try d02 here
lon_d02 = c(readRDS('lon_grid_d01')) #can also try d02 here
d02_grid = cbind(lon_d02, lat_d02)
total_d02 = dim(d02_grid)[1]
filter1 = d02_grid[,1] >= -98.5 & d02_grid[,1] <= -95
filter2 = d02_grid[,2] >= 31.5 & d02_grid[,2] <= 34
d02_filter = filter1 & filter2
d02_grid = d02_grid[d02_filter,]


#specify/print number of cores
options(cores = 24)
inla.getOption("num.threads")
threads_list = c("24:1", "12:1", "8:1")
threads = threads_list[index]


#set inla mode
inla.setOption(inla.mode="experimental")



####################
### FANVOA Model ###
####################



#start inla and SPDE code
coord_sam = d02_grid
mesh1<-inla.mesh.create(loc=as.matrix(coord_sam))
SPDE = inla.spde2.matern(mesh1,
                         alpha = 2) #I can set priors here for range, etc. What should they be?



#Run time points i.i.d
start = proc.time()


#declare locations
locations = cbind(fanova_dat$lon, fanova_dat$lat)

#declare variables
#x1 = as.double(wanted_data[[t]]$combined)
x1 = as.double(fanova_dat$cityflag)
x2 = as.double(fanova_dat$emissions)
Y = as.numeric(log(fanova_dat$precipitation+0.01)) #maybe take a log here



#build INLA model components
A1 = inla.spde.make.A(SPDE, loc = locations, repl = fanova_dat$repl, n.repl = replicates, weights = x1, group = fanova_dat$time)
A2 = inla.spde.make.A(SPDE, loc = locations, repl = fanova_dat$repl, n.repl = replicates, weights = x2, group = fanova_dat$time)
spatial.idx.1 = inla.spde.make.index("b1.field", n.spde=SPDE$n.spde, n.repl = replicates, n.group = max(fanova_dat$time))
spatial.idx.1 = lapply(spatial.idx.1, as.numeric)
spatial.idx.2 = inla.spde.make.index("b2.field", n.spde=SPDE$n.spde, n.repl = replicates, n.group = max(fanova_dat$time))
spatial.idx.2 = lapply(spatial.idx.2, as.numeric)
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

stk.data=inla.stack(stk.b1,stk.b2)
#stk.data = inla.stack(stk.b1)





#declare formula --- try both model = "iid" and model = SPDE ---- WHAT SHOULD WE CHOOSE FOR THE LATENT MODEL


formula = y ~ -1 +
  b1.intercept +
  f(b1.field, model = 'iid', replicate = b1.field.repl) +
  b2.intercept +
  f(b2.field, model = 'iid', replicate = b2.field.repl)
  
  
#formula = y ~ -1 +
#b1.intercept +
#f(b1.field, model = 'iid', replicate = b1.field.repl)


#run inla
result = inla(formula,
              family="normal", #swicth to normal
              data=inla.stack.data(stk.data),
              control.predictor = list(A=inla.stack.A(stk.data),compute=TRUE),
              control.compute = list(dic=TRUE, openmp.strategy="huge"),
              verbose = FALSE,
              num.threads = threads)

inla_output = result$summary.random
  

proc.time() - start



#save output
save(inla_output, file = 'FANOVA_Results_DFW_March_PostCityLASTHOURPRCP_IID_D01Mesh_LOG.RData')
saveRDS(inla_output, 'FANOVA_Results_DFW_March_PostCityLASTHOURPRCP_IID_D01Mesh_LOG.rds')

