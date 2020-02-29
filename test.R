source("RFoccu.R")


set.seed(42)

nsite = 300
occu_x = rnorm(nsite)


occu_rate = exp(-0.5 + 2*occu_x)/(1+exp(-0.5+2*occu_x))

Z=runif(nsite)<=occu_rate

#Z = occu_x>=-0.5 & occu_x<=0.8

nperiod = 5

det_x = list()
dethist = matrix(NA,nsite,nperiod)

for(i in 1:nperiod){
  det_x_temp=rnorm(nsite)
  
  det_rate = exp(-0.5 + det_x_temp)/(1+exp(-0.5+det_x_temp))
  
  dethist[,i] = (runif(nsite)<=det_rate) * Z 
  
  det_x[[i]] = cbind(occu_x,det_x_temp)
}


test_RFoccu = RFoccu(y=dethist, matrix(occu_x), det_x, occu_x_new=NULL, 
                  burn_in=1000,n_sample = 10000,thin_by=10)
  


