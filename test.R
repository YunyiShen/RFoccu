source("RFoccu.R")


set.seed(42)
nsite = 500
occu_x =cbind(1, rnorm(nsite))
Z = 1*( occu_x>=-0.5 & occu_x<=0.8)


nperiod = 3

det_x = list()
dethist = matrix(NA,nsite,nperiod)

for(i in 1:nperiod){
  det_x_temp=rnorm(nsite)
  
  dethist[,i] = (det_x_temp<=-.5) * Z 
  
  det_x[[i]] = cbind(occu_x,det_x_temp)
}


test_RFoccu = RFoccu(y=dethist, matrix(occu_x), det_x, occu_x_new=NULL, 
                  burn_in=1000,n_sample = 10000,thin_by=1)
  


