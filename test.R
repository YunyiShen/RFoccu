source("RFoccu.R")
source("GAMoccu.R")
require(coda)
require(mgcv)
require(randomForest)
require(svMisc)

set.seed(42)

nsite = 100
occu_x = runif(nsite,-3,3)


occu_rate = exp( 1 - 4*occu_x^2+3*occu_x)/(1+exp(1-4*occu_x^2+3*occu_x))
plot(occu_x,occu_rate)


Z=runif(nsite)<=occu_rate

#Z = occu_x>=-0.5 & occu_x<=0.8

nperiod = 3

det_x = list()
dethist = matrix(NA,nsite,nperiod)

for(i in 1:nperiod){
  det_x_temp=rnorm(nsite)
  
  det_rate = exp(-1+ det_x_temp - occu_x)/(1+exp(-1+det_x_temp - occu_x))
  
  dethist[,i] = (runif(nsite)<=det_rate) * Z 
  
  det_x[[i]] = cbind(occu_x,det_x_temp)
}


test_RFoccu = RFoccu(y=dethist, matrix(occu_x), 
                     det_x, det.formula = y~ti(det_x_temp,occu_x),
                     occu_x_new=NULL, 
                     burn_in=1000,n_sample = 5000,thin_by=1)
  
test_GAMoccu = GAMoccu(y=dethist, 
                       data.frame(occu_x), 
                       occu.formular = y~ti(occu_x),
                       det_x, 
                       det.formula = y~ti(det_x_temp,occu_x),
                       occu_x_new=NULL, 
                       burn_in=1000,n_sample = 5000,thin_by=10)

Baseline_GAM = gam(y~ti(occu_x),data = data.frame(y=1*(rowSums(dethist)>0),occu_x),family = binomial())
Baseline_true_GAM = gam(y~ti(occu_x),data = data.frame(y=Z,occu_x),family = binomial())


plot(occu_x[order(occu_x)],occu_rate[order(occu_x)],type = "line",col="blue")
lines(occu_x[order(occu_x)],Baseline_true_GAM$fitted.values[order(occu_x)],col="green")
points(occu_x,colMeans(test_GAMoccu$psi),col="red")
points(occu_x,Baseline_GAM$fitted.values)
