source("RFoccu.R")
source("GAMoccu.R")
require(coda)
require(mgcv)
require(randomForest)
require(svMisc)

set.seed(114514)

nsite = 100
occu_x = rnorm(2*nsite,0,1)
occu_x = matrix(occu_x,ncol = 2)
colnames(occu_x) = c("x1","x2")

occu_rate = exp( 1 - 4*occu_x[,1]^2+3*occu_x[,2])/(1+exp(1-4*occu_x[,1]^2+3*occu_x[,2]))
plot(occu_x,occu_rate)


Z=runif(nsite)<=occu_rate

#Z = occu_x>=-0.5 & occu_x<=0.8

nperiod = 3

det_x = list()
dethist = matrix(NA,nsite,nperiod)

for(i in 1:nperiod){
  det_x_temp=rnorm(nsite)
  
  det_rate = exp( -0.5+ det_x_temp+occu_x[,1]- occu_x[,2]^2)/(1+exp(-0.5+det_x_temp+occu_x[,1]- occu_x[,2]^2))
  
  dethist[,i] = (runif(nsite)<=det_rate) * Z 
  
  det_x[[i]] = cbind(occu_x,det_x_temp)
}


test_RFoccu = RFoccu(y=dethist, data.frame(occu_x), 
                     det_x, det.formula = y~ti(det_x_temp)+ti(x1)+ti(x2),
                     occu_x_new=NULL, 
                     burn_in=1000,n_sample = 5000,thin_by=1)
  
test_GAMoccu = GAMoccu(y=dethist, 
                       data.frame(occu_x), 
                       occu.formular = y~ti(x1)+ti(x2),
                       det_x, 
                       det.formula = y~ti(det_x_temp)+ti(x1)+ti(x2),
                       occu_x_new=NULL, 
                       burn_in=1000,n_sample = 5000,thin_by=10)

Baseline_GAM = gam(y~ti(x1)+ti(x2),data = data.frame(y=1*(rowSums(dethist)>0),occu_x),family = binomial())
#Baseline_GAM = gam(y~ti(occu_x),data = data.frame(y=1*(rowSums(dethist)>0),occu_x),family = binomial())
Baseline_true_GAM = gam(y~ti(x1)+ti(x2),data = data.frame(y=Z,occu_x),family = binomial())


plot(occu_rate[order(occu_rate)],occu_rate[order(occu_rate)],type = "line",col="blue")
lines(occu_rate[order(occu_rate)],Baseline_true_GAM$fitted.values[order(occu_rate)],col="green")
points(occu_rate,colMeans(test_GAMoccu$psi),col="red")
points(occu_rate,colMeans(test_RFoccu$psi),col="darkred")
points(occu_rate,Baseline_GAM$fitted.values)

roc_RF = roc(Z,colMeans(test_RFoccu$psi))
roc_GAM = roc(Z,colMeans(test_GAMoccu$psi))
roc_naiveGAM = roc(Z,Baseline_GAM$fitted.values)
roc_trueGAM = roc(Z,Baseline_true_GAM$fitted.values)




