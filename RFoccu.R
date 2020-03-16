RFoccu = function(y, occu_x, det_x,det.formular, occu_x_new=NULL, 
                  burn_in=1000,n_sample = 10000,thin_by=10,...){
  # @y: a detection histroy matrix, row as site and col as period, no NA for now
  # @x_occu x_det: 
  #  design matrix for occupancy and detection, detection part should be a list
  # @x_new: new design matrix to project on, occupancy only for now
  # @burn_in burn in 
  # @n_sample: how many posterior sample needed
  # @thin_by: thinning
  # @...: ... for randomforest
  require(randomForest)
  require(mgcv)
  require(coda)
  
  psi = mcmc(matrix(NA,floor(n_sample/thin_by),nrow(y)),thin = thin_by)
  Z = psi
  if(!is.null(occu_x_new)) psi_new = mcmc(matrix(NA,floor(n_sample/thin_by),nrow(x_new)),thin = thin_by)
  p = mcmc(matrix(NA,floor(n_sample/thin_by),nrow(y)*ncol(y)),thin = thin_by)
  
  Z_curr = 1 * (rowSums(y)>0)
  p_curr_mat = y
  psi_curr = Z_curr
  Z_abs = Z_curr
  det_x_long = Reduce(rbind,det_x)
  k = 1
  for (i in 1:(burn_in+n_sample)){
    svMisc::progress(((i-1)/(n_sample + burn_in))*100,progress.bar = T)
    # find psi:
    RF_occu_curr = randomForest(x=occu_x,y=Z_curr,...)
    psi_curr = predict(RF_occu_curr,newdata = occu_x)
    psi_curr = psi_curr * (psi_curr>=0)
    if(!is.null(occu_x_new)) psi_new_curr = predict(RF_occu_curr,newdata = occu_x_new)
    Z_curr = runif(length(psi_curr))<=psi_curr
    
    
    # now, find posterior of p:
    red_det_x = lapply(det_x,function(w,Z_curr){
      w[Z_curr==1,]
    },Z_curr) # only those site occupied was considered in training
    
    red_det_his = y[Z_curr==1,] # only get the detection history when Z=1 (red for reduced)
    red_det_x_long = Reduce(rbind,red_det_x) # get the long version design matrix

    red_det_his_long = matrix((red_det_his),length(red_det_his)) # long response 
    # is by period i.e. first row was site1 period1 and second row was site 2 period 1
    rm(red_det_x,red_det_his)
    # fit detection discrimant model, logistic here:
    #RF_det_curr = randomForest(x=red_det_x_long,y=red_det_his_long,...)
    #glm_det_curr = summary( glm(red_det_his_long~red_det_x_long-1,family = binomial))$coefficients
    gam_data_temp = data.frame(y=red_det_his_long,red_det_x_long)
    gam_det_curr = gam(det.formular,family = binomial,data = gam_data_temp)
    
    # laplacian approximation with flat prior:
    p_curr_pred = predict(gam_det_curr,newdata = data.frame(det_x_long),se.fit = T)
    logit_p_sample_laplace = rnorm(nrow(det_x_long),p_curr_pred$fit,p_curr_pred$se.fit)
    p_curr = 1/(1+exp(-logit_p_sample_laplace))
    # current detection probability estimation:
    #p_curr = predict(glm_det_curr,newdata = det_x_long)
    p_curr_mat = (matrix(p_curr,ncol = ncol(y)))  # matrix version 
    
    #likelihood of seeing the detection history
    joint_p_curr_mat = exp(rowSums(y*log(p_curr_mat+1e-12)+(1-y)*log(1-p_curr_mat+1e-12)))
    
    # now sample posterior Z
    ## remember this later one only work when no detection at such site
    psi_posterior=joint_p_curr_mat*psi_curr/
                  ((1-psi_curr)+joint_p_curr_mat*psi_curr)
    Z_temp = runif(nrow(y))<psi_posterior
    Z_curr[Z_abs==0] = Z_temp[Z_abs==0]
    
    
    # saving current sample
    if(i>burn_in & (i-burn_in)%% thin_by==0) {
      p[k,] = p_curr
      psi[k,]=psi_curr
      Z[k,]=Z_curr
      if(!is.null(occu_x_new)) psi_new[k,]=psi_new_curr
      k=k+1
    }
  }
  
  if(is.null(occu_x_new)) psi_new=NULL
  list(psi=psi, p=p, Z=Z, psi_new=psi_new)
}