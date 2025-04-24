conf.aft.fun = function(x,y,delta,dat,
                        test.t, test.x,dat_test,
                        alpha, side, dist="lognormal",B=1000){
  
  ## Check the dimension of the covariates
  if(is.null(dim(x)[1])){
    len_x <- length(x)
    p <- 1
  }else{
    len_x <- dim(x)[1]
    p <- dim(x)[2]
  }
  xnames <- paste0("x",1:p)
  colnames(x) = xnames  
  
  ## prepare model formula
  fmla <- as.formula(paste("Surv(y, delta) ~ ", paste(xnames, collapse= "+")))
  
  # KM estimate for the censoring, no covariates scenario
  cox.censor=survfit(Surv(y,(1-delta))~1, data=dat) 
  
  # estimate weights for the failed subjects
  dat.event = dat[dat$delta==1,]
  n1 = nrow(dat.event)
  tto = dat.event$y  
  p.tto = summary(cox.censor,tto)$surv
  mass.tto = 1/p.tto
  mass.tto = mass.tto/sum(mass.tto)  
  
  id0=1:n1
  
  ########################################
  #### conformal AFT
  ########################################
  
  z=rep(NA,B)
  fit.aft = survreg(fmla, data=dat,dist=dist) 
  
  for(b in 1:B)
  {
    boot.dat = dat[sample(nrow(dat), replace=TRUE),]
    boot.aft = survreg(fmla, data=boot.dat,dist=dist) 
    
    # new data point sampling from the uncensored population
    tid=sample(id0,size=1,replace=T,prob=mass.tto)
    new.dat = dat.event[tid,]
    
    new.xb = predict(boot.aft, newdata=new.dat, type='lp')
    z[b]= 1-pnorm((log(new.dat$y)-new.xb)/boot.aft$scale)
    
  }
  
  z= na.omit(z)
  z=sort(z)
  
  ## the fitted quantile for the testing data 
  if(side=="two"){
    u = quantile(z,1-alpha/2)
    l = quantile(z,alpha/2)
    
    test.xb = predict(fit.aft, newdata=dat_test, type='lp')
    temp = 1-pnorm((log(test.t)-test.xb)/fit.aft$scale)
    cp = mean(ifelse(temp<=u&temp>=l,1,0))
    
  }
  
  if(side=="one"){
    u=z[round(B*(1-alpha))]
    
    test.xb = predict(fit.aft, newdata=dat_test, type='lp')
    temp = 1-pnorm((log(test.t)-test.xb)/fit.aft$scale)
    cp = mean(ifelse(temp<=u,1,0))
    
  }
  
  ## return coverage probability
  return(cp)
}
   


conf.cox.fun = function(x,y,delta,dat,
                        test.t, test.x,
                        alpha, side,B=1000){
  
  ## Check the dimension of the covariates
  if(is.null(dim(x)[1])){
    len_x <- length(x)
    p <- 1
  }else{
    len_x <- dim(x)[1]
    p <- dim(x)[2]
  }
  xnames <- paste0("x",1:p)
  colnames(x) = xnames  
  
  ## prepare model formula
  fmla <- as.formula(paste("Surv(y, delta) ~ ", paste(xnames, collapse= "+")))
  
  # KM estimate for the censoring, no covariates scenario
  cox.censor=survfit(Surv(y,(1-delta))~1, data=dat) 
  
  # estimate weights for the failed subjects
  dat.event = dat[dat$delta==1,]
  n1 = nrow(dat.event)
  tto = dat.event$y
  p.tto = summary(cox.censor,tto)$surv
  mass.tto = 1/p.tto
  mass.tto = mass.tto/sum(mass.tto)
  
  id0=1:n1
  
  z=rep(NA,B)
  
  # cox model from the original data
  cox1 = coxph(fmla, data=dat)
  beta.1=cox1$coefficients
  bh <- basehaz(cox1)
  
  bh = basehaz(cox1)
  sfun0 = stepfun(bh$time, c(0,bh$hazard))

  for(b in 1:B)
  {
    boot.dat = dat[sample(nrow(dat), replace=TRUE),]
    boot.cox1 = coxph(fmla, data=boot.dat)
    boot.coef=boot.cox1$coefficients
    boot.bh <- basehaz(boot.cox1)
    
    boot.sfun0 = stepfun(boot.bh$time, c(0,boot.bh$hazard))
    
    # new data point sampling from the censoring distribution
    tid=sample(id0,size=1,replace=T,prob=mass.tto)
    new.dat = dat.event[tid,]
    new.xb = t(boot.coef) %*% t(new.dat[,xnames])   # note p is used here
    z[b] = exp(-boot.sfun0(new.dat$y)*exp(new.xb)) 
  }
  z= na.omit(z)
  z=sort(z)
  
  ## the fitted quantile for the testing data 
  if(side=="two"){
    u = quantile(z,1-alpha/2)
    l = quantile(z,alpha/2)
    
    test.xb = as.numeric(t(beta.1) %*% t(test.x)) 
    
    temp =exp(-sfun0(test.t)*exp(test.xb)) 
    cp = mean(ifelse(temp<=u&temp>=l,1,0))
  }
  
  if(side=="one"){
    u=z[round(B*(1-alpha))]
    test.xb = as.numeric(t(beta.1) %*% t(test.x)) 
    temp =exp(-sfun0(test.t)*exp(test.xb)) 
    cp = mean(ifelse(temp<=u,1,0))
  }
  
  ## return coverage probability
  return(cp)
}




conf.weibull.fun = function(x,y,delta,dat,
                            test.t, test.x,dat_test,
                            alpha, side, dist="weibull",B=1000){
  
  ## Check the dimension of the covariates
  if(is.null(dim(x)[1])){
    len_x <- length(x)
    p <- 1
  }else{
    len_x <- dim(x)[1]
    p <- dim(x)[2]
  }
  xnames <- paste0("x",1:p)
  colnames(x) = xnames  
  
  ## prepare model formula
  fmla <- as.formula(paste("Surv(y, delta) ~ ", paste(xnames, collapse= "+")))
  
  # KM estimate for the censoring, no covariates scenario
  cox.censor=survfit(Surv(y,(1-delta))~1, data=dat) 
  
  # estimate weights for the failed subjects
  dat.event = dat[dat$delta==1,]
  n1 = nrow(dat.event)
  tto = dat.event$y
  p.tto = summary(cox.censor,tto)$surv
  mass.tto = 1/p.tto
  mass.tto = mass.tto/sum(mass.tto)
  
  id0=1:n1
  
  z=rep(NA,B)
  fit.weibull = survreg(  fmla , data=dat,dist="weibull") 
  
  for(b in 1:B)
  {
    boot.dat = dat[sample(nrow(dat), replace=TRUE),]
    boot.weibull = survreg(fmla, data=boot.dat,dist="weibull") #  
    
    # new data point sampling from the censoring distribution
    tid=sample(id0,size=1,replace=T,prob=mass.tto)
    new.dat = dat.event[tid,]
    new.xb = predict(boot.weibull, newdata=new.dat, type='lp')
    
    z[b]=1-pweibull(new.dat$y,shape=1/boot.weibull$scale,scale=exp(new.xb))
  }
  z= na.omit(z)
  z=sort(z)
  
  ## the fitted quantile for the testing data 
  if(side=="two"){
    u = quantile(z,1-alpha/2)
    l = quantile(z,alpha/2)
    test.xb = predict(fit.weibull, newdata=dat_test, type='lp')
    temp = 1-pweibull(test.t,shape=1/fit.weibull$scale,scale=exp(test.xb))
    cp = mean(ifelse(temp<=u&temp>=l,1,0))
  }
  
  if(side=="one"){
    u=z[round(B*(1-alpha))]
    test.xb = predict(fit.weibull, newdata=dat_test, type='lp')
    temp = 1-pweibull(test.t,shape=1/fit.weibull$scale,scale=exp(test.xb))
    cp = mean(ifelse(temp<=u,1,0))
  }
  
  ## return coverage probability
  return(cp)
}


