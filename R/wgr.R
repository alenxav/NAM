
wgr = function(y,X,
               it=1500,bi=500,th=1,
               bag=1,rp=FALSE,
               iv=FALSE,de=FALSE,pi=0,
               df=5,R2=0.5,
               eigK=NULL,VarK=0.95,
               verb=FALSE){
  
  if(de) iv=TRUE
  anyNA = function(x) any(is.na(x))
  # Imputation with expectation
  if(anyNA(X)){
    cat('Imputing missing genotypes\n')
    imp = function(x){
      x[is.na(x)] = mean(x,na.rm=TRUE)
      x[is.nan(x)] = 0
      return(x)}
    X = apply(X,2,imp)}
  if(bag!=1) df = df/(bag^2)
  gen0 = X
  # Polygenic term
  if(!is.null(eigK)){
    V = eigK$values
    pk = which.max((cumsum(V)/length(V))>VarK)
    U0 = U = eigK$vectors[,1:pk]
    V = V[1:pk]
    H = h = rep(0,pk)
    dh = rep(0,pk)
    Vk = rep(1,pk)
    xxK = rep(bag,pk)
  }
  # Remove missing y's
  if(anyNA(y)){
    mis = which(is.na(y))
    y = y[-mis]
    X = X[-mis,]
    if(!is.null(eigK)) U = U[-mis,]
  }
  # MCMC settings
  post = seq(bi,it,th)
  mc = length(post)
  # Preparing inputs
  n = nrow(X)
  p = ncol(X)
  xx = apply(X,2,crossprod)*bag
  b = rep(0,p)
  d = rep(1,p)
  mu = mean(y - X%*%b)
  e = y - mu
  MSx = sum(apply(X,2,var,na.rm=T))
  Va = MSx
  Vb = rep(Va,p)
  Ve = 1
  L = Vb/Ve
  vy = var(y,na.rm=T)
  # Priors 10/15/17
  Sb = (R2)*df*vy/MSx;
  if((pi>0)&(!iv)){Sb = Sb/(1-pi)}
  Se = (1-R2)*df*vy;
  if(!is.null(eigK)) Sk = R2*var(y,na.rm=T)*df
  # Storing Posterior
  B0 = VA = VE = VP = S = 0
  VB = D = B = rep(0,p)
  #RUN
  if(verb) pb = txtProgressBar(style = 3)
  for(i in 1:it){
    # Resampling
    if(bag!=1) Use = sort(sample(n,n*bag,rp))-1
    # Update polygenic term and regression coefficients
    if(!is.null(eigK)){
      Lk = Ve/(V*Vk)
      if(bag!=1){
        update = KMUP2(U,Use,h,dh,xxK,e,Lk,Ve,0)
      }else{
        update = KMUP(U,h,dh,xxK,e,Lk,Ve,0)
      }
      h = update[[1]]
      e = update[[3]]
      if(bag!=1){update = KMUP2(X,Use,b,d,xx,e,L,Ve,pi)}else{update = KMUP(X,b,d,xx,e,L,Ve,pi)}
      if(pi>0) d = update[[2]]
      b = update[[1]]
      e = update[[3]]
    }else{
      # Update regression coefficients without polygene
      if(bag!=1){update = KMUP2(X,Use,b,d,xx,e,L,Ve,pi)}else{update = KMUP(X,b,d,xx,e,L,Ve,pi)}
      if(pi>0) d = update[[2]]
      b = update[[1]]
      e = update[[3]]
    }
    # Update genetic variance
    if(iv){
      # Variable selection?
      if(pi>0){
        q = d; q[q==0]=pi
        if(de){
          # Laplace?
          Vb = sqrt( b^2*Ve/MSx )
        }else{
          # T?
          Vb = (Sb+(b)^2)/rchisq(p,df+1)
        }
        # All-in?
      }else{
        # Laplace?
        if(de){
          Vb = sqrt( b^2*Ve/MSx )
        }else{
          # T?
          Vb = (Sb+b^2)/rchisq(p,df+1)
        }
      }
    }else{
      Va = (crossprod(b)+Sb)/rchisq(1,df+p)
      Vb = rep(Va,p)
    }
    if(!is.null(eigK)){
      Vp = (sum(h^2/V)+Sk)/rchisq(1,df+pk)
      Vk = rep(Vp,pk)
    }
    # Residual variance
    Ve = (crossprod(e)+Se)[1,1]/rchisq(1,n*bag+df)
    L = Ve/Vb;
    # Intercept
    if(!is.null(eigK)){e = y-mu-X%*%b-U%*%h}else{e = y-mu-X%*%b}
    mu0 = rnorm(1,mean(e), Ve/n )
    mu = mu+mu0
    e = e-mu0
    # Save posterior
    if(i%in%post){
      B0 = B0+mu;
      B = B+b
      D = D+d
      VE = VE+Ve
      if(iv){VB = VB+Vb}else{VA = VA+Va}
      if(!is.null(eigK)){H = H+h; VP = VP+Vp}
    }
    if(verb) setTxtProgressBar(pb, i/it)
  }
  if(verb) close(pb)
  # Posterior mean
  B0 = B0/mc
  D = D/mc
  B = B/mc/mean(D)
  VE = VE/mc
  if(iv){VB = VB/mc}else{VA = VA/mc}
  if(!is.null(eigK)){H = H/mc; VP = VP/mc}
  # Fitted values
  if(!is.null(eigK)){
    poly = U0 %*% H
    HAT = B0 + gen0 %*% B + poly
  }else{
    HAT = B0 + gen0 %*% B
  }
  # Output
  if(!is.null(eigK)){
    if(iv){
      final = list('mu'=B0,'b'=B,'Vb'=VB,'d'=D,'Ve'=VE,'hat'=HAT,'u'=poly,'Vk'=VP,'cxx'=mean(xx))
    }else{
      final = list('mu'=B0,'b'=B,'Vb'=VA,'d'=D,'Ve'=VE,'hat'=HAT,'u'=poly,'Vk'=VP,'cxx'=mean(xx))
    }
  }else{
    if(iv){
      final = list('mu'=B0,'b'=B,'Vb'=VB,'d'=D,'Ve'=VE,'hat'=HAT,'cxx'=mean(xx))
    }else{
      final = list('mu'=B0,'b'=B,'Vb'=VA,'d'=D,'Ve'=VE,'hat'=HAT,'cxx'=mean(xx))
    }
  }
  return(final)
}

ben = function(y,gen,it=750,bi=250,th=1,bag=0.80,alpha=0.5,wpe=50,MH=FALSE,verb=TRUE){
  
  X=gen;
  rm(gen);
  
  # Function to update beta
  upB = function(e,mu,X,b,l,a,xx,p,E2,X2,bag,pi,wpe,ve,O,mh=TRUE){
    xx = xx*bag
    pi = (0.5+a*99)/100
    if(mh){
      a_new = rbeta(1,200*pi,200*(1-pi))
      b1 = b2 = rep(NA,p)
      e1 = e2 = e
      L1_1 = (1-a)*l/(2*xx)
      L2_1 = (a)*(xx/(xx*l))
      L1_2 = (1-a_new)*l/(2*xx)
      L2_2 = (1-a_new)*(xx/(xx*l))
      # Update regression coefficients
      for(j in O){
        # Old alpha
        xy1 = (crossprod(e1,X[,j])+b[j]*xx[j])/xx[j]
        s1 = sign(xy1)
        beta = abs(xy1)-L1_1[j]
        ben = s1*beta*L2_1[j]
        b1[j] = ifelse(beta>0,rnorm(1,ben,ve/(xx[j])),0)
        e1 = e1 - X[,j]*(b1[j]-b[j])
        # New alpha
        xy2 = (crossprod(e2,X[,j])+b[j]*xx[j])/xx[j]
        s2 = sign(xy2)
        beta = abs(xy2)-L1_2[j]
        ben = s2*beta*L2_2[j]
        b2[j] = ifelse(beta>0,rnorm(1,ben,ve/(xx[j])),0)
        e2 = e2 - X[,j]*(b2[j]-b[j])
      }
      # Loss function
      SSPE_1 = sum(as.vector(tcrossprod(b1,X2)-E2)^2)
      SSPE_2 = sum(as.vector(tcrossprod(b2,X2)-E2)^2)
      LOSS1 = wpe*SSPE_1+crossprod(e1)+l*(0.5*crossprod(b1)*(1-a)+sum(abs(b1))*a)
      LOSS2 = wpe*SSPE_2+crossprod(e2)+l*(0.5*crossprod(b2)*(1-a_new)+sum(abs(b2))*a_new)
      LR = LOSS2/LOSS1
      if(is.na(LR)|is.nan(LR)) LR=0
      if(LR>1){
        P=list('b'=b2,'a'=a_new,'e'=e2,'oob'=SSPE_2)
      }else{
        if(MH){
          # Metropolis-Hastings
          if(LR>runif(1)){
            P=list('b'=b2,'a'=a_new,'e'=e2,'oob'=SSPE_2)
          }else{
            P=list('b'=b1,'a'=a,'e'=e1,'oob'=SSPE_1)
          }
        }else{
          # Acceptance-Rejection
          P=list('b'=b1,'a'=a,'e'=e1,'oob'=SSPE_1)
        }
      }
    }else{
      b1 = rep(NA,p)
      e1 = e
      L1_1 = (1-a)*l/(2*xx)
      L2_1 = (a)*(xx/(xx*l))
      # Update regression coefficients
      for(j in O){
        # Old alpha
        xy1 = (crossprod(e1,X[,j])+b[j]*xx[j])/xx[j]
        s1 = sign(xy1)
        beta = abs(xy1)-L1_1[j]
        ben = s1*beta*L2_1[j]
        b1[j] = ifelse(beta>0,rnorm(1,ben,ve/(xx[j])),0)
        e1 = e1 - X[,j]*(b1[j]-b[j])}
      # Loss function
      SSPE = sum(as.vector(tcrossprod(b1,X2)-E2)^2)
      P=list('b'=b1,'a'=a,'e'=e1,'oob'=SSPE)
    }
    return(P)
  }
  
  
  # Missing
  if(anyNA(y)){
    mis = which(is.na(y))
    Y = y
    XX = X[mis,]
    y = y[-mis]
    X = X[-mis,]
    MISS = TRUE
  }else{
    MISS = FALSE
  }
  # Data
  xx = apply(X,2,function(x)crossprod(x))
  b0 = crossprod(X,y)[,1]/xx
  O = order(b0^2,decreasing = TRUE)
  n = nrow(X)
  p = ncol(X)
  bn = round(n*bag)
  MCMC = seq(bi,it,th)
  MC = length(MCMC)
  # Parameters
  mu = mean(y)
  e = y-mu
  b = rep(0,p)
  a = alpha
  l = 1
  ve = 0.1
  # Store posterior
  B = rep(0,p)
  MU = A = L = SSPE = 0
  if(verb) pb = txtProgressBar(style = 3)
  # Loop
  for(i in 1:it){
    
    # Bagging
    s = sort(sample(1:n,n-bn,replace=FALSE))
    # UPDATE
    UP = upB(e[-s],mu,X[-s,],b,l,a,xx*bag,p,e[s],
             X[s,],bag,pi=a,wpe,ve,O=O,mh=i%%10==0)
    b = UP[[1]]
    a = UP[[2]]
    e = UP[[3]]
    
    mu = mu + mean(e)
    df_prior = 2+rpois(1,3)
    Se = runif(1,0,1)
    Sb = runif(1,0,1)*df_prior
    ve = (crossprod(e)+Se)/(n+2)
    vb = (crossprod(b)+Sb)/(p+df_prior)
    l = ve/vb
    e = as.vector(y-(mu+tcrossprod(b,X)))
    
    # STORE
    if(i%in%MCMC){
      B = B+b
      MU = MU+mu
      A = A+a
      L = L+l
      if(bag<1) SSPE = SSPE+UP$oob/(n-bn)
    }
    if(verb) setTxtProgressBar(pb, i/it)
  }
  if(verb) close(pb)
  # Posterior
  Bhat = B/MC
  MUhat = MU/MC
  Ahat = A/MC
  Lhat = L/MC
  MSPEout = mean(SSPE)/MC
  # Prediction
  if(MISS){
    Yhat = Y
    Yhat[-mis] = as.vector(mu+tcrossprod(Bhat,X))
    Yhat[mis] = as.vector(mu+tcrossprod(Bhat,XX))
  }else{
    Yhat = as.vector(mu+tcrossprod(Bhat,X))
  }
  # OUTPUT
  LIST = list('hat'=Yhat,'coef'=Bhat,'b0'=MUhat,
              'alp'=Ahat,'lmb'=Lhat,'MSPEoob'=MSPEout)
  
}
