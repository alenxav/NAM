
gmm = function(y,gen,dta=NULL,it=75,bi=25,th=1,model="BRR",...){
  
  # models: BRR - BayesA - GBLUP - RKHS - RF
  # spline: dta must include a column called "ID Block Row Col"
  
  if(is.null(dta)&nrow(gen)==length(y)){
    n = length(y)
    ids = sort(gsub('\\..\\.','',paste('geno',round(1:n/n,5)+0.000001,sep='.')))
    dta = data.frame("ID"=ids)
    rownames(gen) = ids
    rm(n,ids)
  }
  
  # Dealing with Genotypic information and Missing loci
  gen = data.matrix(gen)
  if(anyNA(gen)){
    cat("Removing markers with more than 80% missing \n")
    m = apply(gen,2,function(x) mean(is.na(x)))
    if(any(m>0.8)) gen = gen[,-which(m>0.8)]
    cat("Imputing missing markers with expectations \n")
    expectation = function(x){
      if(anyNA(x)) x[is.na(x)] = mean(x)
      return(x)}
    E = apply(gen,2,expectation)
    dim(E) = dim(gen)
    gen = E; rm(E)
  } 
  
  # Checking if the model is acceptable
  if(! model %in% c("BRR","BayesA","GBLUP","RKHS","RF")) stop("Model must be either 'BRR', 'BayesA','GBLUP','RKHS' or 'RF' ")
  if(model=="GBLUP"|model=="RKHS"){
    KERN = TRUE
    if(model=="GBLUP"){
      idnames = rownames(gen)
      gen = GRM(gen)
      cat("Eigendecomposing GRM for GBLUP \n")
      gen = eigen(gen,TRUE)
      V = gen$values
      gen = gen$vectors
      rownames(gen) = idnames
      colnames(gen) = paste('PC',1:ncol(gen),sep='')
    }
    if(model=="RKHS"){
      idnames = rownames(gen)
      gen = GAU(gen)
      cat("Eigendecomposing Gaussian Kernel for RKHS \n")
      gen = eigen(gen,TRUE)
      V = gen$values
      gen = gen$vectors
      rownames(gen) = idnames
      colnames(gen) = paste('PC',1:ncol(gen),sep='')
    }
    # Reduction of dimensionality
    ww = which(cumsum(V)/length(V)<.98);gen = gen[,ww];V = V[ww];rm(ww)
  }else{
    KERN = FALSE
  }
  
  cat("Checking for missing and setting up parameters \n")
  # Dealing with missing information
  if(anyNA(y)){
    w = which(is.na(y))
    y = y[-w]
    dtanm = names(dta)
    dta = data.frame(dta[-w,])
    names(dta) = dtanm; rm(dtanm)
  }
  n = length(y)
  
  # Isolating spatial information
  if(sum(c("Block","Row","Col")%in%names(dta))==3){
    spline = TRUE
    spdta = dta[,c("Block","Row","Col")]
  }else{
    spline = FALSE
    cat("No spline: 'dta' does not have columns Block-Row-Col \n")
  }
  
  # Isolating genotypic information
  id = dta[,'ID']
  
  # Isolating of fixed effects
  if(any(!names(dta)%in%c("Block","Row","Col","ID"))){
    fixed = TRUE
    dta2 = dta[,-which(names(dta)%in%c("Block","Row","Col","ID"))]
    dta2 = data.frame(dta2)
    X = model.matrix(~.-1,data=dta2)
  }else{
    fixed = FALSE
    cat("No additional effects other than splines and genomics \n")
  }
  
  # Search for NN plots
  NNsrc = function(sp=NULL,rho=1,dist=3){
    # OCTREE search function 
    NN = function(a,sp){
      # shared environment
      s0 = which(sp[,1]==a[1])
      # row neighbors
      s1 = s0[which(abs(sp[s0,2]-a[2])<(rho))]
      s2 = s1[which(abs(sp[s1,3]-a[3])<=(rho*dist))]
      # col neighbors
      s3 = s0[which(abs(sp[s0,3]-a[3])<(rho*dist*0.5))]
      s4 = s3[which(abs(sp[s3,2]-a[2])<=(rho))]
      # Neighbors only
      s5 = union(s2,s4)
      wo = which(sp[s5,1]==a[1]&sp[s5,2]==a[2]&sp[s5,3]==a[3])
      s5 = s5[-wo]
      return(s5)}
    # Set the data format
    rownames(sp) = 1:nrow(sp)
    if(is.data.frame(sp)) sp = data.matrix(sp)
    # Search for neighbors within environment
    Local_NN = apply(sp,1,NN,sp=sp)
    cat(paste('Average n# of components:',round(mean(sapply(Local_NN,length)),2)),'\n')
    # RETURN
    return(Local_NN)
  }
  
  # Create covariate from neighbor plots
  SPcov = function(NN,y) sapply(X = NN,FUN = function(x,y) mean(y[x],na.rm=TRUE),y = y,USE.NAMES = FALSE)
  
  if(isTRUE(spline)){
    # Preparing kNN matrix
    cat("Preparing sparse mapping matrix (spline) \n")
    R = NNsrc(spdta,...)
  }
  
  # Preparing genotypic mapping matrix
  cat("Generating the incidence matrix of genotypes \n")
  Z = model.matrix(~id-1)
  colnames(Z) = gsub('^id','',colnames(Z))
  Weight = colSums(Z)
  
  # Dealing with unphenotyped material
  if(any(Weight==0)){
    Z = Z[,-which(Weight==0)]
    Weight = Weight[-which(Weight==0)]
  }
  
  # Dealing with ungenotyped material
  if(any(!colnames(Z)%in%rownames(gen))){
    w = which(!colnames(Z)%in%rownames(gen))
    if(fixed){
      X = cbind(X,Z[,w])
      Z = Z[,-w]
      Weight = Weight[-w]
    }else{
      fixed = TRUE
      X = Z[,w]
      Z = Z[,-w]
      Weight = Weight[-w]
    }
  }
  
  # Getting the first round fitted
  # (a) Fixed
  fitSM = function(y,X,rdg=1){
    XX = crossprod(X)
    diag(XX) = diag(XX)+rdg
    Xy = crossprod(X,y)
    b = solve(XX,Xy)
    hat = tcrossprod(t(b),X)  
    e = y - as.vector(hat)
    ret = list('b'=as.vector(b),'e'=e,'y'=hat)
  }
  # intercept
  mu = mean(y)
  e = y-mu
  # other "fixed" effects with minimal ridge
  if(isTRUE(fixed)){
    LE = fitSM(e,X)
    e = LE$e
    b = LE$b
  }
  
  # (b) Preparing genotypic data
  gen0 = gen
  gen = gen[colnames(Z),]
  MAP = function(e,Z,Weight) as.vector(crossprod(Z,e))/Weight
  
  # Setting random forest parameters
  if(model=="RF"){
    RF = TRUE
    RFReg = function(x,y,NTR=50,MTR=50,SSZ=80,NSZ=4,
                     MNO=10,IMP=TRUE,RPL=TRUE)
          randomForest::randomForest(x=x,y=y,ntree=NTR,
                    mtry=MTR,sampsize=SSZ,nodesize=NSZ,
                    maxnodes=MNO,importance=IMP,
                    replace=RPL)
  }else{
    RF = FALSE
  }
  
  if(!RF){
    cat("Computing priors and genomic parameters \n")
    xx = apply(gen,2,crossprod)
    p = ncol(gen)
    g = rep(0,p)
    L = rep(1,p)
    Ve = 1
  }
  
  # (c) Running first round
  if(RF){
    E = MAP(e,Z,Weight)
    # FIT
    u=RFReg(gen,E,...)
    # PREDICT
    bv_pred = as.vector(predict(u,gen))
    bv = u$predicted
    bvs = as.vector(tcrossprod(bv,Z))
    e = e-bvs
    
  }else{
    
    # KERNEL AND REGRESSION METHODS
    d = rep(1,ncol(gen))
    E = MAP(e,Z,Weight) 
    update = KMUP(X=gen,b=g,d=d,xx=xx,e=E,L=L,Ve=Ve,pi=0)
    # First round of WGR: Setting priors
    if(KERN){
      R2 = 0.5
      df_prior = 5
      Sk_prior = R2 * var(y, na.rm = T) * (df_prior + 2)
    }else{
      R2 = 0.5
      df_prior = 3
      MSx = sum(apply(gen, 2, var, na.rm = T))
      S_prior = R2 * var(y, na.rm = T) * (df_prior + 2)/MSx
      shape_prior = 1.1
      rate_prior = (shape_prior - 1)/S_prior
    }
    Se_prior = (1-R2) * var(y, na.rm = T) * (df_prior + 2)
    
    # BV
    g = update$b
    bv = tcrossprod(g,gen)
    bvs = as.vector(tcrossprod(bv,Z))
    e = e-bvs
    
    if(KERN){
      Vb = (sum(g^2/V) + Sk_prior)/rchisq(1, df_prior + p)
      Ve = (crossprod(e)+Se_prior)/rchisq(1,n+df_prior)
      L = Ve/(Vb*V)
    }else{
      # Update BayesA variance components
      S_conj = rgamma(1, p * df_prior/2 + shape_prior, sum(1/(g^2))/2 + rate_prior)
      # For every round
      Vb = (S_conj + g^2)/rchisq(p, df_prior + 1)
      S_conj = rgamma(1, p * df_prior/2 + shape_prior,sum(1/Vb)/2 + rate_prior)  
      # Update Ve and Lambda
      Ve = (crossprod(e)+Se_prior)/rchisq(1,n+df_prior)
      L = Ve/Vb
    }
    
  }
  

  # (d) Adding Splines to X
  if(isTRUE(spline)){
    spEff = SPcov(R,e)
    Q = qr.solve(spEff,e)*spEff
    e = e - Q
    sp = Q
  }
  
  # What to STORE
  MCMC = seq(bi,it,th)
  mc = length(MCMC)
  if(RF){
    Mu = 0
    BV = rep(0,nrow(gen))
    BV_pred = rep(0,nrow(gen0))
  }else{
    G = VE = Mu = 0
    VB = rep(0,p)
  }
  
  if(isTRUE(spline)) SP = rep(0,n) 
  if(isTRUE(fixed)) B = rep(0,ncol(X))
  
  # LOOP STARTS HERE
  cat("Fitting the model \n")
  pb = txtProgressBar(style = 3)
  for(i in 1:it){
    
    # (a) Intercept
    muE = mean(e)
    mu = mu+muE
    e = e-muE
    
    # (b) Other fixed
    if(isTRUE(fixed)){
      LE = fitSM(e,X)
      b = LE$b + b
      e = LE$e
    }
    
    # (c) Update genetic components
    if(RF){
      
      bvs0 = bvs
      E = MAP(e,Z,Weight)+bv
      # FIT
      u=RFReg(gen,E,...)
      
      # PREDICT
      bv_pred = as.vector(predict(u,gen0))
      bv = u$predicted
      bvs = as.vector(tcrossprod(bv,Z))
      e = e+bvs0-bvs
      
    }else{
      
      bvs0 = bvs
      E = MAP(e,Z,Weight)
      update = KMUP(X=gen,b=g,d=d,xx=xx,e=E,L=L,Ve=Ve,pi=0)
      g = update$b
      bv = tcrossprod(g,gen)
      bvs = as.vector(tcrossprod(Z,bv))
      e = e+bvs0-bvs
      
      # (d) Update VC
      if(model=="BRR"){
        Va = (sum(g^2) + S_prior)/rchisq(1, df_prior + p)
        Vm = rep(Va,p)
        Ve = (crossprod(e)+Se_prior)/rchisq(1,n+df_prior)
        L = c(Ve/Vm)
      }
      if(model=="BayesA"){
        Vm = (S_conj + g^2)/rchisq(p, df_prior + 1)
        S_conj = rgamma(1, p * df_prior/2 + shape_prior,sum(1/Vb)/2 + rate_prior)  
        Ve = (crossprod(e)+Se_prior)/rchisq(1,n+df_prior)
        L = c(Ve/Vm)
      }
      if(KERN){
        Va = (sum(g^2/V) + Sk_prior)/rchisq(1, df_prior + p)
        Vm = rep(Va,p)
        Ve = (crossprod(e)+Se_prior)/rchisq(1,n+df_prior)
        L = Ve/(Vm*V)
      }
    }
    
    # (e) Update Splines
    if(isTRUE(spline)){
      spEff = SPcov(R,e)
      Q = qr.solve(spEff,e)*spEff
      e = e - Q
      sp = Q + sp
    }
    
    # (f) Storage
    if(i %in% MCMC){
      Mu = Mu+mu
      if(RF){
        BV = BV + bv
        BV_pred = BV_pred+bv_pred
      }else{
        G = G+g
        VB = VB+Vm
        VE = VE+Ve
      }
      if(isTRUE(fixed)) B = B+b
      if(isTRUE(spline)) SP = SP+sp
    }
    
    # (g) Update status bar
    setTxtProgressBar(pb, i/it)
  }
  close(pb)
  
  # Posterior means
  if(RF){
    BV = BV/mc
    BV_pred = BV_pred/mc
  }else{
    Mu = Mu/mc
    G = G/mc
    VB = VB/mc
    VE = VE/mc
    if(model=="BRR"|KERN) VB=mean(VB)
  }
  
  # Fitted
  if(RF){
    RND = as.vector(tcrossprod(BV,Z))
    names(RND) = rownames(gen)
    GEBV = BV_pred
    names(GEBV) = rownames(gen0)
  }else{
    RND = as.vector(tcrossprod(Z,tcrossprod(G,gen)))
    names(RND) = rownames(gen)
    GEBV = as.vector(tcrossprod(G,gen0))
    names(GEBV) = rownames(gen0)
  }
  
  # Hat
  hat = Mu + RND
  
  # Adding terms
  if(isTRUE(spline)){
    SP = SP/mc
    hat = hat + SP
  } 
  
  if(isTRUE(fixed)){
    B = B/mc; names(B)=colnames(X)
    FIX = as.vector(tcrossprod(B,X))
    hat = hat + FIX
  }
  
  if(RF){
    FINAL = list('hat'=hat,'obs'=y,
                 "mu"=Mu,"EBV"=GEBV)
  }else{
    FINAL = list('hat'=hat,'obs'=y,
                 "mu"=Mu,"Z"=Z,"g"=G,
                 "Vg"=VB,"Ve"=VE,
                 "EBV"=GEBV)
  }
  
  if(model=="BRR"|KERN){
    FINAL = c(FINAL,list('cxx'=mean(xx)))
  } 
  
  if(isTRUE(spline)){
    FINAL = c(FINAL,list("sp"=SP))
  }
  
  if(isTRUE(fixed)){
    FINAL = c(FINAL,list("X"=X,"b"=B))
  }
  
  class(FINAL) = "SSM"
  return(FINAL) 
  
}
