# Classical implementation
gwas = function(y,gen,fam=NULL,chr=NULL,window=NULL,fixed=FALSE){

  ##################
  ## INTRODUCTION ##
  ##################
  
  # REMOVAL OF MISSING Y's  
  anyNA = function(x) any(is.na(x))
  if(any(is.na(y))){
    w=which(is.na(y))
    y=y[-w];gen=gen[-w,];fam=fam[-w]
  }
  
  if(is.null(window)){method="RH"}else{method="EB"}
  SNPs=colnames(gen)
  
  # SETTING THE FIXED EFFECT, CHROMOSOME AND FAMILY WHEN IT IS NULL
  covariate=matrix(1,length(y),1)
  if(is.null(fam)){fam=rep(1,length(y))} else{if(any(is.na(fam))) stop("Family vector must have no NA's")} # calc
  if(is.null(chr)){chr=ncol(gen)} # calc
  
  # ORDERING DATA BY FAMILY
  
  cat("Ordering Data",'\n')
  FAM0 = as.vector(fam)
  organ=function(fam,y,covariate,Z){
    a=cbind(fam,y,covariate,Z);a=a[order(a[,1]),]
    b=array(summary(factor(fam)))
    c=c();for(i in 1:length(b)){d=rep(i,b[i]);c=c(c,d)}
    y=a[,2];fam=c;covariate=a[,3];Z=a[,-c(1:3)];return(list(fam,y,covariate,Z))}
  f=organ(fam,y,covariate,gen);fam=f[[1]];y=f[[2]];covariate=f[[3]];gen=f[[4]];rm(f);
  covariate=matrix(covariate,ncol=1)
  
  # MARKER IMPUT FUNCTION
  gen[gen==5]=NA
  IMPUT=function(X){
    X[is.na(X)]=5
    doEXP=function(X){
      Mode=function(x){x=na.omit(x);ux=unique(x);
      ux[which.max(tabulate(match(x,ux)))]}
      exp=Mode(X[,1])
      return(exp)}
    EXP=doEXP(X)
    N=length(X[1,])
    IMP=t(apply(X, MARGIN=1, FUN=inputRow, n = N, exp = EXP))
    return(IMP)}
  if(any(is.na(gen))) gen=IMPUT(gen) # calc
  
  # SPARSE DESING MATRIX FUNCTION
  cat("Generating Design Matrix",'\n')
  SDM = function(gen,fam){
    m = ncol(gen)
    n = nrow(gen)
    u = unique(fam)
    f = length(u)+1
    gg=matrix(0,n,m*f)
    gg[,(1:m-1)*f+1]=gen
    gen=-(gen-2)
    pb=txtProgressBar(style=3)
    for(i in 1:(f-1)){
      w = which(fam==i)
      gg[w,(1:m-1)*f+1+i]=gen[w,]
      setTxtProgressBar(pb,i/(f-1))
    };close(pb)
    return(t(gg))}
  GEN=SDM(gen,fam) # calc
  
  # GENOMIC RELATIONSHIP MATRIX
  Gmat=function(gen){
    ZZ=crossprod(gen) #individuals are columns
    lambda=sum(diag(ZZ))/nrow(ZZ)
    G=ZZ/lambda
    return(list(G,lambda))}
  
  INPUT=function(gen,fam,chr){
    if(sum(chr)!=ncol(gen)) stop("Number of markers and chromosomes don't match")
    if(length(fam)!=nrow(gen)) stop("Family and genotype don't match")
    MAP=function(gen,fam){ # MGD = Map of Genetic Distance
      ORDER=function(gen,fam){
        Matrix=cbind(fam,gen)
        Matrix=SORT(Matrix)
        mat=Matrix[,-1]
        return(mat)}
      SORT=function(foo){(foo[order(foo[,1]),])}
      # Orgenazing matrices
      Mat=data.frame(ORDER(gen,fam))
      Fam=SORT(matrix(fam,ncol=1))
      nF=dim(array((summary(factor(Fam)))))
      # Functions
      GD=function(snpA,snpB){ # genetic distance
        kosambi = function(r) min(.25*log((1+2*r)/(1-2*r)),.5)
        a0=which(snpA==0)
        a2=which(snpA==2)
        b0=which(snpB==0)
        b2=which(snpB==2)
        NR=length(intersect(a0,b0))+length(intersect(a2,b2))
        RE=length(intersect(a0,b2))+length(intersect(a2,b0))
        if(RE<NR){r=RE/(NR+RE)}else{r=NR/(NR+RE)}
        r = r/(2*(1-r)) # relation r=f(R) in RILs
        k = kosambi(r)
        return(k)}
      AR=function(SNP){ # able to recombine
        ar=function(snp){
          uni=unique(snp)
          int=intersect(uni,c(0,2))
          two=length(int)
          if(two<2) result=NA else result=TRUE
          return(result)}
        ar2=apply(SNP,MARGIN=2,FUN=ar)
        return(ar2)}
      # FUNCTION TO MAP DISTANCE BEFORE SPLITTING
      DOIT=function(Gen){
        nSNP=ncol(Gen)
        GDst=rep(0,nSNP)
        for(i in 2:nSNP){GDst[i]=GD(Gen[,(i-1)],Gen[,i])}
        Good=AR(Gen)
        Good2=c(1,Good[1:(length(Good)-1)])
        GDst=GDst*Good*Good2
        return(GDst)}
      # FUNCTION DOIT BY FAMILY
      DBF=function(gen,fam){ # DOIT by Family
        nF=dim(array((summary(factor(fam)))))
        LIST=split(gen,fam)
        LIST2=lapply(LIST,DOIT)
        A=matrix(0,nF,ncol(gen))
        for(i in 1:nF){A[i,]=unlist(LIST2[i])}
        Final=colMeans(A,na.rm=T)
        return(Final)}
      ccM=function(DBFoutput){
        vect=as.vector(DBFoutput)
        nSNP=length(vect)
        CGD=rep(0,nSNP) # Cumulative gen. dist.
        for(i in 2:nSNP){CGD[i]=vect[i]+CGD[i-1]}
        return(CGD)}
      A=DBF(Mat,Fam)
      A[which(is.na(A))]=0
      for(i in 1:length(A)){if(A[i]>0.5){A[i]=A[i]-0.5}}
      B=ccM(A)
      return(B)}
    # CHROMOSOME VECTOR
    CHROM=function(vector){ 
      len=length(vector);total=sum(vector);
      Vec1=c();
      for(i in 1:len){a=rep(i,vector[i]);Vec1=c(Vec1,a)};
      Vec2=c()
      for(i in 1:len){b=(1:vector[i]);Vec2=c(Vec2,b)};
      Final=cbind(Vec1,Vec2)
      colnames(Final)=c("chr","bin")
      return(Final)}
    a=CHROM(chr)
    ccM=suppressWarnings(round(MAP(gen,fam),3))
    begin=tapply(ccM,a[,1],function(x){X=x-x[1];return(X)},simplify = T)
    XX=c(); for(i in 1:length(chr)){XX=c(XX,unlist(begin[i]))}; begin=XX; rm(XX)
    end=tapply(begin,a[,1],function(x){X=x[-1];
    X=c(X,(X[length(X)]+0.5));return(X)},simplify = T)
    XX=c(); for(i in 1:length(chr)){XX=c(XX,unlist(end[i]))}; end=XX; rm(XX)
    size=end-begin
    final=cbind(a,begin,end,size,ccM)
    return(final)}
  MAP=INPUT(gen,fam,chr) # calc
  
  # Shizhong's MIXED function
  mixed<-function(x,y,kk){
    loglike<-function(theta){
      lambda<-exp(theta)
      logdt<-sum(log(lambda*delta+1))
      h<-1/(lambda*delta+1)
      yy<-sum(yu*h*yu)
      yx<-matrix(0,q,1)
      xx<-matrix(0,q,q)
      for(i in 1:q){
        yx[i]<-sum(yu*h*xu[,i])
        for(j in 1:q){
          xx[i,j]<-sum(xu[,i]*h*xu[,j])}}
      loglike<- -0.5*logdt-0.5*(n-q)*log(yy-t(yx)%*%solve(xx)%*%yx)-0.5*log(det(xx))
      return(-loglike)}
    fixed<-function(lambda){
      h<-1/(lambda*delta+1)
      yy<-sum(yu*h*yu)
      yx=timesVec(yu,h,xu,q)
      xx=timesMatrix(xu,h,xu,q,q)
      beta<-solve(xx,yx)
      sigma2<-(yy-t(yx)%*%solve(xx)%*%yx)/(n-q)
      var<-diag(solve(xx)*sigma2)
      stderr<-sqrt(var)
      return(c(beta,stderr,sigma2))}
    qq<-eigen(as.matrix(kk),symmetric=T)
    delta<-qq[[1]]
    uu<-qq[[2]]
    q<-ncol(x)
    n<-ncol(kk)
    vp<-var(y)
    yu<-t(uu)%*%y
    xu<-t(uu)%*%x
    theta<-0
    parm<-optim(par=theta,fn=loglike,hessian = TRUE,method="L-BFGS-B",lower=-5,upper=5)
    lambda<-exp(parm$par)
    conv<-parm$convergence
    fn1<-parm$value
    fn0<-loglike(-Inf)
    lrt<-2*(fn0-fn1)
    hess<-parm$hessian
    parmfix<-fixed(lambda)
    beta<-parmfix[1:q]
    stderr<-parmfix[(q+1):(2*q)]
    sigma2<-parmfix[2*q+1]
    lod<-lrt/4.61
    p_value<-1-pchisq(lrt,1)
    sigma2g<-lambda*sigma2
    goodness<-(vp-sigma2)/vp
    par<-data.frame(conv,fn1,fn0,lrt,goodness,beta,stderr,sigma2,lambda,sigma2g,lod,p_value)
    return(par)}
  
  # Shizhong's BLUP function
  blup<-function(gen,map,fam,x,y,kk,beta,lambda,cc){
    qq<-eigen(as.matrix(kk),symmetric=T)
    delta<-qq[[1]]
    uu<-qq[[2]]
    yu<-t(uu)%*%y
    xu<-t(uu)%*%x
    h<-1/(delta*lambda+1)
    r<- max(fam)+1
    m<-nrow(map)
    rr<-yu-xu%*%beta
    gamma<-matrix(0,m,r)
    for(k in 1:m){
      sub<-seq(((k-1)*r+1),((k-1)*r+r))
      z<-as.matrix(gen[sub,])
      zu<-t(uu)%*%t(z)
      zr<-matrix(0,r,1)
      for(i in 1:r){zr[i,1]=sum(rr*h*zu[,i])}
      gamma[k,]<-t(lambda*zr)/cc}
    return(gamma)}
  
  #############
  ## METHODS ##
  #############
  
  RANDOMsma = function (GEN=GEN,MAP=MAP,fam=fam,chr=chr,y=y,COV=covariate,SNPnames=SNPs){
    if(is.vector(y)!=TRUE) stop("Phenotypes: 'y' must be a vector")
    if(is.vector(chr)!=TRUE) stop("Chromosome: 'chr' must be a vector")
    if(is.vector(fam)!=TRUE) stop("Family: 'fam' must be a vector")
    if(is.matrix(gen)!=TRUE) stop("Genotypes: 'gen' must be a matrix")
    if(sum(chr)!=ncol(gen)) stop("Number of markers and chromosomes do not match")
    if(length(fam)!=nrow(gen)) stop("Family and genotype don't match")
    if(length(y)!=length(fam)) stop("Family and phenotypes must have the same length")
    if(length(y)!=nrow(gen)) stop("Dimensions of genotypes and phenotypes do not match")
    map=MAP
    gen=GEN
    cat("Calculating G matrix",'\n')
    G=Gmat(gen)
    kk=(t(G)[[1]]);cc=(t(G)[[2]])
    m<-nrow(map)
    r<-max(fam)+1
    s<-1
    ccM<-map[,6]
    n<-nrow(kk)
    x<-COV
    cat("Solving polygenic model",'\n')
    parm<-mixed(x,y,kk)
    lambda<-parm$lambda
    beta<-parm$beta
    cat("Starting Eigendecomposition",'\n')
    qq<-eigen(as.matrix(kk),symmetric=T)
    delta<-qq[[1]]
    uu<-qq[[2]]
    h<-1/(delta*lambda+1)
    x<-matrix(1,n,1)
    xu<-t(uu)%*%x
    yu<-t(uu)%*%y
    xx<-matrix(0,s,s)
    for(i in 1:s){for(j in 1:s){xx[i,j]<-sum(xu[,i]*h*xu[,j])}}
    loglike<-function(theta){
      xi<-exp(theta)
      tmp0<-zz*xi+diag(r)
      tmp<-xi*solve(tmp0)
      yHy<-yy-t(zy)%*%tmp%*%zy
      yHx<-yx-zx%*%tmp%*%zy
      xHx<-xx-zx%*%tmp%*%t(zx)
      logdt2<-log(det(tmp0))
      loglike<- -0.5*logdt2-0.5*(n-s)*log(yHy-t(yHx)%*%solve(xHx)%*%yHx)-0.5*log(det(xHx))
      return(-loglike)}
    fixed<-function(xi){
      tmp0<-zz*xi+diag(r)
      tmp<-xi*solve(tmp0)
      yHy<-yy-t(zy)%*%tmp%*%zy
      yHx<-yx-zx%*%tmp%*%zy
      xHx<-xx-zx%*%tmp%*%t(zx)
      zHy<-zy-zz%*%tmp%*%zy
      zHx<-zx-zx%*%tmp%*%zz
      zHz<-zz-zz%*%tmp%*%zz
      beta<-solve(xHx,yHx)
      tmp2<-solve(xHx)
      sigma2<-(yHy-t(yHx)%*%tmp2%*%yHx)/(n-s)
      gamma<-xi*zHy-xi*t(zHx)%*%tmp2%*%yHx
      var<-abs((xi*diag(r)-xi*zHz*xi)*as.numeric(sigma2))
      stderr<-sqrt(diag(var))
      result<-list(gamma,stderr,beta,sigma2)
      return(result)}
    parr<-numeric()
    blupp<-numeric()
    cat("Starting Marker Analysis",'\n')
    pb=txtProgressBar(style=3)
    for(k in 1:m){
      sub<-seq(((k-1)*r+1),((k-1)*r+r))
      genk<-gen[sub,]
      zu<-t(uu)%*%t(genk)
      yy<-sum(yu*h*yu)
      zu=as.matrix(zu)
      zx=timesMatrix(xu,h,zu,s,r)
      yx=timesVec(yu,h,xu,s)
      zy=timesVec(yu,h,zu,r)
      zz=timesMatrix(zu,h,zu,r,r)
      theta<-c(0)
      parm<-optim(par=theta,fn=loglike,hessian = TRUE,method="L-BFGS-B",lower=-5,upper=5)
      xi<-exp(parm$par)
      conv<-parm$convergence
      fn1<-parm$value
      fn0<-loglike(-Inf)
      hess<-parm$hessian
      parmfix<-fixed(xi)
      gamma<-parmfix[[1]][1:r]
      stderr<-parmfix[[2]][1:r]
      beta<-parmfix[[3]]
      sigma2<-parmfix[[4]]
      lam_k<-xi
      tau_k<-lam_k*sigma2
      lrt<-2*(fn0-fn1)
      if(lrt<0){
        lrt = 0
        lod = 0
        pval = 0
      }else{
        lod = lrt/4.61
        #pval = round(-log(dchisq(lrt,0.5),base = 10),2)
        pval = round(-log(pchisq(lrt,0.5,lower.tail=FALSE),base = 10),2)
        if(pval<0) pval = 0
      }
      if(r==2){
        names(gamma) = c("eff.A","eff.a")  
        #names(stderr) = c("sd.eff.A","sd.eff.a")
      }else{
        names(gamma) = c("std.eff",paste("eff.founder.",1:(r-1),sep=""))  
        #names(stderr) = c("sd.eff.std",paste("sd.eff.founder.",1:(r-1),sep=""))
      }
      gamma = t(data.frame(gamma))
      #stderr = t(data.frame(stderr))
      sigma2g = sigma2*(lambda+lam_k)
      h2 = sigma2g / (sigma2g+sigma2)
      par<-data.frame(conv,fn1,fn0,lod,pval,lrt,sigma2g,sigma2,h2,lam_k,"var.snp"=tau_k,"intercept"=beta,gamma)#,stderr)
      
      blup<-c(gamma)#,stderr)
      parr<-rbind(parr,par)
      blupp<-rbind(blupp,blup)
      setTxtProgressBar(pb,k/m)
    };close(pb)
    cat("Calculations were performed successfully",'\n')
    rownames(parr) = paste('SNP',1:m,sep='')
    return(list("PolyTest"=parr,"Method"="Empirical Bayes","MAP"=MAP,"SNPs"=SNPnames))}
  
  RandomCIM = function (GEN=GEN,MAP=MAP,fam=fam,chr=chr,y=y,COV=covariate,WIN=window,SNPnames=SNPs){
    if(is.vector(y)!=TRUE) stop("Phenotypes: 'y' must be a vector")
    if(is.vector(chr)!=TRUE) stop("Chromosome: 'chr' must be a vector")
    if(is.vector(fam)!=TRUE) stop("Family: 'fam' must be a vector")
    if(is.matrix(gen)!=TRUE) stop("Genotypes: 'gen' must be a matrix")
    if(sum(chr)!=ncol(gen)) stop("Number of markers and chromosomes do not match")
    if(length(fam)!=nrow(gen)) stop("Family and genotype don't match")
    if(length(y)!=length(fam)) stop("Family and phenotypes must have the same length")
    if(length(y)!=nrow(gen)) stop("Dimensions of genotypes and phenotypes do not match")
    gen=GEN
    map=MAP
    cat("Calculating G matrix",'\n')
    G=Gmat(gen)
    kk=(t(G)[[1]]);cc=(t(G)[[2]])
    m<-nrow(map)
    r<-max(fam)+1
    s<-1
    ccM<-map[,6]
    n<-nrow(kk)
    x<-COV
    cat("Solving polygenic model",'\n')
    parm<-mixed(x,y,kk)
    lambda<-parm$lambda
    beta<-parm$beta
    cat("Starting Eigendecomposition",'\n')
    qq<-eigen(as.matrix(kk),symmetric=T)
    delta<-qq[[1]]
    uu<-qq[[2]]
    h<-1/(delta*lambda+1)
    xu<-t(uu)%*%x
    yu<-t(uu)%*%y
    yy<-sum(yu*h*yu) 
    yx=timesVec(yu,h,xu,s)
    xx=timesMatrix(xu,h,xu,s,s)
    loglike<-function(theta){
      xi<-exp(theta)
      q<-length(xi)
      psi<-NULL
      for(i in 1:q){psi<-c(psi,rep(xi[i],r))}
      psi<-diag(psi)
      tmp0<-zz%*%psi+diag(r*q)
      tmp<-psi%*%solve(tmp0)
      yHy<-yy-t(zy)%*%tmp%*%zy
      yHx<-yx-zx%*%tmp%*%zy 
      xHx<-xx-zx%*%tmp%*%t(zx)
      logdt2<-log(det(tmp0))
      loglike<- -0.5*logdt2-0.5*(n-s)*log(yHy-t(yHx)%*%solve(xHx)%*%yHx)-0.5*log(det(xHx))
      return(-loglike)}
    psi<-NULL
    fixed<-function(xi){q<-length(xi);
    for(i in 1:q){psi<-c(psi,rep(xi[i],r))}
    psi<-diag(psi)
    tmp0<-zz%*%psi+diag(r*q)
    tmp<-psi%*%solve(tmp0)
    yHy<-yy-t(zy)%*%tmp%*%zy
    yHx<-yx-zx%*%tmp%*%zy
    xHx<-xx-zx%*%tmp%*%t(zx)
    zHy<-zy-zz%*%tmp%*%zy
    zHx<-zx-zx%*%tmp%*%zz
    zHz<-zz-zz%*%tmp%*%zz
    beta<-solve(xHx,yHx)
    tmp2<-solve(xHx)
    sigma2<-(yHy-t(yHx)%*%tmp2%*%yHx)/(n-s)
    gamma<-psi%*%zHy-psi%*%t(zHx)%*%tmp2%*%yHx
    var<-abs((psi%*%diag(r*q)-psi%*%zHz%*%psi)*as.numeric(sigma2))
    stderr<-sqrt(diag(var))
    result<-list(gamma,stderr,beta,sigma2)
    return(result)}
    parr<-numeric()
    blupp<-numeric()
    cat("Starting Marker Analysis",'\n')  
    pb=txtProgressBar(style=3)
    window<-WIN
    for(k in 1:m){
      sub<-seq(((k-1)*r+1),((k-1)*r+r))
      gen2<-gen[sub,]
      c1<-ccM[k]-0.5*window
      c3<-ccM[k]+0.5*window
      k1<-suppressWarnings(max(which(ccM<=c1)))
      k3<-suppressWarnings(min(which(ccM>=c3)))
      if(k1==-Inf){gen1<-gen[sub,sample(1:n)]
      } else {
        sub1<-seq(((k1-1)*r+1),((k1-1)*r+r))
        gen1<-gen[sub1,]}
      if(k3==Inf){gen3<-gen[sub,sample(1:n)]
      } else {
        sub3<-seq(((k3-1)*r+1),((k3-1)*r+r))
        gen3<-gen[sub3,]}
      gen1=as.matrix(gen1);gen2=as.matrix(gen2);gen3=as.matrix(gen3);
      ggg=t(rbind(gen1,gen2,gen3));
      zu<-t(uu)%*%ggg
      r3<-3*r
      zx=timesMatrix(xu,h,zu,s,r3)
      zy=timesVec(yu,h,zu,r3)
      zz=timesMatrix(zu,h,zu,r3,r3)
      theta<-c(0,0,0)
      parm<-optim(par=theta,fn=loglike,hessian = TRUE,method="L-BFGS-B",lower=-5,upper=5)
      xi<-exp(parm$par)
      conv<-parm$convergence
      fn1<-parm$value
      hess<-parm$hessian
      parmfix<-fixed(xi)
      gamma<-parmfix[[1]][(r+1):(2*r)]
      stderr<-parmfix[[2]][(r+1):(2*r)]
      beta<-parmfix[[3]]
      sigma2<-parmfix[[4]]
      lam_k<-xi[2]
      tau_k<-lam_k*sigma2
      zx<-matrix(zx[,-seq(r+1,2*r)],s,2*r)
      zy<-matrix(zy[-seq(r+1,2*r),],2*r,1)
      zz<-zz[-seq(r+1,2*r),-seq(r+1,2*r)]
      theta<-c(0,0)
      parm<-optim(par=theta,fn=loglike,hessian = TRUE,method="L-BFGS-B",lower=-5,upper=5)
      fn0<-parm$value
      lrt<-2*(fn0-fn1)
      if(lrt<0){
        lrt = 0
        lod = 0
        pval = 0
      }else{
        lod = lrt/4.61
        #pval = round(-log(dchisq(lrt,0.5),base = 10),2)
        pval = round(-log(pchisq(lrt,0.5,lower.tail=FALSE),base = 10),2)
        if(pval<0) pval = 0
      }    
      if(r==2){
        names(gamma) = c("eff.A","eff.a")  
        #names(stderr) = c("sd.eff.A","sd.eff.a")
      }else{
        names(gamma) = c("std.eff",paste("eff.founder.",1:(r-1),sep=""))  
        #names(stderr) = c("sd.eff.std",paste("sd.eff.founder.",1:(r-1),sep=""))
      }
      gamma = t(data.frame(gamma))
      #stderr = t(data.frame(stderr))
      sigma2g = sigma2*(lambda+lam_k)
      h2 = sigma2g / (sigma2g+sigma2)
      par<-data.frame(conv,fn1,fn0,lod,pval,lrt,sigma2g,sigma2,h2,lam_k,"var.snp"=tau_k,"intercept"=beta,gamma)#,stderr)
      blup<-c(gamma)#,stderr)
      parr<-rbind(parr,par)
      blupp<-rbind(blupp,blup)
      setTxtProgressBar(pb,k/m)
    };close(pb)
    cat("Calculations were performed successfully",'\n')
    rownames(parr) = paste('SNP',1:m,sep='')
    return(list("PolyTest"=parr,"Method"="Empirical Bayes with moving-window strategy","MAP"=MAP,"SNPs"=SNPnames))}
  
  FIXEDsma = function (GEN=GEN,MAP=MAP,fam=fam,chr=chr,y=y,COV=covariate,SNPnames=SNPs){
    if(is.vector(y)!=TRUE) stop("Phenotypes: 'y' must be a vector")
    if(is.vector(chr)!=TRUE) stop("Chromosome: 'chr' must be a vector")
    if(is.vector(fam)!=TRUE) stop("Family: 'fam' must be a vector")
    if(is.matrix(gen)!=TRUE) stop("Genotypes: 'gen' must be a matrix")
    if(sum(chr)!=ncol(gen)) stop("Number of markers and chromosomes do not match")
    if(length(fam)!=nrow(gen)) stop("Family and genotype don't match")
    if(length(y)!=length(fam)) stop("Family and phenotypes must have the same length")
    if(length(y)!=nrow(gen)) stop("Dimensions of genotypes and phenotypes do not match")
    map=MAP
    gen=GEN
    cat("Calculating genomic kinship",'\n')
    G=Gmat(gen)
    kk=(t(G)[[1]]);cc=(t(G)[[2]])
    m<-nrow(map)
    r<-max(fam)+1
    q<-1
    ccM<-map[,6]
    n<-nrow(kk)
    x<-COV
    cat("Solving polygenic model (P3D)",'\n')
    parm<-mixed(x,y,kk)
    lambda<-parm$lambda
    beta<-parm$beta
    cat("Starting Eigendecomposition",'\n')
    qq<-eigen(as.matrix(kk),symmetric=T)
    delta<-qq[[1]]
    uu<-qq[[2]]
    h<-1/(delta*lambda+1)
    xu<-t(uu)%*%x
    yu<-t(uu)%*%y
    yy<-sum(yu*h*yu)
    yx = timesVec(yu,h,xu,q)
    xx = timesMatrix(xu,h,xu,q,q)
    beta<-solve(xx,yx)
    parr<-numeric()
    blupp<-numeric()
    cat("Starting Single Marker Analysis",'\n')
    pb=txtProgressBar(style=3)
    for(k in 1:m){
      sub<-seq(((k-1)*r+1),((k-1)*r+r))
      zu<-t(uu)%*%t(gen[sub,])
      zu = as.matrix(zu)
      zx = timesMatrix(xu,h,zu,q,r)
      zy = timesVec(yu,h,zu,r)
      zz = timesMatrix(zu,h,zu,r,r)
      diag(zz) = diag(zz)+1e-8
      zzi<-solve(zz)
      b<-solve(zz,zy)
      s2<-(yy-t(zy)%*%zzi%*%zy)/(n-r-q);
      v<-zzi*drop(s2);
      g<-b-mean(b)
      gamma<-g
      wald<-t(g)%*%solve(v)%*%g
      stderr<-sqrt(diag(v))
      sigma2<-s2
      if(r==2){
        rownames(gamma) = c("eff.A","eff.a")  
        #names(stderr) = c("sd.eff.A","sd.eff.a")
      }else{
        rownames(gamma) = c("std.eff",paste("eff.",1:(r-1),sep=""))  
        #names(stderr) = c("sd.eff.std",paste("sd.eff.founder.",1:(r-1),sep=""))
      }
      gamma = t(data.frame(gamma))
      #stderr = t(data.frame(stderr))
      par<-data.frame(wald,sigma2,gamma)#,stderr)
      parr<-rbind(parr,par)
      setTxtProgressBar(pb,k/m)
    };close(pb)
    rownames(parr) = paste('SNP',1:m,sep='')
    cat("Calculations were performed successfully",'\n')
    return(list("PolyTest"=parr,"Method"="P3D","MAP"=MAP,"SNPs"=SNPnames))
    
  }
  
  
  #########
  ## RUN ##
  #########
  
  if(fixed==TRUE){
    fit=FIXEDsma(GEN=GEN,MAP=MAP,fam=fam,chr=chr,y=y,COV=covariate,SNPnames=SNPs)
  }else{
    if(method=="RH"){fit=RANDOMsma(GEN=GEN,MAP=MAP,fam=fam,chr=chr,y=y,COV=covariate,SNPnames=SNPs)}
    if(method=="EB"){fit=RandomCIM(GEN=GEN,MAP=MAP,fam=fam,chr=chr,y=y,COV=covariate,WIN=window,SNPnames=SNPs)}
    neg = which(fit$PolyTest$lrt<0);fit$PolyTest$lrt[neg]=0
  }
  
  # Fixing output
  NullEff = which(colMeans(fit$PolyTest)==0)
  whEff = grep('eff.',names(fit$PolyTest))
  NullCol = intersect(NullEff,whEff)
  if(any(NullCol)) fit$PolyTest = fit$PolyTest[,-NullCol]
  if(any(whEff)){ whEff = grep('eff.',names(fit$PolyTest))
  names(fit$PolyTest)[whEff] = paste('eff.',unique(FAM0),sep='')} 
  #rownames(fit$PolyTest) = fit$SNPs

class(fit) <- "NAM"

return(fit)}

# Optimized implementation
gwas2 = function(y,gen,fam=NULL,chr=NULL,fixed=FALSE,EIG=NULL,cov=NULL){
  
  ##################
  ## INTRODUCTION ##
  ##################
  
  # SETTING THE FIXED EFFECT, CHROMOSOME, COVARIATE AND FAMILY WHEN IT IS NULL
  if(is.null(fam)){fam=rep(1,length(y))};
  if(is.null(chr)){chr=ncol(gen)}
  if(is.null(cov)){covariate=matrix(1,length(y),1)}else{y=y-mean(y,na.rm=T); covariate=matrix(cov,ncol=1)}
  
  # REMOVAL OF MISSING Y's  
  anyNA = function(x) any(is.na(x))
  if(any(is.na(y))){
    wMIS=which(is.na(y))
    y=y[-wMIS]
    gen=gen[-wMIS,]
    fam=fam[-wMIS]
    covariate=covariate[-wMIS,]
  }else{wMIS=NULL}
  
  method="RH"
  fx=fixed
  SNPs=colnames(gen)
  
  # ORDERING DATA BY FAMILY
  if(is.data.frame(gen)) gen=data.matrix(gen)
  FAM0 = as.vector(fam)
  if(mean(order(fam)==1:length(fam))!=1){
    cat("Ordering Data",'\n')
    W = order(fam)
    fam = fam[W]
    y = y[W]
    covariate = covariate[W]
    gen = gen[W,]
    if(!is.null(EIG)){EIG$vectors = EIG$vectors[W,]}
  }
  covariate=matrix(covariate,ncol=1)  
  
  # MARKER IMPUT FUNCTION
  gen[gen==5]=NA
  IMPUT=function(X){
    X[is.na(X)]=5
    doEXP=function(X){
      Mode=function(x){x=na.omit(x);ux=unique(x);
      ux[which.max(tabulate(match(x,ux)))]}
      exp=Mode(X[,1])
      return(exp)}
    EXP=doEXP(X)
    N=length(X[1,])
    IMP=t(apply(X, MARGIN=1, FUN=inputRow, n = N, exp = EXP))
    return(IMP)}
  if(any(is.na(gen))){gen=IMPUT(gen)} # calc
  
  # acceleration
  GGG=function(G,fam){
    f=max(fam)+1
    POP = function(gfa){
      J = rep(0,f);g = gfa[1];fa = gfa[2]
      if(g==2){J[1]=2}else{if(g==1)
      {J[1]=1;J[fa+1]=1}else{J[fa+1]=2}}
      return(J)}
    gfa = cbind(G,fam)
    return(unlist(apply(gfa,1,POP)))}
  
  if(any(fam!=1)){
    Gmat = function(gen,fam){
      # common parent linear kernel
      g1 = tcrossprod(gen)
      g1 = g1/mean(diag(g1))  
      # founder parent linear kernel
      g2 = ((gen-1)*-1)+1
      g2 = tcrossprod(g2)
      g2 = g2/mean(diag(g2))
      # adjusting intra-family relationship
      for(i in unique(fam)){
        nam = which(fam==i)
        g1[nam,nam]=g2[nam,nam]+g1[nam,nam]}
      # Final estimates
      lambda = mean(diag(g1))
      G = g1/lambda
      return(list(G,lambda))}
    
  }else{
    Gmat = function(gen,fam){
      # common parent linear kernel
      g1 = tcrossprod(gen)
      # Final estimates
      lambda = mean(diag(g1))
      G = g1/lambda
      return(list(G,lambda))}
  }
  
  INPUT=function(gen,fam,chr){
    if(sum(chr)!=ncol(gen)) stop("Number of markers and chromosomes don't match")
    if(length(fam)!=nrow(gen)) stop("Family and genotype don't match")
    MAP=function(gen,fam){ # MGD = Map of Genetic Distance
      ORDER=function(gen,fam){
        Matrix=cbind(fam,gen)
        Matrix=SORT(Matrix)
        mat=Matrix[,-1]
        return(mat)}
      SORT=function(foo){(foo[order(foo[,1]),])}
      # Organazing matrices
      Mat=data.frame(ORDER(gen,fam))
      Fam=SORT(matrix(fam,ncol=1))
      nF=dim(array((summary(factor(Fam)))))
      # Functions
      GD=function(snpA,snpB){ # genetic distance
        kosambi = function(r) min(.25*log((1+2*r)/(1-2*r)),.5)
        a0=which(snpA==0)
        a2=which(snpA==2)
        b0=which(snpB==0)
        b2=which(snpB==2)
        NR=length(intersect(a0,b0))+length(intersect(a2,b2))
        RE=length(intersect(a0,b2))+length(intersect(a2,b0))
        if(RE<NR){r=RE/(NR+RE)}else{r=NR/(NR+RE)}
        r = r/(2*(1-r)) # relation r=f(R) in RILs
        d = kosambi(r)
        return(d)}
      AR=function(SNP){ # able to recombine
        ar=function(snp){
          uni=unique(snp)
          int=intersect(uni,c(0,2))
          two=length(int)
          if(two<2) result=NA else result=TRUE
          return(result)}
        ar2=apply(SNP,MARGIN=2,FUN=ar)
        return(ar2)}
      # FUNCTION TO MAP DISTANCE BEFORE SPLITTING
      DOIT=function(Gen){
        nSNP=ncol(Gen)
        GDst=rep(0,nSNP)
        for(i in 2:nSNP){GDst[i]=GD(Gen[,(i-1)],Gen[,i])}
        Good=AR(Gen)
        Good2=c(1,Good[1:(length(Good)-1)])
        GDst=GDst*Good*Good2
        return(GDst)}
      # FUNCTION DOIT BY FAMILY
      DBF=function(gen,fam){ # DOIT by Family
        nF=dim(array((summary(factor(fam)))))
        LIST=split(gen,fam)
        LIST2=lapply(LIST,DOIT)
        A=matrix(0,nF,ncol(gen))
        for(i in 1:nF){A[i,]=unlist(LIST2[i])}
        Final=colMeans(A,na.rm=T)
        return(Final)}
      ccM=function(DBFoutput){
        vect=as.vector(DBFoutput)
        nSNP=length(vect)
        CGD=rep(0,nSNP) # Cumulative gen. dist.
        for(i in 2:nSNP){CGD[i]=vect[i]+CGD[i-1]}
        return(CGD)}
      A=DBF(Mat,Fam)
      A[which(is.na(A))]=0
      for(i in 1:length(A)){if(A[i]>0.5){A[i]=A[i]-0.5}}
      B=ccM(A)
      return(B)}
    # CHROMOSOME VECTOR
    CHROM=function(vector){ 
      len=length(vector);total=sum(vector);
      Vec1=c();
      for(i in 1:len){a=rep(i,vector[i]);Vec1=c(Vec1,a)};
      Vec2=c()
      for(i in 1:len){b=(1:vector[i]);Vec2=c(Vec2,b)};
      Final=cbind(Vec1,Vec2)
      colnames(Final)=c("chr","bin")
      return(Final)}
    a=CHROM(chr)
    ccM=round(MAP(gen,fam),4)
    begin=tapply(ccM,a[,1],function(x){X=x-x[1];return(X)},simplify = T)
    XX=c(); for(i in 1:length(chr)){XX=c(XX,unlist(begin[i]))}; begin=XX; rm(XX)
    end=tapply(begin,a[,1],function(x){X=x[-1];
    X=c(X,(X[length(X)]+0.5));return(X)},simplify = T)
    XX=c(); for(i in 1:length(chr)){XX=c(XX,unlist(end[i]))}; end=XX; rm(XX)
    size=end-begin
    final=cbind(a,begin,end,size,ccM)
    return(final)}
  MAP=INPUT(gen,fam,chr) # calc
  
  # Shizhong's MIXED function
  mixed<-function(x,y,kk){
    
    loglike<-function(theta){
      lambda<-exp(theta)
      logdt<-sum(log(lambda*delta+1))
      h<-1/(lambda*delta+1)
      yy<-sum(yu*h*yu)
      yx<-matrix(0,q,1)
      xx<-matrix(0,q,q)
      for(i in 1:q){
        yx[i]<-sum(yu*h*xu[,i])
        for(j in 1:q){
          xx[i,j]<-sum(xu[,i]*h*xu[,j])}}
      loglike<- -0.5*logdt-0.5*(n-q)*log(yy-t(yx)%*%solve(xx)%*%yx)-0.5*log(det(xx))
      return(-loglike)}
    
    fixed<-function(lambda){
      h<-1/(lambda*delta+1)
      yy<-sum(yu*h*yu)
      yx=timesVec(yu,h,xu,q)
      xx=timesMatrix(xu,h,xu,q,q)
      beta<-solve(xx,yx)
      sigma2<-(yy-t(yx)%*%solve(xx)%*%yx)/(n-q)
      var<-diag(solve(xx)*sigma2)
      stderr<-sqrt(var)
      return(c(beta,stderr,sigma2))}
    
    if(!is.null(EIG)){
      qq=EIG
      if(!is.null(wMIS)){
        len = length(qq$values)-length(wMIS)
        qq$values = qq$values[1:len]*(len/length(qq$values))
        qq$vectors = qq$vectors[-wMIS,]
        CE = which(order(abs(cor(qq$vectors,y)))<=len)
        qq$vectors = qq$vectors[,CE]
      }
    }else{
      qq=eigen(as.matrix(kk),symmetric=T)
    }
    
    delta<-qq[[1]]
    uu<-qq[[2]]
    q<-ncol(x)
    n<-ncol(uu)
    vp<-var(y)
    yu<-crossprod(uu,y)
    xu<-crossprod(uu,x)
    theta<-0
    parm<-optim(par=theta,fn=loglike,hessian = TRUE,method="L-BFGS-B",lower=-5,upper=5)
    lambda<-exp(parm$par)
    conv<-parm$convergence
    fn1<-parm$value
    fn0<-loglike(-Inf)
    lrt<-2*(fn0-fn1)
    hess<-parm$hessian
    parmfix<-fixed(lambda)
    beta<-parmfix[1:q]
    stderr<-parmfix[(q+1):(2*q)]
    sigma2<-parmfix[2*q+1]
    lod<-lrt/4.61
    p_value<-1-pchisq(lrt,1)
    sigma2g<-lambda*sigma2
    goodness<-(vp-sigma2)/vp
    par<-data.frame(conv,fn1,fn0,lrt,goodness,beta,stderr,sigma2,lambda,sigma2g,lod,p_value)
    return(par)}
  
  # Shizhong's BLUP function
  blup<-function(gen,map,fam,x,y,kk,beta,lambda,cc){
    if(!is.null(EIG)){
      qq=EIG
      if(!is.null(wMIS)){
        len = length(qq$values)-length(wMIS)
        qq$values = qq$values[1:len]*(len/length(qq$values))
        qq$vectors = qq$vectors[-wMIS,]
        CE = which(order(abs(cor(qq$vectors,y)))<=len)
        qq$vectors = qq$vectors[,CE]
      }
    }else{qq<-eigen(as.matrix(kk),symmetric=T)}
    delta<-qq[[1]]
    uu<-qq[[2]]
    yu<-crossprod(uu,y)
    xu<-crossprod(uu,x)
    h<-1/(delta*lambda+1)
    r<- max(fam)+1
    m<-nrow(map)
    rr<-yu-xu%*%beta
    gamma<-matrix(0,m,r)
    for(k in 1:m){
      sub<-seq(((k-1)*r+1),((k-1)*r+r))
      z<-as.matrix(gen[sub,])
      zu<-t(uu)%*%t(z)
      zr<-matrix(0,r,1)
      for(i in 1:r){zr[i,1]=sum(rr*h*zu[,i])}
      gamma[k,]<-t(lambda*zr)/cc}
    return(gamma)}
  
  #############
  ## METHODS ##
  #############
  
  RANDOMsma = function (GEN=gen,MAP=MAP,fam=fam,chr=chr,y=y,COV=covariate,SNPnames=SNPs){
    if(is.vector(y)!=TRUE) stop("Phenotypes: 'y' must be a vector")
    if(is.vector(chr)!=TRUE) stop("Chromosome: 'chr' must be a vector")
    if(is.vector(fam)!=TRUE) stop("Family: 'fam' must be a vector")
    if(is.matrix(gen)!=TRUE) stop("Genotypes: 'gen' must be a matrix")
    if(sum(chr)!=ncol(gen)) stop("Number of markers and chromosomes do not match")
    if(length(fam)!=nrow(gen)) stop("Family and genotype don't match")
    if(length(y)!=length(fam)) stop("Family and phenotypes must have the same length")
    if(length(y)!=nrow(gen)) stop("Dimensions of genotypes and phenotypes do not match")
    map=MAP
    gen=GEN
    
    if(is.null(EIG)){
      cat("Calculating genomic kinship\n")
      G=Gmat(gen,fam)
      kk=(t(G)[[1]])
      cc=(t(G)[[2]])
    }else{
      kk=1
      cc=2
    }
    
    m<-nrow(map)
    r<-max(fam)+1
    s<-1
    n<-length(y)
    x<-COV
    cat("Solving polygenic model",'\n')
    parm<-mixed(x,y,kk)
    lambda<-parm$lambda
    beta<-parm$beta
    if(!is.null(EIG)){
      qq=EIG
      if(!is.null(wMIS)){
        len = length(qq$values)-length(wMIS)
        qq$values = qq$values[1:len]*(len/length(qq$values))
        qq$vectors = qq$vectors[-wMIS,]
        CE = which(order(abs(cor(qq$vectors,y)))<=len)
        qq$vectors = qq$vectors[,CE]
      }
    }else{
      cat("Starting Eigendecomposition",'\n')
      qq<-eigen(as.matrix(kk),symmetric=T)}
    delta<-qq[[1]]
    uu<-qq[[2]]
    h<-1/(delta*lambda+1)
    yu<-crossprod(uu,y)
    xu<-crossprod(uu,x)
    xx<-matrix(0,s,s)
    for(i in 1:s){for(j in 1:s){xx[i,j]<-sum(xu[,i]*h*xu[,j])}}
    loglike<-function(theta){
      xi<-exp(theta)
      tmp0<-zz*xi+diag(r)
      tmp<-xi*solve(tmp0)
      yHy<-yy-t(zy)%*%tmp%*%zy
      yHx<-yx-zx%*%tmp%*%zy
      xHx<-xx-zx%*%tmp%*%t(zx)
      logdt2<-log(det(tmp0))
      loglike<- -0.5*logdt2-0.5*(n-s)*log(yHy-t(yHx)%*%solve(xHx)%*%yHx)-0.5*log(det(xHx))
      return(-loglike)}
    fixed<-function(xi){
      tmp0<-zz*xi+diag(r)
      tmp<-xi*solve(tmp0)
      yHy<-yy-t(zy)%*%tmp%*%zy
      yHx<-yx-zx%*%tmp%*%zy
      xHx<-xx-zx%*%tmp%*%t(zx)
      zHy<-zy-zz%*%tmp%*%zy
      zHx<-zx-zx%*%tmp%*%zz
      zHz<-zz-zz%*%tmp%*%zz
      beta<-solve(xHx,yHx)
      tmp2<-solve(xHx)
      sigma2<-(yHy-t(yHx)%*%tmp2%*%yHx)/(n-s)
      gamma<-xi*zHy-xi*t(zHx)%*%tmp2%*%yHx
      var<-abs((xi*diag(r)-xi*zHz*xi)*as.numeric(sigma2))
      stderr<-sqrt(diag(var))
      result<-list(gamma,stderr,beta,sigma2)
      return(result)}
    if(max(fam)==1) r=1
    parr<-matrix(NA,nrow = m,ncol = (12+r))
    cat("Starting Marker Analysis",'\n')
    pb=txtProgressBar(style=3)
    for(k in 1:m){
      if(r==2){
        genk<-t(gen[,k])
      }else{genk<-GGG(gen[,k],fam)}
      zu<-t(uu)%*%t(genk)
      yy<-sum(yu*h*yu)
      zu=as.matrix(zu)
      zx=timesMatrix(xu,h,zu,s,r)
      yx=timesVec(yu,h,xu,s)
      zy=timesVec(yu,h,zu,r)
      zz=timesMatrix(zu,h,zu,r,r)
      theta<-c(0)
      parm<-optim(par=theta,fn=loglike,hessian = TRUE,method="L-BFGS-B",lower=-5,upper=5)
      xi<-exp(parm$par)
      conv<-parm$convergence
      fn1<-parm$value
      fn0<-loglike(-Inf)
      hess<-parm$hessian
      parmfix<-fixed(xi)
      gamma<-parmfix[[1]][1:r]
      stderr<-parmfix[[2]][1:r]
      beta<-parmfix[[3]]
      sigma2<-parmfix[[4]]
      lam_k<-xi
      tau_k<-lam_k*sigma2
      lrt<- 2*(fn0-fn1)
      if(lrt<0){
        lrt = 0
        lod = 0
        pval = 0
      }else{
        lod = lrt/4.61
        #pval = round(-log(dchisq(lrt,0.5),base = 10),2)
        pval = round(-log(pchisq(lrt,0.5,lower.tail=FALSE),base = 10),2)
        if(pval<0) pval = 0
      }
      if(r>2) names(gamma) = c("std.eff",paste("eff.founder.",1:(r-1),sep=""))
      gamma = t(data.frame(gamma))
      sigma2g = sigma2*(lambda+lam_k)
      h2 = sigma2g / (sigma2g+sigma2)
      par<-data.frame(conv,fn1,fn0,lod,pval,lrt,sigma2g,sigma2,h2,lam_k,"var.snp"=tau_k,"intercept"=beta,gamma)
      parr[k,] = as.numeric(par)
      setTxtProgressBar(pb,k/m)
    }
    parr=data.frame(parr)
    names(parr) = names(par)
    close(pb)
    cat("Calculations were performed successfully",'\n')
    if(max(fam)==1){names(parr)[13] = c("eff")}
    rownames(parr) = paste('SNP',1:m,sep='')
    return(list("PolyTest"=parr,"Method"="Empirical Bayes","MAP"=MAP,"SNPs"=SNPnames))}
  
  FIXEDsma = function (GEN=gen,MAP=MAP,fam=fam,chr=chr,y=y,COV=covariate,SNPnames=SNPs){
    if(is.vector(y)!=TRUE) stop("Phenotypes: 'y' must be a vector")
    if(is.vector(chr)!=TRUE) stop("Chromosome: 'chr' must be a vector")
    if(is.vector(fam)!=TRUE) stop("Family: 'fam' must be a vector")
    if(is.matrix(gen)!=TRUE) stop("Genotypes: 'gen' must be a matrix")
    if(sum(chr)!=ncol(gen)) stop("Number of markers and chromosomes do not match")
    if(length(fam)!=nrow(gen)) stop("Family and genotype don't match")
    if(length(y)!=length(fam)) stop("Family and phenotypes must have the same length")
    if(length(y)!=nrow(gen)) stop("Dimensions of genotypes and phenotypes do not match")
    map=MAP
    gen=GEN
    
    if(is.null(EIG)){
      cat("Calculating genomic kinship \n")
      G=Gmat(gen,fam)
      kk=(t(G)[[1]])
      cc=(t(G)[[2]])
    }else{
      kk=1
      cc=2
    }
    
    m<-nrow(map)
    r<-max(fam)+1
    q<-1
    ccM<-map[,6]
    n<-length(y)
    x<-COV
    cat("Solving polygenic model",'\n')
    parm<-mixed(x,y,kk)
    lambda<-parm$lambda
    beta<-parm$beta
    if(!is.null(EIG)){
      qq=EIG
      if(!is.null(wMIS)){
        len = length(qq$values)-length(wMIS)
        qq$values = qq$values[1:len]*(len/length(qq$values))
        qq$vectors = qq$vectors[-wMIS,]
        CE = which(order(abs(cor(qq$vectors,y)))<=len)
        qq$vectors = qq$vectors[,CE]
      }
    }else{
      cat("Starting Eigendecomposition",'\n')
      qq<-eigen(as.matrix(kk),symmetric=T)}
    delta<-qq[[1]]
    uu<-qq[[2]]
    h<-1/(delta*lambda+1)
    yu<-crossprod(uu,y)
    xu<-crossprod(uu,x)
    yy<-sum(yu*h*yu)
    yx = timesVec(yu,h,xu,q)
    xx = timesMatrix(xu,h,xu,q,q)
    beta<-solve(xx,yx)
    if(max(fam)==1) r=1
    parr<-matrix(NA,nrow=m,ncol=(2+r))
    cat("Starting Single Marker Analysis",'\n')
    pb=txtProgressBar(style=3)
    for(k in 1:m){
      if(r==1){
        genk<-t(gen[,k])
      }else{
        genk<-GGG(gen[,k],fam)
        }
      zu<-t(uu)%*%t(genk)
      zu = as.matrix(zu)
      zx = timesMatrix(xu,h,zu,q,r)
      zy = timesVec(yu,h,zu,r)
      zz = timesMatrix(zu,h,zu,r,r)
      diag(zz) = diag(zz)+1e-12
      zzi<-solve(zz)
      b<-solve(zz,zy)
      s2<-(yy-t(zy)%*%zzi%*%zy)/(n-r-q)
      v<-zzi*drop(s2);
      if(r==1){
        g<-b; gamma<-g
        wald<-t(b)%*%chol2inv(v)%*%b
      }else{
        g<-b-mean(b); gamma<-g
        wald<-t(g)%*%chol2inv(v)%*%g}
      sigma2<-s2
      if(r>1){
        rownames(gamma) = c("std.eff",paste("eff.",1:(r-1),sep=""))
        gamma = t(data.frame(gamma))
      }   
      par = data.frame(wald,sigma2,gamma)
      parr[k,] = as.numeric(par)
      setTxtProgressBar(pb,k/m)
    }
    close(pb)
    parr = data.frame(parr)
    names(parr) = names(par)
    if(r==1) names(parr)[3] = c("eff")
    cat("Calculations were performed successfully",'\n')
    rownames(parr) = paste('SNP',1:m,sep='')
    return(list("PolyTest"=parr,"Method"="P3D","MAP"=MAP,"SNPs"=SNPnames))}
  
  #########
  ## RUN ##
  #########
  
  y = as.vector(y)
  if(fx==1){
    fit=FIXEDsma(GEN=gen,MAP=MAP,fam=fam,chr=chr,y=y,COV=covariate,SNPnames=SNPs)
  }else{
    fit=RANDOMsma(GEN=gen,MAP=MAP,fam=fam,chr=chr,y=y,COV=covariate,SNPnames=SNPs)
    moda=function(x){
      it=5;ny=length(x);k=ceiling(ny/2)-1; while(it>1){
        y=sort(x); inf=y[1:(ny-k)]; sup=y[(k+1):ny]
        diffs=sup-inf; i=min(which(diffs==min(diffs)))
        M=median(y[i:(i+k)]); it=it-1}; return(M)}
    neg = which(fit$PolyTest$lrt<0);fit$PolyTest$lrt[neg]=0
    mo = moda(fit$PolyTest$lrt)
    mos = which(fit$PolyTest$lrt==mo);fit$PolyTest$lrt[mo]=0
  }
  
  # Fixing output
  NullEff = which(colMeans(fit$PolyTest)==0)
  if(!fx){
  whEff = grep('eff.',names(fit$PolyTest))
    }else{whEff = c()}
  NullCol = intersect(NullEff,whEff)
  if(any(NullCol)) fit$PolyTest = fit$PolyTest[,-NullCol]
  if(any(whEff)){ whEff = grep('eff.',names(fit$PolyTest))
    names(fit$PolyTest)[whEff] = paste('eff.',unique(FAM0),sep='')} 
  rownames(fit$PolyTest) = fit$SNPs

  class(fit) <- "NAM"
  return(fit)}

# GxE meta-analysis implementation
gwasGE = function(Phe,gen,fam,chr=NULL,cov=NULL,ge=FALSE,ammi=1){
  
  if(is.matrix(Phe)!=TRUE) stop("Phenotypes: 'Phe' must be a matrix")
  if(nrow(Phe)<2) stop("matrix Phe must have at least 2 columns (environments)")
  if(is.vector(fam)!=TRUE) stop("Family: 'fam' must be a vector")
  if(is.matrix(gen)!=TRUE) stop("Genotypes: 'gen' must be a numeric matrix")
  if(length(fam)!=nrow(gen)) stop("Family and genotype don't match")
  if(nrow(Phe)!=length(fam)) stop("Family and phenotypes must have the same length")
  if(nrow(Phe)!=nrow(gen)) stop("Dimensions of genotypes and phenotypes do not match")
  
  gwas = function(y,gen,fam,chr=NULL,fixed=FALSE,EIG=NULL,cov=NULL){
    
    ##################
    ## INTRODUCTION ##
    ##################
    
    # Centralize when there is a cov
    if(!is.null(cov)){
      tm = tapply(y,fam,mean,na.rm=T)
      cnt = tm[as.character(fam)]
      y=y-cnt}
    
    # REMOVAL OF MISSING Y's  
    anyNA = function(x) any(is.na(x))
    if(any(is.na(y))){
      wMIS=which(is.na(y))
      y=y[-wMIS];gen=gen[-wMIS,];fam=fam[-wMIS];cov=cov[-wMIS]
      #if(!is.null(EIG)){EIG$vectors = EIG$vectors[-wMIS,]}
    }else{wMIS=NULL}
    
    method="RH"
    fx=fixed
    SNPs=colnames(gen)
    
    # SETTING THE FIXED EFFECT, CHROMOSOME AND FAMILY WHEN IT IS NULL
    if(is.null(cov)) covariate=matrix(1,length(y),1) else covariate=matrix(cov,ncol=1)
    if(is.null(chr)){chr=ncol(gen)} # calc
    
    # ORDERING DATA BY FAMILY
    
    cat("Ordering Data",'\n')
    organ=function(fam,y,covariate,Z){
      a=cbind(fam,y,covariate,Z)
      a=a[order(a[,1]),]
      b=array(summary(factor(fam)))
      c=c()
      for(i in 1:length(b)){d=rep(i,b[i]);c=c(c,d)}
      y=a[,2];fam=c;covariate=a[,3];Z=a[,-c(1:3)];return(list(fam,y,covariate,Z))}
    f=organ(fam,y,covariate,gen)
    fam=f[[1]];y=f[[2]];covariate=f[[3]];gen=f[[4]];rm(f);
    covariate=matrix(covariate,ncol=1)
    
    # MARKER IMPUT FUNCTION
    gen[gen==5]=NA
    IMPUT=function(X){
      X[is.na(X)]=5
      doEXP=function(X){
        Mode=function(x){x=na.omit(x);ux=unique(x);
        ux[which.max(tabulate(match(x,ux)))]}
        exp=Mode(X[,1])
        return(exp)}
      EXP=doEXP(X)
      N=length(X[1,])
      IMP=t(apply(X, MARGIN=1, FUN=inputRow, n = N, exp = EXP))
      return(IMP)}
    if(any(is.na(gen))){gen=IMPUT(gen)} # calc
    
    # acceleration
    GGG=function(G,fam){
      f=max(fam)+1
      POP = function(gfa){
        J = rep(0,f);g = gfa[1];fa = gfa[2]
        if(g==2){J[1]=2}else{if(g==1)
        {J[1]=1;J[fa+1]=1}else{J[fa+1]=2}}
        return(J)}
      gfa = cbind(G,fam)
      return(unlist(apply(gfa,1,POP)))}
    if(any(fam!=1)){
      Gmat = function(gen,fam){
        # common parent linear kernel
        g1 = tcrossprod(gen)
        g1 = g1/mean(diag(g1))  
        # founder parent linear kernel
        g2 = ((gen-1)*-1)+1
        g2 = tcrossprod(g2)
        g2 = g2/mean(diag(g2))
        # adjusting intra-family relationship
        for(i in unique(fam)){
          nam = which(fam==i)
          g1[nam,nam]=g2[nam,nam]+g1[nam,nam]}
        # Final estimates
        lambda = mean(diag(g1))
        G = g1/mean(diag(g1))
        return(list(G,lambda))}
    }else{
      Gmat = function(gen,fam){
        # common parent linear kernel
        g1 = tcrossprod(gen)
        g1 = g1/mean(diag(g1))  
        # Final estimates
        lambda = mean(diag(g1))
        G = g1/mean(diag(g1))
        return(list(G,lambda))}
    }
    
    INPUT=function(gen,fam,chr){
      if(sum(chr)!=ncol(gen)) stop("Number of markers and chromosomes don't match")
      if(length(fam)!=nrow(gen)) stop("Family and genotype don't match")
      MAP=function(gen,fam){ # MGD = Map of Genetic Distance
        ORDER=function(gen,fam){
          Matrix=cbind(fam,gen)
          Matrix=SORT(Matrix)
          mat=Matrix[,-1]
          return(mat)}
        SORT=function(foo){(foo[order(foo[,1]),])}
        # Orgenazing matrices
        Mat=data.frame(ORDER(gen,fam))
        Fam=SORT(matrix(fam,ncol=1))
        nF=dim(array((summary(factor(Fam)))))
        # Functions
        GD=function(snpA,snpB){ # genetic distance
          kosambi = function(r) min(.25*log((1+2*r)/(1-2*r)),.5)
          a0=which(snpA==0)
          a2=which(snpA==2)
          b0=which(snpB==0)
          b2=which(snpB==2)
          NR=length(intersect(a0,b0))+length(intersect(a2,b2))
          RE=length(intersect(a0,b2))+length(intersect(a2,b0))
          if(RE<NR){r=RE/(NR+RE)}else{r=NR/(NR+RE)}
          r = r/(2*(1-r)) # relation r=f(R) in RILs
          d = kosambi(r)
          return(d)}
        AR=function(SNP){ # able to recombine
          ar=function(snp){
            uni=unique(snp)
            int=intersect(uni,c(0,2))
            two=length(int)
            if(two<2) result=NA else result=TRUE
            return(result)}
          ar2=apply(SNP,MARGIN=2,FUN=ar)
          return(ar2)}
        # FUNCTION TO MAP DISTANCE BEFORE SPLITTING
        DOIT=function(Gen){
          nSNP=ncol(Gen)
          GDst=rep(0,nSNP)
          for(i in 2:nSNP){GDst[i]=GD(Gen[,(i-1)],Gen[,i])}
          Good=AR(Gen)
          Good2=c(1,Good[1:(length(Good)-1)])
          GDst=GDst*Good*Good2
          return(GDst)}
        # FUNCTION DOIT BY FAMILY
        DBF=function(gen,fam){ # DOIT by Family
          nF=dim(array((summary(factor(fam)))))
          LIST=split(gen,fam)
          LIST2=lapply(LIST,DOIT)
          A=matrix(0,nF,ncol(gen))
          for(i in 1:nF){A[i,]=unlist(LIST2[i])}
          Final=colMeans(A,na.rm=T)
          return(Final)}
        ccM=function(DBFoutput){
          vect=as.vector(DBFoutput)
          nSNP=length(vect)
          CGD=rep(0,nSNP) # Cumulative gen. dist.
          for(i in 2:nSNP){CGD[i]=vect[i]+CGD[i-1]}
          return(CGD)}
        A=DBF(Mat,Fam)
        A[which(is.na(A))]=0
        for(i in 1:length(A)){if(A[i]>0.5){A[i]=A[i]-0.5}}
        B=ccM(A)
        return(B)}
      # CHROMOSOME VECTOR
      CHROM=function(vector){ 
        len=length(vector);total=sum(vector);
        Vec1=c();
        for(i in 1:len){a=rep(i,vector[i]);Vec1=c(Vec1,a)};
        Vec2=c()
        for(i in 1:len){b=(1:vector[i]);Vec2=c(Vec2,b)};
        Final=cbind(Vec1,Vec2)
        colnames(Final)=c("chr","bin")
        return(Final)}
      a=CHROM(chr)
      ccM=round(MAP(gen,fam),3)
      begin=tapply(ccM,a[,1],function(x){X=x-x[1];return(X)},simplify = T)
      XX=c(); for(i in 1:length(chr)){XX=c(XX,unlist(begin[i]))}; begin=XX; rm(XX)
      end=tapply(begin,a[,1],function(x){X=x[-1];
      X=c(X,(X[length(X)]+0.5));return(X)},simplify = T)
      XX=c(); for(i in 1:length(chr)){XX=c(XX,unlist(end[i]))}; end=XX; rm(XX)
      size=end-begin
      final=cbind(a,begin,end,size,ccM)
      return(final)}
    MAP=INPUT(gen,fam,chr) # calc
    
    # Shizhong's MIXED function
    mixed<-function(x,y,kk){
      loglike<-function(theta){
        lambda<-exp(theta)
        logdt<-sum(log(lambda*delta+1))
        h<-1/(lambda*delta+1)
        yy<-sum(yu*h*yu)
        yx<-matrix(0,q,1)
        xx<-matrix(0,q,q)
        for(i in 1:q){
          yx[i]<-sum(yu*h*xu[,i])
          for(j in 1:q){
            xx[i,j]<-sum(xu[,i]*h*xu[,j])}}
        loglike<- -0.5*logdt-0.5*(n-q)*log(yy-t(yx)%*%solve(xx)%*%yx)-0.5*log(det(xx))
        return(-loglike)}
      fixed<-function(lambda){
        h<-1/(lambda*delta+1)
        yy<-sum(yu*h*yu)
        yx=timesVec(yu,h,xu,q)
        xx=timesMatrix(xu,h,xu,q,q)
        beta<-solve(xx,yx)
        sigma2<-(yy-t(yx)%*%solve(xx)%*%yx)/(n-q)
        var<-diag(solve(xx)*sigma2)
        stderr<-sqrt(var)
        return(c(beta,stderr,sigma2))}
      
      # NEW FEATURE
      if(!is.null(EIG)){
        qq=EIG
        if(!is.null(wMIS)){
          len = length(qq$values)-length(wMIS)
          qq$values = qq$values[1:len]*(len/length(qq$values))
          qq$vectors = qq$vectors[-wMIS,]
          CE = which(order(abs(cor(qq$vectors,y)))<=len)
          qq$vectors = qq$vectors[,CE]
        }
      }else{qq<-eigen(as.matrix(kk),symmetric=T)}
      
      delta<-qq[[1]]
      uu<-qq[[2]]
      q<-ncol(x)
      n<-ncol(kk)
      vp<-var(y)
      yu<-crossprod(uu,y)
      xu<-crossprod(uu,x)
      
      theta<-0
      parm<-optim(par=theta,fn=loglike,hessian = TRUE,method="L-BFGS-B",lower=-5,upper=5)
      lambda<-exp(parm$par)
      conv<-parm$convergence
      fn1<-parm$value
      fn0<-loglike(-Inf)
      lrt<-2*(fn0-fn1)
      hess<-parm$hessian
      parmfix<-fixed(lambda)
      beta<-parmfix[1:q]
      stderr<-parmfix[(q+1):(2*q)]
      sigma2<-parmfix[2*q+1]
      lod<-lrt/4.61
      p_value<-1-pchisq(lrt,1)
      sigma2g<-lambda*sigma2
      goodness<-(vp-sigma2)/vp
      par<-data.frame(conv,fn1,fn0,lrt,goodness,beta,stderr,sigma2,lambda,sigma2g,lod,p_value)
      return(par)}
    
    # Shizhong's BLUP function
    blup<-function(gen,map,fam,x,y,kk,beta,lambda,cc){
      if(!is.null(EIG)){
        qq=EIG
        if(!is.null(wMIS)){
          len = length(qq$values)-length(wMIS)
          qq$values = qq$values[1:len]*(len/length(qq$values))
          qq$vectors = qq$vectors[-wMIS,]
          CE = which(order(abs(cor(qq$vectors,y)))<=len)
          qq$vectors = qq$vectors[,CE]
        }
      }else{qq<-eigen(as.matrix(kk),symmetric=T)}
      delta<-qq[[1]]
      uu<-qq[[2]]
      yu<-crossprod(uu,y)
      xu<-crossprod(uu,x)
      
      h<-1/(delta*lambda+1)
      r<- max(fam)+1
      m<-nrow(map)
      rr<-yu-xu%*%beta
      gamma<-matrix(0,m,r)
      for(k in 1:m){
        sub<-seq(((k-1)*r+1),((k-1)*r+r))
        z<-as.matrix(gen[sub,])
        zu<-t(uu)%*%t(z)
        zr<-matrix(0,r,1)
        for(i in 1:r){zr[i,1]=sum(rr*h*zu[,i])}
        gamma[k,]<-t(lambda*zr)/cc}
      return(gamma)}
    
    #############
    ## METHODS ##
    #############
    
    RANDOMsma = function (GEN=gen,MAP=MAP,fam=fam,chr=chr,y=y,COV=covariate,SNPnames=SNPs){
      map=MAP
      gen=GEN
      cat("Calculating G matrix",'\n')
      G=Gmat(gen,fam)
      kk=(t(G)[[1]]);cc=(t(G)[[2]])
      m<-nrow(map)
      r<-max(fam)+1
      s<-1
      ccM<-map[,6]
      n<-nrow(kk)
      x<-COV
      cat("Solving polygenic model",'\n')
      parm<-mixed(x,y,kk)
      lambda<-parm$lambda
      beta<-parm$beta
      if(!is.null(EIG)){
        qq=EIG
        if(!is.null(wMIS)){
          len = length(qq$values)-length(wMIS)
          qq$values = qq$values[1:len]*(len/length(qq$values))
          qq$vectors = qq$vectors[-wMIS,]
          CE = which(order(abs(cor(qq$vectors,y)))<=len)
          qq$vectors = qq$vectors[,CE]
        }
      }else{
        cat("Starting Eigendecomposition",'\n')
        qq<-eigen(as.matrix(kk),symmetric=T)}
      delta<-qq[[1]]
      uu<-qq[[2]]
      h<-1/(delta*lambda+1)
      yu<-crossprod(uu,y)
      xu<-crossprod(uu,x)
      xx<-matrix(0,s,s)
      for(i in 1:s){for(j in 1:s){xx[i,j]<-sum(xu[,i]*h*xu[,j])}}
      loglike<-function(theta){
        xi<-exp(theta)
        tmp0<-zz*xi+diag(r)
        tmp<-xi*solve(tmp0)
        yHy<-yy-t(zy)%*%tmp%*%zy
        yHx<-yx-zx%*%tmp%*%zy
        xHx<-xx-zx%*%tmp%*%t(zx)
        logdt2<-log(det(tmp0))
        loglike<- -0.5*logdt2-0.5*(n-s)*log(yHy-t(yHx)%*%solve(xHx)%*%yHx)-0.5*log(det(xHx))
        return(-loglike)}
      fixed<-function(xi){
        tmp0<-zz*xi+diag(r)
        tmp<-xi*solve(tmp0)
        yHy<-yy-t(zy)%*%tmp%*%zy
        yHx<-yx-zx%*%tmp%*%zy
        xHx<-xx-zx%*%tmp%*%t(zx)
        zHy<-zy-zz%*%tmp%*%zy
        zHx<-zx-zx%*%tmp%*%zz
        zHz<-zz-zz%*%tmp%*%zz
        beta<-solve(xHx,yHx)
        tmp2<-solve(xHx)
        sigma2<-(yHy-t(yHx)%*%tmp2%*%yHx)/(n-s)
        gamma<-xi*zHy-xi*t(zHx)%*%tmp2%*%yHx
        var<-abs((xi*diag(r)-xi*zHz*xi)*as.numeric(sigma2))
        stderr<-sqrt(diag(var))
        result<-list(gamma,stderr,beta,sigma2,var)
        return(result)}
      if(max(fam)==1) r=1
      parr<-matrix(NA,nrow = m,ncol = (12+r))
      VAR = array(NA,c(r,r,m))
      cat("Starting Marker Analysis",'\n')
      pb=txtProgressBar(style=3)
      for(k in 1:m){
        if(r==2){
          genk<-t(gen[,k])
        }else{genk<-GGG(gen[,k],fam)}
        zu<-t(uu)%*%t(genk)
        yy<-sum(yu*h*yu)
        zu=as.matrix(zu)
        zx=timesMatrix(xu,h,zu,s,r)
        yx=timesVec(yu,h,xu,s)
        zy=timesVec(yu,h,zu,r)
        zz=timesMatrix(zu,h,zu,r,r)
        theta<-c(0)
        parm<-optim(par=theta,fn=loglike,hessian = TRUE,method="L-BFGS-B",lower=-5,upper=5)
        xi<-exp(parm$par)
        conv<-parm$convergence
        fn1<-parm$value
        fn0<-loglike(-Inf)
        hess<-parm$hessian
        parmfix<-fixed(xi)
        gamma<-parmfix[[1]][1:r]
        stderr<-parmfix[[2]][1:r]
        beta<-parmfix[[3]]
        sigma2<-parmfix[[4]]
        VAR[,,k] = parmfix[[5]] # ERROR OF THE NEW MODEL
        lam_k<-xi
        tau_k<-lam_k*sigma2
        lrt<- 2*(fn0-fn1)
        if(lrt<0){
          lrt = 0
          lod = 0
          pval = 0
        }else{
          lod = lrt/4.61
          #pval = round(-log(dchisq(lrt,0.5),base = 10),2)
          pval = round(-log(pchisq(lrt,0.5,lower.tail=FALSE),base = 10),2)
          if(pval<0) pval = 0
        }
        if(r>2) names(gamma) = c("eff.std",paste("eff.founder.",1:(r-1),sep=""))
        gamma = t(data.frame(gamma))
        sigma2g = sigma2*(lambda+lam_k)
        h2 = sigma2g / (sigma2g+sigma2)
        par<-data.frame(conv,fn1,fn0,lod,pval,lrt,sigma2g,sigma2,h2,lam_k,"var.snp"=tau_k,"intercept"=beta,gamma)
        parr[k,] = as.numeric(par)
        setTxtProgressBar(pb,k/m)
      }
      parr=data.frame(parr)
      names(parr) = names(par)
      close(pb)
      cat("Calculations were performed successfully",'\n')
      if(max(fam)==1){names(parr)[13] = c("eff")}
      rownames(parr) = paste('SNP',1:m,sep='')
      return(list("PolyTest"=parr,"Method"="Empirical Bayes","MAP"=MAP,"SNPs"=SNPnames,"VAR"=VAR))}
    
    #########
    ## RUN ##
    #########
    
    
    fit=RANDOMsma(GEN=gen,MAP=MAP,fam=fam,chr=chr,y=y,COV=covariate,SNPnames=SNPs)
    moda=function(x){
      it=5;ny=length(x);k=ceiling(ny/2)-1; while(it>1){
        y=sort(x); inf=y[1:(ny-k)]; sup=y[(k+1):ny]
        diffs=sup-inf; i=min(which(diffs==min(diffs)))
        M=median(y[i:(i+k)]); it=it-1}; return(M)}
    neg = which(fit$PolyTest$lrt<0);fit$PolyTest$lrt[neg]=0
    mo = moda(fit$PolyTest$lrt)
    mos = which(fit$PolyTest$lrt==mo);fit$PolyTest$lrt[mo]=0
    
    class(fit) <- "NAM"
    
    return(fit)}
  
  ##### Meta-analysis FUNCTION STARTS HERE, WITH Phe, gen, fam, chr
  
  ByEnv = list()
  
  eigX = function(gen,fam){
    if(anyNA(gen)){
      mis = apply(gen,2,function(x) mean(is.na(x)))
      gen = gen[,mis<.8]
      gen = markov(gen,ncol(gen))
    }
    gm=function(gen,fam){
      g1 = tcrossprod(gen)
      g1 = g1/mean(diag(g1))
      g2 = ((gen-1)*-1)+1
      g2 = tcrossprod(g2)
      g2 = g2/mean(diag(g2))
      for (i in unique(fam)){
        nam = which(fam == i)
        g1[nam, nam] = g2[nam, nam] + g1[nam, nam]
      }
      lambda = mean(diag(g1))
      G = g1/lambda
      return(G)}
    eig = eigen(gm(gen,fam),symmetric=TRUE)
    return(eig)
  }
  
  if(!any(is.na(Phe))){
    cat('\n Pre-step: Eigendecomposition of kinship\n')
    eG = eigX(gen,fam)
  } 
  E = ncol(Phe)
  for(i in 1:E){
    cat('\n Starting environment',i,'\n')
    if(!any(is.na(Phe))){
      ByEnv[[i]] = gwas(Phe[,i],gen,fam,chr,EIG=eG)
    }else{
      ByEnv[[i]] = gwas(Phe[,i],gen,fam,chr) 
    }
    names(ByEnv)[i] = paste('E',i,sep='')
    save(ByEnv,file='GE.gwas.RData')
  }
  
  ################################
  ###       META-ANALYSIS      ###
  ################################
  
  cat('\n Performing meta-analysis\n')
  
  META = function(ByEnv,Phe,fam,GEonly = TRUE,PLOT=FALSE,ammi=1){
    
    ML = function(par){ # Woodbury solution
      vg = par[1]; ve = par[2];  vi = par[3]
      U = cbind(z*vg,w*ve,u*vi)
      v = zz*(vg^2)+ww*(ve^2)+zw*(vi^2)+kk
      RU = crossprod(ikk,U)
      Capacitance = crossprod(U,RU)
      diag(Capacitance) = diag(Capacitance)+1
      detCAP = determinant.matrix(Capacitance)$modulus[1]
      detK = determinant.matrix(kk)$modulus[1]
      d = detCAP+detK # logDet
      Capacitance = chol2inv(Capacitance)
      vv = ikk - tcrossprod(RU,Capacitance)%*%crossprod(U,ikk)
      xpx = sum(vv)
      xpy = sum(crossprod(vv,y))
      mu = xpy/xpx
      r = y-mu
      like = -0.5*( d + sum(crossprod(vv,r)*r) )
      return(-like)
    }
    
    MLtest = function(par){ # Woodbury solution
      vg = par[1]; ve = par[2];  vi = par[3]; mu = par[4];
      U = cbind(z*vg,w*ve,u*vi)
      v = zz*(vg^2)+ww*(ve^2)+zw*(vi^2)+kk
      RU = crossprod(ikk,U)
      Capacitance = crossprod(U,RU)
      diag(Capacitance) = diag(Capacitance)+1
      detCAP = determinant.matrix(Capacitance)$modulus[1]
      detK = determinant.matrix(kk)$modulus[1]
      d = detCAP+detK # logDet
      Capacitance = chol2inv(Capacitance)
      vv = ikk - tcrossprod(RU,Capacitance)%*%crossprod(U,ikk)
      r = y-mu
      like = -0.5*( d + sum(crossprod(vv,r)*r) )
      return(-like)
    }
    
    fixed<-function(par){
      vg = par[1]; ve = par[2];  vi = par[3]
      U = cbind(z*vg,w*ve,u*vi)
      v = zz*(vg^2)+ww*(ve^2)+zw*(vi^2)+kk
      RU = crossprod(ikk,U)
      Capacitance = crossprod(U,RU)
      diag(Capacitance) = diag(Capacitance)+1
      Capacitance = chol2inv(Capacitance)
      vv = ikk - tcrossprod(RU,Capacitance)%*%crossprod(U,ikk)
      xpx = sum(vv)
      xpy = sum(crossprod(vv,y))
      mu = xpy/xpx
      vb<-1/xpx
      return(c(mu,vb))
    }
    
    ### TESTING MODELS
    
    ML_MM = function(vi){ # Woodbury solution
      U = u*sqrt(vi)
      v = zw*vi+kk
      RU = crossprod(ikk,U)
      Capacitance = crossprod(U,RU)
      diag(Capacitance) = diag(Capacitance)+1
      detCAP = determinant.matrix(Capacitance)$modulus[1]
      detK = determinant.matrix(kk)$modulus[1]
      d = detCAP+detK # logDet
      Capacitance = chol2inv(Capacitance)
      vv = ikk - tcrossprod(RU,Capacitance)%*%crossprod(U,ikk)
      xpx = sum(vv)
      xpy = sum(crossprod(vv,y))
      mu = xpy/xpx
      r = y-mu
      like = -0.5*( d + sum(crossprod(vv,r)*r) )
      return(-like)
    }
    
    MLtest_MM = function(par){ # Woodbury solution
      vi = par[1]; mu = par[2];
      U = u*sqrt(vi)
      v = zw*vi+kk
      RU = crossprod(ikk,U)
      Capacitance = crossprod(U,RU)
      diag(Capacitance) = diag(Capacitance)+1
      detCAP = determinant.matrix(Capacitance)$modulus[1]
      detK = determinant.matrix(kk)$modulus[1]
      d = detCAP+detK # logDet
      Capacitance = chol2inv(Capacitance)
      vv = ikk - tcrossprod(RU,Capacitance)%*%crossprod(U,ikk)
      r = y-mu
      like = -0.5*( d + sum(crossprod(vv,r)*r) )
      return(-like)
    }
    
    fixed_MM = function(vi){
      U = u*vi
      v = zw*(vi^2)+kk
      RU = crossprod(ikk,U)
      Capacitance = crossprod(U,RU)
      diag(Capacitance) = diag(Capacitance)+1
      Capacitance = chol2inv(Capacitance)
      vv = ikk - tcrossprod(RU,Capacitance)%*%crossprod(U,ikk)
      xpx = sum(vv)
      xpy = sum(crossprod(vv,y))
      mu = xpy/xpx
      vb<-1/xpx
      return(c(mu,vb))
    }
    
    m = length(ByEnv[[1]]$SNPs)
    E = length(ByEnv)
    
    # Dealing with environments with not all families
    FAM = fam*!is.na(Phe); FAM[FAM==0]=NA
    Vi = rep(0,E); for(j in 1:E) Vi[j] = dim(ByEnv[[j]]$VAR[,,1])[1]
    cVi = c(0,cumsum(Vi))
    dimK = sum(Vi)
    kk = ikk = matrix(0,dimK,dimK)
    p = list()
    levF=c('0',sort(unique(as.vector(FAM[!is.na(FAM)]))))
    for(j in 1:E){
      f = factor( FAM[,j] ,levF)
      p[[j]]= c(1,sort(unique(as.numeric(f[!is.na(f)]))))
    } 
    
    # Preparing inputs, output matrix and design matrices
    p = as.factor(unlist(p))
    e = c(); for(j in 1:E) e = c(e,rep(j,Vi[j]))
    e = factor(e)
    y = rep(NA,length(p))
    z = model.matrix(~p-1)
    w = model.matrix(~e-1)
    zz= tcrossprod(z)
    ww= tcrossprod(w)
    x = matrix(1,dimK,1)
    u = model.matrix(~p:e-1)
    zw= tcrossprod(u)
    
    # Starting GxE gwas
    ee = as.numeric(e)
    pp = as.numeric(p)
    maxP = length(levels(p))
    
    ammi2 = function(y,pc){
      yy = lm(y~p+e)$residuals
      ge = matrix(0,maxP,E)
      for(k in 1:length(y)) ge[pp[k],ee[k]] = yy[k]
      S = svd(ge,pc,pc)
      SVD = matrix(0,length(y),pc)
      for(o in 1:pc){
        SS = tcrossprod(S$u[,o],S$v[,o])
        for(k in 1:length(y)) SVD[k,o] = SS[pp[k],ee[k]]
      }
      return(SVD)}
    
    ammi3 = function(y,pc){
      yy = y
      ge = matrix(0,maxP,E)
      for(k in 1:length(y)) ge[pp[k],ee[k]] = yy[k]
      S = svd(ge,pc,pc)
      SVD = matrix(0,length(y),pc)
      for(o in 1:pc){
        SS = tcrossprod(S$u[,o],S$v[,o])
        for(k in 1:length(y)) SVD[k,o] = SS[pp[k],ee[k]]
      }
      return(SVD)}
    
    ################
    ###   LOOP   ###  
    ################
    
    if(GEonly){
      
      output = matrix(NA,m,6)
      colnames(output) = c('mu','vb','vi','l1','l0','lrt')
      
      pb=txtProgressBar(style=3)
      for(i in 1:m){
        # generate known residual covariance matrix
        for(j in 1:E){
          k1 = ByEnv[[j]]$VAR[,,i]
          k2 = chol2inv(k1)
          kk[(cVi[j]+1):cVi[j+1],(cVi[j]+1):cVi[j+1]]=k1
          ikk[(cVi[j]+1):cVi[j+1],(cVi[j]+1):cVi[j+1]]=k2
        }
        # collect the regression coefficient of marker effect
        for(j in 1:E) y[(cVi[j]+1):cVi[j+1]] = ByEnv[[j]]$PolyTest[i,-c(1:12)]
        y = as.matrix(as.numeric(y),ncol=1)
        u = ammi3(y,ammi)
        zw = tcrossprod(u)
        # output functions (ie. Shizhong's lrt file)
        l0 = -MLtest_MM(rep(0,2))
        fit = optim(par=1,fn=ML_MM,method="L-BFGS-B",lower=0,upper=1e8)
        beta = fixed_MM(fit$par)
        mu = beta[1]
        vb = beta[2]
        l1 = -MLtest_MM(c(fit$par,mu))
        lrt = 2*(l1-l0)

        output[i,]=c(mu,vb,fit$par,l1,l0,lrt)
        setTxtProgressBar(pb,i/m)
      }
      close(pb)
      
    }else{
      
      output = matrix(NA,m,8)
      colnames(output) = c('mu','vb','vg','ve','vi','l1','l0','lrt')
      
      pb=txtProgressBar(style=3)
      for(i in 1:m){
        # generate known residual covariance matrix
        for(j in 1:E){
          k1 = ByEnv[[j]]$VAR[,,i]
          k2 = chol2inv(k1)
          kk[(cVi[j]+1):cVi[j+1],(cVi[j]+1):cVi[j+1]]=k1
          ikk[(cVi[j]+1):cVi[j+1],(cVi[j]+1):cVi[j+1]]=k2
        } 
        # collect the regression coefficient of marker effect
        for(j in 1:E) y[(cVi[j]+1):cVi[j+1]] = ByEnv[[j]]$PolyTest[i,-c(1:12)]
        y = as.matrix(as.numeric(y),ncol=1)
        u = ammi2(y,ammi)
        zw = tcrossprod(u)
        # output functions (ie. Shizhong's lrt file)
        l0 = -MLtest(rep(0,4))
        fit = optim(par=rep(0,3),fn=ML,method="L-BFGS-B",lower=0,upper=1e8)
        beta = fixed(fit$par)
        mu = beta[1]
        vb = beta[2]
        l1 = -MLtest(c(fit$par,mu))
        lrt = 2*(l1-l0)
        output[i,]=c(mu,vb,fit$par,l1,l0,lrt)
        setTxtProgressBar(pb,i/m)
      }
      close(pb)
      
    }
    
    
    if(PLOT){
      plot(ByEnv[[1]],pch=20,ylim=c(0,20),alpha=0.05/m)
      for(a in 2:length(ByEnv)) {
        par(new=TRUE)
        plot(ByEnv[[a]],pch=20,ylim=c(0,20),alpha=0.05/m)
      }
    }
    
    return(output)
  }
  
  
  output = META(ByEnv=ByEnv,Phe=Phe,fam=fam,GEonly=ge)
  
  FINAL_OUTPUT = list('PolyTest'=data.frame(output),
                      'GwasByEnv'=ByEnv,
                      'MAP'=ByEnv[[1]]$MAP,
                      'SNPs'=ByEnv[[1]]$SNPs,
                      'Method'="MetaGWA")
  class(FINAL_OUTPUT) = 'NAM'
  
  return(FINAL_OUTPUT)
  
}

# Independent populations
gwas3 = function(y,gen,fam=NULL,chr=NULL,EIG=NULL,cov=NULL){
  
  ##################
  ## INTRODUCTION ##
  ##################
  
  method="RH"
  fx=fixed=NULL
  
  # SETTING THE FIXED EFFECT, CHROMOSOME, COVARIATE AND FAMILY WHEN IT IS NULL
  if(is.null(fam)){fam=rep(1,length(y))};
  if(is.null(chr)){chr=ncol(gen)}
  if(is.null(cov)){covariate=matrix(1,length(y),1)}else{y=y-mean(y,na.rm=T); covariate=matrix(cov,ncol=1)}
  FAM0 = as.vector(fam)
  
  # REMOVAL OF MISSING Y's  
  anyNA = function(x) any(is.na(x))
  if(anyNA(y)&!is.null(EIG)) stop("No missing allowed when EIG is provided")
  if(anyNA(y)){
    wMIS=which(is.na(y))
    y=y[-wMIS];gen=gen[-wMIS,];fam=fam[-wMIS];covariate=covariate[-wMIS]
  }else{
    wMIS=NULL}
  
  # MARKER QUALITY CONTROLS FUNCTIONS
  max75 = function(X,chr){
    nas = function(x) mean(is.na(x))
    m = apply(X,2,nas)
    if(any(m>0.75)){
      w = which(m>0.75)
      X = X[,-w]
      Chr = c()
      for(i in 1:length(chr)) Chr = c(Chr, rep(i,chr[i]))
      Chr = Chr[-w]
      chr = data.frame(table(Chr))[,2]}
    return(list("gen"=X,"chr"=chr))}
  IMPUTATION = function(X){
    dn = dimnames(X)
    imp = function(x){
      if(any(is.na(x))){
        x[is.na(x)] = mean(x,na.rm = TRUE)}
      return(x)}
    X = apply(X,2,imp)
    dimnames(X) = dn
    X = data.matrix(X)
    return(X)}
  # NA's?
  if(any(is.na(gen))){
    cat('Imputing missing loci \n')
    gen = max75(gen,chr)
    chr = gen$chr
    gen = IMPUTATION(gen$gen)
  }
  SNPs=colnames(gen)
  
  ## Acceleration GGG defines allele coding
  GGG=function(G,fam) t(model.matrix(~G:factor(fam)-1))
  
  if(any(fam!=1)){
    Gmat = function(gen,fam){
      # common parent linear kernel
      g1 = tcrossprod(gen)
      g1 = g1/mean(diag(g1))  
      # adjusting intra-family relationship
      for(i in unique(fam)){
        nam = which(fam==i)
        g1[nam,nam]=2*g1[nam,nam]}
      # Final estimates
      lambda = mean(diag(g1))
      G = g1/mean(diag(g1))
      return(list(G,lambda))}
  }else{
    Gmat = function(gen,fam){
      # common parent linear kernel
      g1 = tcrossprod(gen)
      g1 = g1/mean(diag(g1))  
      # Final estimates
      lambda = mean(diag(g1))
      G = g1/mean(diag(g1))
      return(list(G,lambda))}
  }
  
  INPUT=function(gen,fam,chr){
    if(sum(chr)!=ncol(gen)) stop("Number of markers and chromosomes don't match")
    if(length(fam)!=nrow(gen)) stop("Family and genotype don't match")
    MAP=function(gen,fam){ # MGD = Map of Genetic Distance
      ORDER=function(gen,fam){
        Matrix=cbind(fam,gen)
        Matrix=SORT(Matrix)
        mat=Matrix[,-1]
        return(mat)}
      SORT=function(foo){(foo[order(foo[,1]),])}
      # Orgenazing matrices
      Mat=data.frame(ORDER(gen,fam))
      Fam=SORT(matrix(fam,ncol=1))
      nF=dim(array((summary(factor(Fam)))))
      # Functions
      GD=function(snpA,snpB){ # genetic distance
        kosambi = function(r) min(.25*log((1+2*r)/(1-2*r)),.5);
        a0=which(snpA==min(snpA))
        a2=which(snpA==max(snpA))
        b0=which(snpB==min(snpB))
        b2=which(snpB==max(snpB))
        NR=length(intersect(a0,b0))+length(intersect(a2,b2))
        RE=length(intersect(a0,b2))+length(intersect(a2,b0))
        if(RE<NR){r=RE/(NR+RE)}else{r=NR/(NR+RE)}
        r = r/(2*(1-r)) # relation r=f(R) in RILs
        d = kosambi(r)
        return(d)}
      AR=function(SNP){ # able to recombine
        ar=function(snp){
          uni=unique(snp)
          int=intersect(uni,c(0,2))
          two=length(int)
          if(two<2) result=NA else result=TRUE
          return(result)}
        ar2=apply(SNP,MARGIN=2,FUN=ar)
        return(ar2)}
      # FUNCTION TO MAP DISTANCE BEFORE SPLITTING
      DOIT=function(Gen){
        nSNP=ncol(Gen)
        GDst=rep(0,nSNP)
        for(i in 2:nSNP){GDst[i]=GD(Gen[,(i-1)],Gen[,i])}
        Good=AR(Gen)
        Good2=c(1,Good[1:(length(Good)-1)])
        GDst=GDst*Good*Good2
        return(GDst)}
      # FUNCTION DOIT BY FAMILY
      DBF=function(gen,fam){ # DOIT by Family
        nF=dim(array((summary(factor(fam)))))
        LIST=split(gen,fam)
        LIST2=lapply(LIST,DOIT)
        A=matrix(0,nF,ncol(gen))
        for(i in 1:nF){A[i,]=unlist(LIST2[i])}
        Final=colMeans(A,na.rm=T)
        return(Final)}
      ccM=function(DBFoutput){
        vect=as.vector(DBFoutput)
        nSNP=length(vect)
        CGD=rep(0,nSNP) # Cumulative gen. dist.
        for(i in 2:nSNP){CGD[i]=vect[i]+CGD[i-1]}
        return(CGD)}
      A=DBF(Mat,Fam)
      A[which(is.na(A))]=0
      for(i in 1:length(A)){if(A[i]>0.5){A[i]=A[i]-0.5}}
      B=ccM(A)
      return(B)}
    # CHROMOSOME VECTOR
    CHROM=function(vector){ 
      len=length(vector);total=sum(vector);
      Vec1=c();
      for(i in 1:len){a=rep(i,vector[i]);Vec1=c(Vec1,a)};
      Vec2=c()
      for(i in 1:len){b=(1:vector[i]);Vec2=c(Vec2,b)};
      Final=cbind(Vec1,Vec2)
      colnames(Final)=c("chr","bin")
      return(Final)}
    a=CHROM(chr)
    ccM=round(MAP(gen,fam),3)
    begin=tapply(ccM,a[,1],function(x){X=x-x[1];return(X)},simplify = T)
    XX=c(); for(i in 1:length(chr)){XX=c(XX,unlist(begin[i]))}; begin=XX; rm(XX)
    end=tapply(begin,a[,1],function(x){X=x[-1];
    X=c(X,(X[length(X)]+0.5));return(X)},simplify = T)
    XX=c(); for(i in 1:length(chr)){XX=c(XX,unlist(end[i]))}; end=XX; rm(XX)
    size=end-begin
    final=cbind(a,begin,end,size,ccM)
    return(final)}
  MAP=INPUT(gen,fam,chr) # calc
  
  # Shizhong's MIXED function
  mixed<-function(x,y,kk){
    loglike<-function(theta){
      lambda<-exp(theta)
      logdt<-sum(log(lambda*delta+1))
      h<-1/(lambda*delta+1)
      yy<-sum(yu*h*yu)
      yx<-matrix(0,q,1)
      xx<-matrix(0,q,q)
      for(i in 1:q){
        yx[i]<-sum(yu*h*xu[,i])
        for(j in 1:q){
          xx[i,j]<-sum(xu[,i]*h*xu[,j])}}
      loglike<- -0.5*logdt-0.5*(n-q)*log(yy-t(yx)%*%solve(xx)%*%yx)-0.5*log(det(xx))
      return(-loglike)}
    fixed<-function(lambda){
      h<-1/(lambda*delta+1)
      yy<-sum(yu*h*yu)
      yx=timesVec(yu,h,xu,q)
      xx=timesMatrix(xu,h,xu,q,q)
      beta<-solve(xx,yx)
      sigma2<-(yy-t(yx)%*%solve(xx)%*%yx)/(n-q)
      var<-diag(solve(xx)*sigma2)
      stderr<-sqrt(var)
      return(c(beta,stderr,sigma2))}
    
    # NEW FEATURE
    if(!is.null(EIG)){
      qq=EIG
      if(!is.null(wMIS)){
        len = length(qq$values)-length(wMIS)
        qq$values = qq$values[1:len]*(len/length(qq$values))
        qq$vectors = qq$vectors[-wMIS,]
        CE = which(order(abs(cor(qq$vectors,y)))<=len)
        qq$vectors = qq$vectors[,CE]
      }
    }else{qq<-eigen(as.matrix(kk),symmetric=T)}
    
    delta<-qq[[1]]
    uu<-qq[[2]]
    q<-ncol(x)
    n<-ncol(kk)
    vp<-var(y)
    yu<-crossprod(uu,y)
    xu<-crossprod(uu,x)
    
    theta<-0
    parm<-optim(par=theta,fn=loglike,hessian = TRUE,method="L-BFGS-B",lower=-5,upper=5)
    lambda<-exp(parm$par)
    conv<-parm$convergence
    fn1<-parm$value
    fn0<-loglike(-Inf)
    lrt<-2*(fn0-fn1)
    hess<-parm$hessian
    parmfix<-fixed(lambda)
    beta<-parmfix[1:q]
    stderr<-parmfix[(q+1):(2*q)]
    sigma2<-parmfix[2*q+1]
    lod<-lrt/4.61
    p_value<-1-pchisq(lrt,1)
    sigma2g<-lambda*sigma2
    goodness<-(vp-sigma2)/vp
    par<-data.frame(conv,fn1,fn0,lrt,goodness,beta,stderr,sigma2,lambda,sigma2g,lod,p_value)
    return(par)}
  
  # Shizhong's BLUP function
  blup<-function(gen,map,fam,x,y,kk,beta,lambda,cc){
    if(!is.null(EIG)){
      qq=EIG
      if(!is.null(wMIS)){
        len = length(qq$values)-length(wMIS)
        qq$values = qq$values[1:len]*(len/length(qq$values))
        qq$vectors = qq$vectors[-wMIS,]
        CE = which(order(abs(cor(qq$vectors,y)))<=len)
        qq$vectors = qq$vectors[,CE]
      }
    }else{qq<-eigen(as.matrix(kk),symmetric=T)}
    delta<-qq[[1]]
    uu<-qq[[2]]
    yu<-crossprod(uu,y)
    xu<-crossprod(uu,x)
    
    h<-1/(delta*lambda+1)
    r<- max(fam)+1
    m<-nrow(map)
    rr<-yu-xu%*%beta
    gamma<-matrix(0,m,r)
    for(k in 1:m){
      sub<-seq(((k-1)*r+1),((k-1)*r+r))
      z<-as.matrix(gen[sub,])
      zu<-t(uu)%*%t(z)
      zr<-matrix(0,r,1)
      for(i in 1:r){zr[i,1]=sum(rr*h*zu[,i])}
      gamma[k,]<-t(lambda*zr)/cc}
    return(gamma)}
  
  #############
  ## METHODS ##
  #############
  
  RANDOMsma = function (gen=gen,MAP=MAP,fam=fam,chr=chr,y=y,COV=covariate,SNPnames=SNPs){
    # GEN=gen;COV=covariate;SNPnames=SNPs
    map=MAP
    cat("Calculating G matrix",'\n')
    G=Gmat(gen,fam)
    kk=(t(G)[[1]]);cc=(t(G)[[2]])
    m<-nrow(map)
    r<-max(fam)
    s<-1
    ccM<-map[,6]
    n<-nrow(kk)
    x<-COV
    cat("Solving polygenic model",'\n')
    parm<-mixed(x,y,kk)
    lambda<-parm$lambda
    beta<-parm$beta
    if(!is.null(EIG)){
      qq=EIG
      if(!is.null(wMIS)){
        len = length(qq$values)-length(wMIS)
        qq$values = qq$values[1:len]*(len/length(qq$values))
        qq$vectors = qq$vectors[-wMIS,]
        CE = which(order(abs(cor(qq$vectors,y)))<=len)
        qq$vectors = qq$vectors[,CE]
      }
    }else{
      cat("Starting Eigendecomposition",'\n')
      qq<-eigen(as.matrix(kk),symmetric=T)}
    delta<-qq[[1]]
    uu<-qq[[2]]
    h<-1/(delta*lambda+1)
    yu<-crossprod(uu,y)
    xu<-crossprod(uu,x)
    xx<-matrix(0,s,s)
    for(i in 1:s){for(j in 1:s){xx[i,j]<-sum(xu[,i]*h*xu[,j])}}
    loglike<-function(theta){
      xi<-exp(theta)
      tmp0<-zz*xi+diag(r)
      tmp<-xi*solve(tmp0)
      yHy<-yy-t(zy)%*%tmp%*%zy
      yHx<-yx-zx%*%tmp%*%zy
      xHx<-xx-zx%*%tmp%*%t(zx)
      logdt2<-log(det(tmp0))
      loglike<- -0.5*logdt2-0.5*(n-s)*log(yHy-t(yHx)%*%solve(xHx)%*%yHx)-0.5*log(det(xHx))
      return(-loglike)}
    fixed<-function(xi){
      tmp0<-zz*xi+diag(r)
      tmp<-xi*solve(tmp0)
      yHy<-yy-t(zy)%*%tmp%*%zy
      yHx<-yx-zx%*%tmp%*%zy
      xHx<-xx-zx%*%tmp%*%t(zx)
      zHy<-zy-zz%*%tmp%*%zy
      zHx<-zx-zx%*%tmp%*%zz
      zHz<-zz-zz%*%tmp%*%zz
      beta<-solve(xHx,yHx)
      tmp2<-solve(xHx)
      sigma2<-(yHy-t(yHx)%*%tmp2%*%yHx)/(n-s)
      gamma<-xi*zHy-xi*t(zHx)%*%tmp2%*%yHx
      var<-abs((xi*diag(r)-xi*zHz*xi)*as.numeric(sigma2))
      stderr<-sqrt(diag(var))
      result<-list(gamma,stderr,beta,sigma2,var)
      return(result)}
    parr<-matrix(NA,nrow = m,ncol = (12+r))
    VAR = array(NA,c(r,r,m))
    cat("Starting Marker Analysis",'\n')
    pb=txtProgressBar(style=3)
    for(k in 1:m){
      
      if(r==1){
        genk <- t(gen[,k])
      }else{
        genk <- GGG(gen[,k],fam)
      }
      
      zu<-t(uu)%*%t(genk)
      yy<-sum(yu*h*yu)
      zu=as.matrix(zu)
      zx=timesMatrix(xu,h,zu,s,r)
      yx=timesVec(yu,h,xu,s)
      zy=timesVec(yu,h,zu,r)
      zz=timesMatrix(zu,h,zu,r,r)
      theta<-c(0)
      parm<-optim(par=theta,fn=loglike,hessian = TRUE,method="L-BFGS-B",lower=-5,upper=5)
      xi<-exp(parm$par)
      conv<-parm$convergence
      fn1<-parm$value
      fn0<-loglike(-Inf)
      hess<-parm$hessian
      parmfix<-fixed(xi)
      gamma<-parmfix[[1]][1:r]
      stderr<-parmfix[[2]][1:r]
      beta<-parmfix[[3]]
      sigma2<-parmfix[[4]]
      VAR[,,k] = parmfix[[5]] # ERROR OF THE NEW MODEL
      lam_k<-xi
      tau_k<-lam_k*sigma2
      lrt<- 2*(fn0-fn1)
      if(lrt<0){
        lrt = 0
        lod = 0
        pval = 0
      }else{
        lod = lrt/4.61
        #pval = round(-log(dchisq(lrt,0.5),base = 10),2)
        pval = round(-log(pchisq(lrt,0.5,lower.tail=FALSE),base = 10),2)
        if(pval<0) pval = 0
      }
      if(r>2) names(gamma) = paste("eff.",1:r,sep="")
      gamma = t(data.frame(gamma))
      sigma2g = sigma2*(lambda+lam_k)
      h2 = sigma2g / (sigma2g+sigma2)
      par<-data.frame(conv,fn1,fn0,lod,pval,lrt,sigma2g,sigma2,h2,lam_k,"var.snp"=tau_k,"intercept"=beta,gamma)
      parr[k,] = as.numeric(par)
      setTxtProgressBar(pb,k/m)
    }
    close(pb)
    cat("Calculations were performed successfully",'\n')
    
    parr=data.frame(parr)
    names(parr) = names(par)
    # Correcting names bug
    rn = rev(names(parr))
    rn[1:r] = paste("eff.",r:1,sep="")
    names(parr) = rev(rn)
    if(max(fam)==1){names(parr)[13] = c("eff")}
    rownames(parr) = paste('SNP',1:m,sep='')
    
    # Return
    return(list("PolyTest"=parr,"Method"="Empirical Bayes","MAP"=MAP,"SNPs"=SNPnames,"VAR"=VAR,"Fam"=fam))}
  
  #########
  ## RUN ##
  #########
  
  y = as.vector(y)
  covariate = matrix(covariate)
  fit = RANDOMsma(gen=gen,MAP=MAP,fam=fam,chr=chr,y=y,COV=covariate,SNPnames=SNPs)
  
  # fixing output bug
  moda=function(x){
    it=5;ny=length(x);k=ceiling(ny/2)-1; while(it>1){
      y=sort(x); inf=y[1:(ny-k)]; sup=y[(k+1):ny]
      diffs=sup-inf; i=min(which(diffs==min(diffs)))
      M=median(y[i:(i+k)]); it=it-1}; return(M)}
  neg = which(fit$PolyTest$lrt<0);fit$PolyTest$lrt[neg]=0
  mo = moda(fit$PolyTest$lrt)
  mos = which(fit$PolyTest$lrt==mo);fit$PolyTest$lrt[mo]=0
  fit$Fam = FAM0
  whEff = grep('eff',names(fit$PolyTest))
  if(any(whEff)) names(fit$PolyTest)[whEff] = paste('eff',unique(FAM0),sep='.')
  rownames(fit$PolyTest) = fit$SNPs
  class(fit) <- "NAM"
  return(fit)}

# Meta analysis gor gwas3
meta3 = function(ByEnv,ammi=1){ 
  
  ML = function(par){ # Woodbury solution
    vg = par[1]; ve = par[2];  vi = par[3]
    U = cbind(z*vg,w*ve,u*vi)
    v = zz*(vg^2)+ww*(ve^2)+zw*(vi^2)+kk
    RU = crossprod(ikk,U)
    Capacitance = crossprod(U,RU)
    diag(Capacitance) = diag(Capacitance)+1
    detCAP = determinant.matrix(Capacitance)$modulus[1]
    detK = determinant.matrix(kk)$modulus[1]
    d = detCAP+detK # logDet
    Capacitance = chol2inv(Capacitance)
    vv = ikk - tcrossprod(RU,Capacitance)%*%crossprod(U,ikk)
    xpx = sum(vv)
    xpy = sum(crossprod(vv,y))
    mu = xpy/xpx
    r = y-mu
    like = -0.5*( d + sum(crossprod(vv,r)*r) )
    return(-like)
  }
  
  MLtest = function(par){ # Woodbury solution
    vg = par[1]; ve = par[2];  vi = par[3]; mu = par[4];
    U = cbind(z*vg,w*ve,u*vi)
    v = zz*(vg^2)+ww*(ve^2)+zw*(vi^2)+kk
    RU = crossprod(ikk,U)
    Capacitance = crossprod(U,RU)
    diag(Capacitance) = diag(Capacitance)+1
    detCAP = determinant.matrix(Capacitance)$modulus[1]
    detK = determinant.matrix(kk)$modulus[1]
    d = detCAP+detK # logDet
    Capacitance = chol2inv(Capacitance)
    vv = ikk - tcrossprod(RU,Capacitance)%*%crossprod(U,ikk)
    r = y-mu
    like = -0.5*( d + sum(crossprod(vv,r)*r) )
    return(-like)
  }
  
  fixed<-function(par){
    vg = par[1]; ve = par[2];  vi = par[3]
    U = cbind(z*vg,w*ve,u*vi)
    v = zz*(vg^2)+ww*(ve^2)+zw*(vi^2)+kk
    RU = crossprod(ikk,U)
    Capacitance = crossprod(U,RU)
    diag(Capacitance) = diag(Capacitance)+1
    Capacitance = chol2inv(Capacitance)
    vv = ikk - tcrossprod(RU,Capacitance)%*%crossprod(U,ikk)
    xpx = sum(vv)
    xpy = sum(crossprod(vv,y))
    mu = xpy/xpx
    vb<-1/xpx
    return(c(mu,vb))
  }
  
  E = length(ByEnv)
  
  # unique set of SNPs
  snps = ByEnv[[1]]$SNPs
  for(j in 2:E) snps = intersect(snps,ByEnv[[j]]$SNPs)
  for(j in 1:E) rownames(ByEnv[[j]]$PolyTest) = ByEnv[[j]]$SNPs
  for(j in 1:E) ByEnv[[j]]$PolyTest = ByEnv[[j]]$PolyTest[snps,]
  for(j in 1:E) ByEnv[[j]]$VAR = ByEnv[[j]]$VAR[,,which(ByEnv[[j]]$SNPs%in%snps)]
  m = length(snps)
  
  map = ByEnv[[1]]$MAP[which(ByEnv[[1]]$SNPs%in%snps),c(1,3)]
  for(j in 2:E) map = map+ByEnv[[j]]$MAP[which(ByEnv[[j]]$SNPs%in%snps),c(1,3)]
  map = round(map/E,4)
  map = data.frame(snps,map)
  names(map) = c("snpID","chr","pos")
  map = map[,c(2,3,1)]
  
  # Dealing with environments with not all families
  
  Vi = rep(0,E); for(j in 1:E) Vi[j] = dim(ByEnv[[j]]$VAR[,,1])[1]
  cVi = c(0,cumsum(Vi))
  dimK = sum(Vi)
  kk = ikk = matrix(0,dimK,dimK)
  p = list()
  
  fams = c()
  for(j in 1:E) fams = c(fams,unique(ByEnv[[j]]$Fam))
  levF=c('0',sort(unique(fams)))
  
  for(j in 1:E){
    f = factor( ByEnv[[j]]$Fam ,levF)
    p[[j]]= sort(unique(as.numeric(f)))
  } 
  
  # Preparing inputs, output matrix and design matrices
  p = as.factor(unlist(p))
  e = c()
  for(j in 1:E) e = c(e,rep(j,Vi[j]))
  e = factor(e)
  y = rep(NA,length(p))
  z = model.matrix(~p-1)
  w = model.matrix(~e-1)
  zz= tcrossprod(z)
  ww= tcrossprod(w)
  x = matrix(1,dimK,1)
  u = model.matrix(~p:e-1)
  zw= tcrossprod(u)
  
  # Starting GxE gwas
  ee = as.numeric(e)
  pp = as.numeric(p)
  maxP = length(levels(p))
  
  ammi2 = function(y,pc=1){
    yy = lm(y~p+e)$residuals
    ge = matrix(0,maxP,E)
    for(k in 1:length(y)) ge[pp[k],ee[k]] = yy[k]
    S = svd(ge,pc,pc)
    SVD = matrix(0,length(y),pc)
    for(o in 1:pc){
      SS = tcrossprod(S$u[,o],S$v[,o])
      for(k in 1:length(y)) SVD[k,o] = SS[pp[k],ee[k]]
    }
    return(SVD)}
  
  ################
  ###   LOOP   ###  
  ################
  
  output = matrix(NA,m,8)
  colnames(output) = c('Mu','Var_Mu','Var_G','Var_E','Var_GE','L1','L0','lrt')
  
  cat(' \n')
  
  pb=txtProgressBar(style=3)
  for(i in 1:m){
    # generate known residual covariance matrix
    for(j in 1:E){
      k1 = ByEnv[[j]]$VAR[,,i]
      k2 = chol2inv(k1)
      kk[(cVi[j]+1):cVi[j+1],(cVi[j]+1):cVi[j+1]]=k1
      ikk[(cVi[j]+1):cVi[j+1],(cVi[j]+1):cVi[j+1]]=k2
    } 
    # collect the regression coefficient of marker effect
    for(j in 1:E) y[(cVi[j]+1):cVi[j+1]] = ByEnv[[j]]$PolyTest[i,-c(1:12)]
    y = as.matrix(as.numeric(y),ncol=1)
    
    u = ammi2(y,ammi)
    zw = tcrossprod(u)
    l0 = -MLtest(rep(0,4))
    fit = optim(par=rep(1e-10,3),fn=ML,method="L-BFGS-B",lower=0,upper=1e-6)
    beta = fixed(fit$par)
    mu = beta[1]
    vb = beta[2]
    l1 = -MLtest(c(fit$par,mu))
    lrt = 2*(l1-l0)
    output[i,]=c(mu,vb,fit$par,l1,l0,lrt)
    setTxtProgressBar(pb,i/m)
  }
  close(pb)
  
  cat('Meta-analysis was performed successfully \n')
  
  output = data.frame(output)
  final = list('PolyTest'=output,'MAP'=map,'Method'="MetaGWA",'SNPs'=snps)
  
  class(final) = "NAM"
  return(final)
}

# NAM-type matrix conversion and eigendecomposition
eigX = function(gen,fam){
  if(any(is.na(gen))){
    mis = apply(gen,2,function(x) mean(is.na(x)))
    gen = gen[,mis<.8]
    gen = markov(gen,ncol(gen))
  }
  gm=function(gen,fam){
    g1 = tcrossprod(gen)
    g1 = g1/mean(diag(g1))
    g2 = ((gen-1)*-1)+1
    g2 = tcrossprod(g2)
    g2 = g2/mean(diag(g2))
    for (i in unique(fam)){nam = which(fam == i)
    g1[nam, nam] = g2[nam, nam] + g1[nam, nam]}
    lambda = mean(diag(g1))
    G = g1/mean(diag(g1))
    return(G)}
  return(eigen(gm(gen,fam),symmetric=TRUE))
}

# NAM plot function
plot.NAM = function(x,...,alpha=0.05,colA=2,colB=4,find=NULL,FDR=NULL,gtz=FALSE,phys=NULL){
  anyNA = function(x) any(is.na(x))
  if(!is.null(FDR)){if(FDR>=1|FDR<0)stop("FRD must be between 0 and 1")}
  gwas=x
  chr=as.numeric(summary(factor(as.numeric(gwas$MAP[,1]))))
  
  if(gwas$Method=="P3D"){
    
    FGWASplot=function(Fgwas,chr,AA,BB,...){
      col=c();for(i in 1:length(chr)){if((i%%2)==0){b=AA}else{b=BB};a=rep(b,chr[i]);col=c(col,a)}
      W = Fgwas$PolyTest$wald
      if(is.null(phys)){ Xaxis = 1:length(W) }else{ Xaxis = cumsum(phys) }
      plot(Xaxis,W,col=col,xlab="Chromosome",ylab="Wald Statistics",xaxt = "n",...)
      return(W)}
    
    # Plot
    pv=FGWASplot(gwas,chr=chr,AA=colA,BB=colB,...)
    
    # Computing where chromosomes start and end
    NumChr = length(chr)
    Ch0 = cumsum(c(1,chr[-NumChr]))-.5
    Ch1 = cumsum(chr)+.5
    
    # Adding Chromosome in X axis
    medians=rep(NA,length(chr))
    if(is.null(phys)){
      for(i in 1:length(chr)) medians[i] = median(which(gwas$MAP[,1]==i))
      axis(1, at=round(medians), labels=1:length(medians))
    }else{
      Xaxis = cumsum(phys)
      #for(i in 1:length(chr)) medians[i] = median(Xaxis[gwas$MAP[,1]==i])
      for(i in 1:length(chr)) medians[i] = mean(range(Xaxis[gwas$MAP[,1]==i]))
      axis(1, at=round(medians), labels=1:length(medians))
    }
    
    # QTL
    if(!is.null(find)){Loc=identify(n=find,x=1:length(pv),y=pv,labels=gwas$SNPs);for(i in Loc) cat(gwas$SNPs[i],'\n')}
    
  }else{
    
    RGWASplot=function(Rgwas,chr,AA,BB,...){
      # Colors
      col=c();for(i in 1:length(chr)){if((i%%2)==0){b=AA}else{b=BB};a=rep(b,chr[i]); col=c(col,a)}
      # Statistics
      
      if(gwas$Method=="MetaGWA"){
        S = Rgwas$PolyTest$lrt
        YLAB = "LRT"
      }else{
        S = Rgwas$PolyTest$pval
        YLAB = "-log(p-value)"
      }
      
      if(is.null(phys)){ Xaxis = 1:length(S) }else{ Xaxis = cumsum(phys) }
      plot(Xaxis,S,col=col,xlab="Chromosome",ylab=YLAB,xaxt = "n",...)
      return(S)
    }
    
    # Plot
    if (is.null(alpha)|gwas$Method=="MetaGWA"){
      
      pv=RGWASplot(gwas,chr=chr,AA=colA,BB=colB,...)
      
    } else {
      
      pv=RGWASplot(gwas,chr=chr,AA=colA,BB=colB,...)
      
      # Computing where chromosomes start and end
      NumChr = length(chr)
      Ch0 = cumsum(c(1,chr[-NumChr]))-.5
      Ch1 = cumsum(chr)+.5
      
      
      if(is.null(FDR)){
        A = 1-alpha
        LRmax = qchisq(A,0.5)
        lim = -log(pchisq(LRmax, 0.5,lower.tail = FALSE),base = 10)
        abline(h=lim,col=1,lty=2)
        
      }else{
        
        # Multiple test correction
        if(gtz==T){
          
          MT = tapply(
            X = gwas$PolyTest$pval,
            INDEX = gwas$MAP[,1],
            FUN = function(x) sum(x>0)
          )
          MT[MT==0]=1
          
        }else{MT=chr}
        
        
        # Thresholds
        for(i in 1:(NumChr)){
          A = 1-alpha/(MT[i]*(1-FDR))
          LRmax = qchisq(A,0.5)
          lim = -log(pchisq(LRmax, 0.5,lower.tail=FALSE),base = 10)
          lines(x = c(Ch0[i],Ch1[i]),y = c(lim,lim))
        }
        
      }
    }
    
    # Adding Chromosome in X axis
    medians=rep(NA,length(chr))
    if(is.null(phys)){
      for(i in 1:length(chr)) medians[i] = median(which(gwas$MAP[,1]==i))
      axis(1, at=round(medians), labels=1:length(medians))
    }else{
      Xaxis = cumsum(phys)
      #for(i in 1:length(chr)) medians[i] = median(Xaxis[gwas$MAP[,1]==i])
      for(i in 1:length(chr)) medians[i] = mean(range(Xaxis[gwas$MAP[,1]==i]))
      axis(1, at=round(medians), labels=1:length(medians))
    }
    
    # QTL
    if(!is.null(find)){Loc=identify(n=find,x=1:length(pv),y=pv,labels=gwas$SNPs);for(i in Loc) cat(gwas$SNPs[i],'\n')}
    
    
  }
  
  
}

# Change reference for GWAS
reference <- function(gen,ref=NULL){
    anyNA = function(x) any(is.na(x))
    if(is.null(ref)){
      maf=function(z){
        Z=z
        z[z==2]="A"; nA=length(which(z=="A"))
        z[z==0]="B"; nB=length(which(z=="B"))
        if(nA==nB){x=Z}
        if(nA>nB){z[z=="A"]=2;z[z=="B"]=0;x=z}
        if(nA<nB){z[z=="B"]=2;z[z=="A"]=0;x=z}
        return(x)}
      W=apply(gen,2,maf);W=(as.numeric(W));W=matrix(W,ncol=ncol(gen))
      colnames(W)=colnames(gen)
      
    }else{
      if(ncol(gen)!=length(ref)) stop("Reference parent and matrix of genotypes display non compatible dimensions")
      if(any(is.na(ref))||length(which(ref==5))>0) {"Reference parent must have no missing values"}
      gen[is.na(gen)]=5    
      CA = which(ref==0) # Changing Alleles
      W = gen-1
      W[,CA] = W[,CA]*-1
      W[W==(-4)]=4
      W = W+1
      W[W==(5)]=NA
    } 
    return(W)}

