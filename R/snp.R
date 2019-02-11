
# Snp Heritability
snpH2=function(gen,K=NULL){
  if(!is.numeric(gen)) stop("Object gen must to be a numeric matrix")
  if(!is.null(K)){ if(!is.numeric(K)) stop("Object K must to be a numeric matrix")  }
  anyNA = function(x) any(is.na(x))
  # Checking for missing and imputing
  if(any(is.na(gen))) gen = IMP(gen)
# Setup
snps = ncol(gen)
obs = nrow(gen)
x = matrix(1:obs,ncol=1)
# If pedigree is not provided
if(is.null(K)){
  K = tcrossprod(gen)
  K = K/mean(diag(K))}
# Defining log-REML
loglike<-function(theta){
  lambda<-exp(theta)
  logdt<-sum(log(lambda*delta+1))
  h<-1/(lambda*delta+1)
  yy<-sum(yu*h*yu)
  yx<-matrix(0,q,1)
  xx<-matrix(0,q,q)
  for(i in 1:q){
    yx[i]<-sum(yu*h*xu[,i])
    for(j in 1:q){xx[i,j]=sum(xu[,i]*h*xu[,j])}}
  loglike = -0.5*logdt-0.5*(n-q)*log(yy-t(yx)%*%solve(xx)%*%yx)-0.5*log(det(xx))
  return(-loglike)}
fixed<-function(lambda){
  h<-1/(lambda*delta+1)
  yy<-sum(yu*h*yu)
  yx=timesVec(yu,h,xu,q)
  xx=timesMatrix(xu,h,xu,q,q)
  beta<-qr.solve(xx,yx)
  sigma2<-(yy-t(yx)%*%solve(xx)%*%yx)/(n-q)
  sigma2 = as.numeric(sigma2)
  var<-diag((solve(xx))*sigma2)
  stderr<-sqrt(var)
  return(c(beta,stderr,sigma2))}
# Eigendecomposition of K
qq<-eigen(as.matrix(K),symmetric=T)
delta<-qq[[1]]
uu<-qq[[2]]
q<-1
n<-ncol(K)
H2 = c()
# LOOP
pb=txtProgressBar(style=3)
for(i in 1:snps){
  y = gen[,i]
  yu<-t(uu)%*%y
  xu<-t(uu)%*%x
  vp<-var(y)
  # Finding lambda through optimization
  parm<-optim(par=10,fn=loglike,method="L-BFGS-B",lower=1e-12,upper=1e12)
  lambda<-exp(parm$par)
  # Results
  parmfix<-fixed(lambda)
  Ve<-parmfix[2*q+1]
  Vg<-lambda*Ve
  h2=Vg/(Vg+Ve)
  # Saving loop
  H2 = c(H2,h2)
  setTxtProgressBar(pb,i/snps)}
close(pb)
names(H2)=colnames(gen)
class(H2)="H2"
return(H2)
}

# Plot H2
plot.H2 = function(x,...,chr=NULL){
  anyNA = function(x) any(is.na(x))
  her = as.numeric(x)
  plot(her,xlab="Genome",ylab="Heritability",main="Gene content")
  rect(par("usr")[1],par("usr")[3],
       par("usr")[2],par("usr")[4],
       col=rgb(0.25,0.25,0.25,0.5))
  lines(her,lwd=4,pch=20,type="p",col="blue")
  lines(her,lwd=1,pch=20,type="p",col="lightblue")
  if(is.null(chr)!=T) abline(v=cumsum(chr[-length(chr)])+0.5,lty=3)
}

# Fixation index functions
Fst = function(gen,fam){ 
  anyNA = function(x) any(is.na(x))
  # Weir and Cockerham (1984): Unbiased Fst from weighted AOV 
  if(anyNA(fam)) stop("No missing fam is allowed")
  if(anyNA(gen)) gen[is.na(gen)]=1 # replacing missing
  # number of markers and populations
  m = ncol(gen)
  r = length(unique(fam))
  n = table(fam)
  nbar = mean(n)
  # function to calculate binomial variance
  Var = function(snp){p=mean(snp,na.rm=T)/2;return(p*(1-p))}
  # function to estimate p and heterozigousity by family
  p_i = function(snp,fam) tapply(snp,fam,mean)/2 
  h_i= function(snp,fam) tapply(snp,fam,function(x)mean(x==1,na.rm=T))
  # Whole-Population parameters
  C2 = var(n)/(mean(n)^2) # squared CV of sample sizes
  nc = nbar*(1-C2/r)
  # Parameters by population
  M = apply(gen,2,p_i,fam=fam) # observed P by family by allele
  N = matrix(n,r,m)
  Pbar = colSums(M*N)/(r*nbar) # expected P
  eP = t(matrix(Pbar,m,r))
  s2 = colSums(N*(M-eP)^2)/((r-1)*nbar) # var over pop
  H = apply(gen,2,h_i,fam=fam) # observed heterozigousity
  Hbar = colSums(H*N)/(r*nbar)
  V = Pbar*(1-Pbar)
  # Three variance components: A, B and C
  Va = nbar/nc*(s2-(1/(nbar-1))*(V-(r-1)/r*s2-Hbar/4))
  Vb = (nbar/(nbar-1))*(V-(r-1)/r*s2-(2*nbar-1)/(4*nbar)*Hbar)
  Vc = Hbar/2
  # F - statistics!
  fst = Va/(Va+Vb+Vc); fst[is.nan(fst)|fst<=0]=0
  fit = 1-Vc/(Va+Vb+Vc); fit[is.nan(fit)|fit<=0]=0
  fis = 1-Vc/(Vb+Vc); fis[is.nan(fis)|fis<=0]=0
  Fstat = data.frame(fst,fit,fis)
  class(Fstat) <- "fst"
  return(Fstat)
}

# Plot FST
plot.fst = function(x,..., p=NULL,chr=NULL){
  anyNA = function(x) any(is.na(x))
  class(x)=NULL; FST=x$fst
  C = 1 - 0.1 * (4 + 6 * abs(FST)/max(abs(FST)))
  
  if(is.null(chr)!=T){
    plot(FST,...,xaxt = "n",main="Fixation Index",pch=20,xlab="Genome",col=(rgb(C,C,C)))
  }else{
    plot(FST,...,main="Fixation Index",pch=20,xlab="Genome",col=(rgb(C,C,C)))
  }
  K = smooth(smooth(FST))
  lines(K,type="l",lwd=1,col=rgb(0,0,1,0.6))
  if(is.null(chr)!=T){
    # Adding Chromosome in X axis
    medians=rep(NA,length(chr))
    CHR=c(); for(i in 1:length(chr)) CHR=c(CHR,rep(i,chr[i]))
    for(i in 1:length(chr)) medians[i] = median(which(CHR==i))
    axis(1, at=round(medians), labels=1:length(medians))
    abline(v=cumsum(chr[-length(chr)])+0.5,lty=3)
  } 
  
  
  # P-value  
  if(is.null(p)!=TRUE){
    par(ask=TRUE); CDF = ecdf(FST);
    Pval = -log(1-CDF(FST)) ; THR = -log(p/length(FST))
    plot(Pval,ylab="Empirical -log probatility",
         xlab="Genome",main="P-values",pch=20,col=1+1*(Pval>THR))
    if(is.null(chr)!=T) abline(v=cumsum(chr[-length(chr)])+0.5,lty=3)
    abline(h=THR,col=rgb(1,0,0,0.4),lty=2);par(ask=FALSE)
    KK = ksmooth(1:length(FST),Pval,"normal",bandwidth = 10)
    lines(KK,type="h",lwd=2,col=rgb(1,0.5,0,0.6))}
}

# function for marker quality control
snpQC=function(gen,psy=1,MAF=0.05,misThr=0.8,remove=TRUE,impute=FALSE){
  anyNA = function(x) any(is.na(x))
  
  # CHECKING REDUNDANT MARKERS
  gen2=gen; redundancy=c(0); for(i in 1:(ncol(gen)-1)){
    a=mean((gen[,i]==gen[,(i+1)]),na.rm=TRUE);redundancy=c(redundancy,a)} 
  a=which(redundancy>psy);b=length(which(redundancy>psy))
  if(b>0){cat("Genotypic data contains",b,"redundant SNPs",'\n')
    if(remove==TRUE){gen2=gen[,-a]}
  }else{cat("No redundant SNPs found",'\n')}
  
  # REMOVE SNP A LOT OF MISSING # NEW!
  msnp = apply(gen2,2,function(x) mean(is.na(x)) )
  noVal = which(msnp>misThr)
  
  # CHECKING MINOR ALLELE FREQUENCY
  if(MAF>0){
    Count_maf = function(x){
      AA=sum(x==2)
      Aa=sum(x==1)
      aa=sum(x==0)
      Total=AA+Aa+aa
      PA=(AA+0.5*Aa)/Total
      Pa=(aa+0.5*Aa)/Total
      lowerAF=min(PA,Pa)
      return(lowerAF)}
    LAF = apply(gen2,2,Count_maf)
    maf=which(LAF<MAF)
    
    # REMOVE SNP A LOT OF MISSING # NEW!
    msnp = apply(gen2,2,function(x) mean(is.na(x)) )
    noVal = which(msnp>misThr)
    
    # hist(LAF,col=3,nclass=50,main="Histogram of MAF",xlab="Minor Allele Frequency")
    if(length(maf)>0){
      cat("There are",length(maf),"markers with MAF below the threshold",'\n')
      if(remove==TRUE){
        gen3=gen2[,-maf]
        if(any(noVal)){gen3=gen3[,-noVal]}
      }else{
        gen3=gen2
      }
      
    }else{cat("No marker below MAF threshold",'\n');
      gen3=gen2
      if(any(noVal)){gen3=gen3[,-noVal]}
      }
  }else{gen3=gen2}
  if(impute){
    rf <- function(xmis){ # Author: D.Stekhoven, stekhoven@stat.math.ethz.ch
      maxiter = 10; ntree = 100
      mtry = floor(sqrt(ncol(xmis)))
      cutoff=NULL;classwt=NULL;strata=NULL;replace=T
      sampsize=NULL;nodesize=NULL;maxnodes=NULL
      n <- nrow(xmis);  p <- ncol(xmis)
      ## perform initial S.W.A.G. on xmis (mean imputation)
      ximp <- xmis; xAttrib <- lapply(xmis, attributes)
      varType <- character(p);  for (t.co in 1:p){
        if (is.null(xAttrib[[t.co]])){ varType[t.co] <- 'numeric'
        ximp[is.na(xmis[,t.co]),t.co] <- mean(xmis[,t.co], na.rm = TRUE)
        }else{ varType[t.co] <- 'factor'
        max.level <- max(table(ximp[,t.co]))
        class.assign <- sample(names(which(max.level == summary(ximp[,t.co]))), 1)
        if (class.assign != "NA's"){ximp[is.na(xmis[,t.co]),t.co] <- class.assign } else {
          while(class.assign=="NA's"){class.assign=sample(names(which(max.level==summary(ximp[,t.co]))), 1)}
          ximp[is.na(xmis[,t.co]),t.co] <- class.assign }}}
      NAloc <- is.na(xmis);  noNAvar <- apply(NAloc, 2, sum) 
      sort.j <- order(noNAvar); sort.noNAvar <- noNAvar[sort.j]
      nzsort.j <- sort.j[sort.noNAvar > 0]; Ximp <- vector('list', maxiter)
      iter <- 0; k <- length(unique(varType));  convNew <- rep(0, k)
      convOld <- rep(Inf, k); OOBerror <- numeric(p); names(OOBerror) <- varType
      if (k==1){if(unique(varType)=='numeric'){names(convNew)=c('numeric')}else{names(convNew)=c('factor')}
        convergence <- c(); OOBerr <- numeric(1)} else {
          names(convNew) <- c('numeric', 'factor');convergence <- matrix(NA, ncol = 2);OOBerr <- numeric(2)}
      stopCriterion <- function(varType, convNew, convOld, iter, maxiter){
        k <- length(unique(varType)); if (k == 1){(convNew < convOld) & (iter < maxiter)} else {
          ((convNew[1] < convOld[1]) | (convNew[2] < convOld[2])) & (iter < maxiter)}}
      while (stopCriterion(varType, convNew, convOld, iter, maxiter)){
        if (iter != 0){convOld <- convNew;OOBerrOld <- OOBerr}
        cat("RF iteration", iter+1, "\n");t.start <- proc.time();ximp.old <- ximp
        for(s in 1:p){varInd=sort.j[s];if(noNAvar[[varInd]]!=0){obsi=!NAloc[,varInd];misi=NAloc[, varInd];
        obsY=ximp[obsi,varInd];obsX=ximp[obsi, seq(1,p)[-varInd]];misX=ximp[misi, seq(1,p)[-varInd]];
        typeY=varType[varInd];RF=randomForest(x=obsX,y=obsY,ntree=ntree,mtry=mtry,replace=T,
                                              sampsize=if(!is.null(sampsize))sampsize[[varInd]] else if (replace) nrow(obsX) else
                                                ceiling(0.632*nrow(obsX)), nodesize = if (!is.null(nodesize)) nodesize[1] else 1,
                                              maxnodes = if (!is.null(maxnodes)) maxnodes else NULL); OOBerror[varInd] <- RF$mse[ntree]
                                              misY <- predict(RF, misX); ximp[misi, varInd] <- misY }}; iter <- iter+1;Ximp[[iter]] <- ximp
        t.co2 <- 1;for (t.type in names(convNew)){ t.ind <- which(varType == t.type)
        convNew[t.co2]=sum((ximp[,t.ind]-ximp.old[,t.ind])^2)/sum(ximp[,t.ind]^2);t.co2=t.co2 + 1}}
      if (iter == maxiter){out <- Ximp[[iter]]}else{out <- Ximp[[iter-1]]}
      return(out)}
    gen=gen3;gen[gen==5]=NA   
    k=100*length(which(is.na(gen)))/length(gen)
    k=round(k,2);cat(k,"% of missing data",'\n')
    cat("Imputations being performed by Random Forest",'\n')
    if(any(is.na(gen))){gen=suppressWarnings(rf(gen));gen=round(gen)}
    gen3=gen}
  return(gen3)}

# function to remove repeated genotypes
cleanREP = function(y,gen,fam=NULL,thr=0.95){
  if(is.vector(y)) y=matrix(y,ncol=1)
  if(is.null(fam)) fam = rep(1,nrow(y))
    
  ibs = function(gen){
  f1 = function(x,gen) apply(gen,1,function(y,x) mean(y==x,na.rm = T),x=x)
  f2 = apply(gen,1,f1,gen=gen)
  return(f2)}  
  
  GG = function(gen, r = 1) {
    a1 = (gen - 1)
    a1[a1 == -1] = 0
    A1 = (tcrossprod(a1))
    a2 = -(gen - 1)
    a2[a2 == -1] = 0
    A2 = (tcrossprod(a2))
    d = round(exp(-abs(gen - 1)))
    D = tcrossprod(d)
    G = A1 + A2 + D
    G = (G/ncol(gen))^r
    return(G)
  }
  cat("solving identity matrix\n")
  G = GG(gen)
  
  rownames(G) = 1:nrow(G)
  lt = G*lower.tri(G) # lower triang
  r = 1* lt>thr # logical matrix: repeatitions
  # starting point of new data
  rownames(gen) = 1:nrow(gen)
  Ny=y;  Nfam=fam;  Ngen=gen
  # summary
  cs = colSums(r) # how many times id will be repeated
  while(any(cs>0)){
    i = which(cs>0)[1]
    cat("indiviual",rownames(gen)[i],"had",cs[i],'duplicate(s)\n')
    w = which(r[,i])
    if(ncol(y)>1){y[i,] = colMeans(y[c(i,w),],na.rm=T)
    }else{y[i] = mean(y[c(i,w)],na.rm=T)}
    if(ncol(y)>1){Ny=Ny[-w,]}else{Ny=Ny[-w]}
    Nfam=Nfam[-w]
    Ngen=Ngen[-w,]
    r = r[-w,]
    cs = colSums(r)
  }
  return(list(y=Ny,gen=Ngen,fam=Nfam))
}

# Some sort of Hidden Markov model for imputation
markov=function(gen,chr=NULL){
  if(is.null(chr)) chr = ncol(gen)
  # vector chr
  CHR=NULL;for(i in 1:length(chr)){CHR=c(CHR,rep(i,chr[i]))}
  # Expectation and Transition Probability
  tr = function(v1,v2){ 
    tp=rep(NA,9) # Transition Probability
    tp[1]=mean(v1==0&v2==0,na.rm=T);tp[2]=mean(v1==0&v2==1,na.rm=T);tp[3]=mean(v1==0&v2==2,na.rm=T)
    tp[4]=mean(v1==1&v2==0,na.rm=T);tp[5]=mean(v1==1&v2==1,na.rm=T);tp[6]=mean(v1==1&v2==2,na.rm=T)
    tp[7]=mean(v1==2&v2==0,na.rm=T);tp[8]=mean(v1==2&v2==1,na.rm=T);tp[9]=mean(v1==2&v2==2,na.rm=T)
    tp[tp==0]=1e-5;tp[1:3]=tp[1:3]/sum(tp[1:3]);tp[4:6]=tp[4:6]/sum(tp[4:6]);tp[7:9]=tp[7:9]/sum(tp[7:9])
    return(tp)}
  # Transition matrix
  TM = function(gen){M = ncol(gen); N = nrow(gen)
  step1 = rbind(gen[,-M],gen[,-1])
  step2 = function(snps) tr(snps[1:N],snps[-c(1:N)])
  step3 = apply(step1,2,step2); rm(step1)
  rownames(step3) = paste(gl(3,3,9,0:2),0:2,sep='to')
  return(step3)}
  # Calculate log-prob of transitions
  mis = gen;  tm=log(TM(gen))
  # Imputation with Expectation
  IE = function(v1,v2,tp){
    exp=rep(NA,3);exp[1]=which.max(tp[1:3])-1;exp[2]=which.max(tp[4:6])-1;exp[3]=which.max(tp[7:9])-1
    w=which(is.na(v2));r=v1[w]+1;v=exp[r];v2[w]=v;return(v2)}
  # Imputing first row (starting point)
  gen[,1] = IE(IE(IE(gen[,4],gen[,3],tm[,3]),gen[,2],tm[,2]),gen[,1],tm[,1])
  if(anyNA(gen[,1]))gen[,1][is.na(gen[,1])]=as.numeric(names(which.max(table(gen[,1],exclude=NA))))
  # Imputing the rest
  for(i in 2:ncol(gen)) gen[,i]=IE(gen[,i-1],gen[,i],tm[,i-1])
  # THE END
  return(gen)
}

# LD matrix function
LD = function(gen){
  # Phasing via EM
  EM=function(A,B,n=10){
    Za = paste(A,B,sep="")
    Z = matrix(0,3,3)
    colnames(Z)=c('A0','A1','A2')
    rownames(Z)=c('B0','B1','B2')
    # Main diagonal
    Z[1,1] = sum(Za=='00')
    Z[2,2] = sum(Za=='11')
    Z[3,3] = sum(Za=='22')
    # Upper diagonal
    Z[1,2] = sum(Za=='01')
    Z[1,3] = sum(Za=='02')
    Z[2,3] = sum(Za=='12')
    # Lower diagonal
    Z[1,2] = sum(Za=='10')
    Z[1,3] = sum(Za=='20')
    Z[2,3] = sum(Za=='21')
    # Initial guess
    COUP = 0.5 ; REPU = 0.5
    # Function to estimate haplotypes
    HapProb=function(Z,Co,Re){
      AB = 2*Z[1,1] + Z[2,1] + Z[1,2] + Co*Z[2,2]
      Ab = 2*Z[3,1] + Z[2,1] + Z[3,2] + Re*Z[2,2]
      aB = 2*Z[1,3] + Z[1,2] + Z[2,3] + Re*Z[2,2]
      ab = 2*Z[3,3] + Z[3,2] + Z[2,3] + Co*Z[2,2]
      props = data.frame(AB,Ab,aB,ab)
      haps=(data.matrix(props)/2)/sum(Z)
      rownames(haps)="Hap"
      return(haps)}
    HHH=matrix(NA,n,6)
    # Loop
    for(i in 1:n){
      H=HapProb(Z,COUP,REPU) 
      H2=matrix(H,2,2)
      dia=H2[1,1]*H2[2,2]
      off=H2[2,1]*H2[1,2]
      tt=dia+off
      COUP = dia/tt
      REPU = off/tt
      oo = c(COUP,REPU)
      oo[is.nan(oo)]=0
      H[is.nan(H)]=0
      HHH[i,]=c(H,oo)
    }
    rownames(HHH)=1:n
    EM=round(tail(HHH,1)[1:4],4)
    names(EM)=colnames(HHH)[1:4]
    return(EM)
  }
  # LD values
  EM_LD = function(A,B){
    phase=EM(A,B)
    if(det(matrix(phase,2,2))<0) phase=phase[c(2,1,4,3)]
    X = matrix(phase,2,2)
    D = X[1,1]-sum(X[,1])*sum(X[1,])
    r=D/sqrt(prod(c(colSums(X),rowSums(X))))
    return(r^2)
  }
  # LD matrix
  LDmat = function(gen){
    n = ncol(gen)
    snps = colnames(gen)
    LD1 = matrix(1,n,n)
    dimnames(LD1) = list(snps,snps)
    for(i in 1:n){for(j in 1:n){if(j>i){
      ld = as.numeric(EM_LD(gen[,i],gen[,j]))
      LD1[i,j]=LD1[j,i]=ld
    }}}
    return(LD1)}
  return(LDmat(gen))
}

# Compute genetic distance
Gdist = function(gen,method=1){
  # from adegenet{dist.genpop}
  nloc = ncol(gen)
  recode = function(gen) cbind(gen,(gen-2)*-1)
  X = recode(gen)  
  # NEI
  if (method == 1) {
    cat("Nei distance\n")
    d = tcrossprod(X)
    vec = sqrt(diag(d))
    d = d/vec[col(d)]
    d = d/vec[row(d)]
    d = -log(d)
    d = as.dist(d)}
  # Edwards
  else if (method == 2) {
    cat("Edwards distance\n")
    X2 = sqrt(X)
    d = tcrossprod(X2)
    d = 1 - d/(nloc*2)
    diag(d) = 0
    d = sqrt(d)
    d = as.dist(d)}
  # Reynolds
  else if (method == 3){
    cat("Reynolds distance\n")
    X = X/2
    denomi = tcrossprod(X)
    vec = rowSums((X)^2)
    d = -2 * denomi + vec[col(denomi)] + vec[row(denomi)]
    diag(d) = 0
    denomi = 2 * (nloc - denomi)
    diag(denomi) = 1
    d = d/denomi
    d = sqrt(d)
    d = as.dist(d)}
  # Rogers
  else if (method == 4) {
    cat("Rogers' distance\n")
    p = ncol(gen)
    d = dist(gen,method = 'euclidean')
    d = d/p}
  # Prevosti
  else if (method == 5){
    cat("Prevosti's distance\n")
    p = ncol(gen)
    d = dist(gen,method = 'manhattan')
    d = d/(2*p)}
  # MRD
  else if (method == 6){
    cat("Modified Rogers' distance\n")
    p = ncol(gen)
    d = dist(gen,method = 'euclidean')
    d = d/sqrt(4*p)}
  return(d)
}

# G2A Kernels
G2A_Kernels = function(gen){
  # Centralizing G2A
  Z = apply(gen,2,function(x) x-mean(x,na.rm = T));Z[is.na(Z)]=0
  gen[gen==2]=0
  W = apply(gen,2,function(x) x-mean(x,na.rm = T));W[is.na(W)]=0
  # Additive
  ZZ = tcrossprod(Z)
  A = ZZ/mean(diag(ZZ))
  # Dominant
  WW = tcrossprod(W)
  D = WW/mean(diag(WW))
  # Epistatic additive
  ZZ = tcrossprod(Z)^2 - tcrossprod(Z^2)
  AA = ZZ/mean(diag(ZZ))
  # Epistatic dominant
  WW = tcrossprod(W)^2 - tcrossprod(W^2)
  DD = WW/mean(diag(WW))
  # Epistatic add x dom
  AD = A*D
  AD = AD/mean(diag(AD))
  # List
  KS = list('A'=A,'D'=D,'AA'=AA,'DD'=DD,'AD'=AD)
  # Return
  return(KS)
}

# IMPORT DATA
Import_data = function(file,type=NULL){
  
  if(is.null(type)) stop("'type' must be specified with one the following: 'GBS','HapMap','VCF'")
  
  # GBS
  GBS = function(path){
    G = read.delim(path)
    rownames(G)=G[,1];G=G[,-1]
    G[G!="A"&G!="C"&G!="G"&G!="T"]=NA
    Recode=function(aa){
      t=table(aa)
      t[t==0]=NA;A=which.max(t)
      a=which.min(t)
      A1=names(t)[A]
      A2=names(t)[a]
      W1=which(aa==A1)
      W2=which(aa==A2)
      aa=as.numeric(aa)
      aa[W1]=0;aa[W2]=2
      return(aa)
    }
    gen=apply(G,2,Recode)
    dimnames(gen)=dimnames(G)
    gen=data.matrix(gen)
    # chr
    chr = NULL; for(i in 1:max(G[,3])) chr=c(chr,sum(G[,3]==i))
    # return
    final = list('gen'=gen,'chr'=chr)
    return(final)
  }
  
  if(type=='GBS'){
    return(GBS(file))
  }
  
  # HAPMAP
  HAP = function(path){
    G = read.delim(path, header=T)
    AA=as.character(G[,2])
    AA=gsub('/','',AA)
    Gen=t(G[,-c(1:11)])
    n=nrow(Gen)
    m=ncol(Gen)
    gen=matrix(NA,n,m)
    for(i in 1:m){
      A1=strsplit(AA[i],'')[[1]][1]
      A2=strsplit(AA[i],'')[[1]][2]
      BB=paste(A1,A1,sep='')
      Bb=paste(A1,A2,sep='')
      bb=paste(A2,A2,sep='')
      M=as.character(Gen[,i])
      gen[M==BB,i]=2
      gen[M==Bb,i]=1
      gen[M==bb,i]=0
    }
    colnames(gen)=paste(G[,3],G[,4],G[,2],sep='.')
    # chr
    chr = NULL; for(i in 1:max(G[,3])) chr=c(chr,sum(G[,3]==i))
    rm(AA,A1,A2,BB,Bb,bb,m,n,Gen)
    # return
    final = list('gen'=gen,'chr'=chr)
    return(final)
  }
  
  if(type=='HapMap'){
    return(HAP(file))
  }
  
  # VCF
  VCF = function(path){
    G = read.table(path,header = TRUE)
    if(any(duplicated(G$POS))) G = G[-which(duplicated(G$POS)),]
    rownames(G) = paste(G$CHROM,G$POS,sep='_')
    G = G[,-c(1:9)]
    tags = dimnames(G)
    recode = function(g){
      g = as.character(g)
      g[g=='1/1']=2
      g[g=='0/1']=1
      g[g=='1/0']=1
      g[g=='0/0']=0
      g = as.numeric(g)
      return(g)
    }
    G = apply(G,2,recode)
    dimnames(G) = tags
    G = t(G)
    chr = as.vector(table(gsub('_.+','',tags[[1]])))
    final = list('gen'=G,'chr'=chr)
    return(final)
  }
  
  if(type=='VCF'){
    return(VCF(file))
  }
  
}
