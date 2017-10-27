
CV_NAM=function(y,gen,k=5,Seeds=1:5,IT=400,BI=100,cl=NULL){
  
  # Cross-validation function
  folds = function(Seed,y,gen,k,it,bi){
    
    require(NAM)
    N = nrow(gen)
    # Begin folds
    set.seed(Seed)
    Nk = round(N/k)
    w=sample(1:N,Nk)
    Y=y; y[w]=NA
    
    # BayesB
    cat('BayesB\n')
    f1=wgr(y,gen,iv = T,pi=0.8,verb=T,it=IT,bi=BI)
    # BayesC
    cat('BayesC\n')
    f2=wgr(y,gen,pi=0.8,verb=T,it=IT,bi=BI)
    # BEN
    cat('BEN\n')
    f3=ben(y,gen,it=IT,bi=BI,bag = 0.8)
    # RF
    cat('RF\n')
    f4=gmm(y,gen,model = 'RF',it=round(IT/100),bi=round(BI/100))
    # BRR
    cat('BRR\n')
    f5=gmm(y,gen,model = 'BRR',it=IT,bi=BI)
    # BayesA
    cat('BayesA\n')
    f6=gmm(y,gen,model = 'BayesA',it=IT,bi=BI)
    # GBLUP
    cat('GBLUP\n')
    f7=gmm(y,gen,model = 'GBLUP',it=IT,bi=BI)
    # RKHS
    cat('RKHS\n')
    f8=gmm(y,gen,model = 'RKHS',it=IT,bi=BI)
    # Ridge
    cat('RR\n')
    f9=emRR(y[-w],gen[-w,])
    # BSR
    cat('BSR\n')
    f10=emBA(y[-w],gen[-w,])
    # SSVS
    cat('BVS\n')
    f11=emBB(y[-w],gen[-w,])
    # Mixture
    cat('Mixture\n')
    f12=emBC(y[-w],gen[-w,])
    # Laplace
    cat('Laplace\n')
    f13=emDE(y[-w],gen[-w,])
    # Mix L1l2
    cat('MixL1L2\n')
    f14=emBL(y[-w],gen[-w,])
    # Elastic net
    cat('ElasticNet\n')
    f15=emEN(y[-w],gen[-w,])
    # BLASSO2
    cat('BLASSO-L2\n')
    f16=wgr(y,gen,it=IT,bi=BI,de=TRUE)
    
    NamesMod = c('wgr.BayesB','wgr.BayesC','ben',
                 'gmm.RF','gmm.BRR','gmm.BayesA','gmm.GBLUP','gmm.RKHS',
                 'emRR','emBA','emBB','emBC','emDE','emBL','emEN','wgr.BL',
                 'OBSERVATION')
    
    M = matrix(NA,Nk,length(NamesMod))
    colnames(M) = NamesMod
    for(i in 1:3) M[,i]=get(paste('f',i,sep=''))$hat[w]
    for(i in c(4:10)) M[,i]=get(paste('f',i,sep=''))$EBV[w]
    for(i in 11:16) M[,i]=gen[w,]%*%get(paste('f',i,sep=''))$b
    M[,17] = Y[w]
    return(M)
  }
  
  # Running cross-validations
  if(is.null(cl)){
    b = lapply(Seeds,FUN=folds,y=y,gen=gen,k=k,it=IT,bi=BI)
    names(b) = paste('CV_',1:length(b),sep='')
  }else{
    b = clusterApply(cl,Seeds,fun=folds,y=y,gen=gen,k=k,it=IT,bi=BI)
    names(b) = paste('CV_',1:length(b),sep='')
  }
  
  return(b)
}

CV_Check=function(cv){
  n = length(cv)
  m = ncol(cv$CV_1)
  dta = matrix(0,0,m)
  for(i in 1:n) dta = rbind(dta,cv[[i]])
  # functions
  MSE = function(A,B) mean((A-B)^2)
  TOPS = function(A,B,TOP=0.2) mean(which(B>quantile(B,1-TOP,na.rm = TRUE))%in%which(A>quantile(A,1-TOP,na.rm = TRUE)))
  BIAS = function(A,B) cov(A,B)/var(A)
  # summary
  PA = sort(cor(dta)[-m,m],decreasing = TRUE)
  Rank = sort(cor(dta,method = 'sp')[-m,m],decreasing = TRUE)
  MSPE = sort(apply(dta[,-m],2,MSE,B=dta[,m]))
  Top20 = sort(apply(dta[,-m],2,TOPS,B=dta[,m]),decreasing = TRUE)
  Bias = apply(dta[,-m],2,BIAS,B=dta[,m])
  Bias = Bias[order(as.matrix(dist(c(1,Bias)))[-1,1])]
  # Model choice
  o1 = o2 = o3 = o4 = o5 = 1:(m-1)
  names(o1) = names(PA); o1 = o1[order(names(o1))]
  names(o2) = names(PA); o2 = o2[order(names(o2))]
  names(o3) = names(PA); o3 = o3[order(names(o3))]
  names(o4) = names(PA); o4 = o4[order(names(o4))]
  names(o5) = names(PA); o5 = o5[order(names(o5))]
  O = o1+o2+o3+o4+o5
  ModelChoice = names(sort(O))[1:5]
  # output
  final = list('PA'=PA,'Rank'=Rank,
               'MSPE'=MSPE,'Bias'=Bias,
               'Top20'=Top20,'ModelChoice'= ModelChoice)
  return(final)
  
}
