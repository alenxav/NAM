
CV_NAM=function(y,gen,k=5,Seeds=1:5,IT=500,BI=100,cl=NULL){
  
  # Cross-validation function
  folds = function(Seed,y,gen,k,it,bi){
    
    require(NAM)
    N = nrow(gen)
    # Begin folds
    set.seed(Seed)
    Nk = round(N/k)
    w=sample(1:N,Nk)
    Y=y; y[w]=NA
    
    cat('\n CROSS-VALIDATION CYCLE \n\n')
    
    # BayesA
    cat('BayesA (wgr)\n')
    f1=wgr(y,gen,iv = T,it=IT,bi=BI,verb=T)
    # BayesB
    cat('BayesB (wgr)\n')
    f2=wgr(y,gen,iv = T,pi=0.8,it=IT,bi=BI,verb=T)
    # BayesC
    cat('BayesC (wgr)\n')
    f3=wgr(y,gen,pi=0.8,it=IT,bi=BI,verb=T)
    # BRR
    cat('BRR (wgr)\n')
    f4=wgr(y,gen,iv = T,it=IT,bi=BI,verb=T)
    # BRR
    cat('BLASSO (wgr)\n')
    f5=wgr(y,gen,de = T,it=IT,bi=BI,verb=T)
    # GBLUP
    cat('GBLUP (gmm)\n')
    f6=gmm(y,gen,model='GBLUP',it=IT,bi=BI)
    # RKHS
    cat('RKHS (gmm)\n')
    f7=gmm(y,gen,model='RKHS',it=IT,bi=BI)
    # Ridge
    cat('emRR\n')
    f8=emRR(y[-w],gen[-w,])
    # Laplace
    cat('emDE\n')
    f9=emDE(y[-w],gen[-w,])
    # Mix L1l2
    cat('emBL\n')
    f10=emBL(y[-w],gen[-w,])
    # Elastic net
    cat('emEN\n')
    f11=emEN(y[-w],gen[-w,])
    # Maximum likelihood
    cat('emML\n')
    f12=emML(y[-w],gen[-w,])
    
    NamesMod = c('BayesA','BayesB','BayesC',
                 'BRR','BLASSO','GBLUP','RKHS',
                 'emRR','emDE','emBL','emEN','emML'
                 'OBSERVATION')
    
    M = matrix(NA,Nk,length(NamesMod))
    colnames(M) = NamesMod
    for(i in 1:5) M[,i]=get(paste('f',i,sep=''))$hat[w]
    for(i in c(6:7)) M[,i]=get(paste('f',i,sep=''))$EBV[w]
    for(i in 8:12) M[,i]=gen[w,]%*%get(paste('f',i,sep=''))$b
    M[,13] = Y[w]
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
