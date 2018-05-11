# Meta-Analysis
META = function(ByEnv,Phe,fam,ammi=1){  
  require(NAM)
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
  for(j in 1:E){ f = factor( FAM[,j] ,levF)
  p[[j]]= c(1,sort(unique(as.numeric(f[!is.na(f)]))))}  
  # Preparing inputs, output matrix and design matrices
  p = as.factor(unlist(p))
  e = c(); for(j in 1:E) e = c(e,rep(j,Vi[j]))
  e = factor(e)
  n = length(p)
  y = rep(NA,n)
  # Weights
  W = apply(Phe,2,function(y) tapply(y, fam, function(x) mean(!is.na(x)) ))
  W = rbind(s=1,W)
  W = as.vector(W)
  W = W[W!=0]  
  # Starting GxE gwas
  ee = as.numeric(e)
  pp = as.numeric(p)
  maxP = length(levels(p))  
  # Interaction term
  ammi2 = function(y,pc){
    yy = lm(y~p)$residuals
    ge = g2 = matrix(0,maxP,E)
    for(k in 1:length(y)) ge[pp[k],ee[k]] = y[k]
    for(k in 1:length(y)) g2[pp[k],ee[k]] = y[k]
    df = sum(colMeans(g2)!=0)+sum(rowMeans(g2)!=0)-3
    S = svd(ge,pc,pc)
    SVD = matrix(0,length(y),pc)
    for(o in 1:pc){ SS = tcrossprod(S$u[,o],S$v[,o])
    for(k in 1:length(y)) SVD[k,o] = SS[pp[k],ee[k]]}
    Load = S$d[1]/sum(S$d)
    HH = list(SVD=SVD,load=Load,df=df)
    return(HH)}  
  # LOOP  
  output = matrix(NA,m,6)
  colnames(output) = c('ld','l0','l1','lrt','df','pv')  
  pb=txtProgressBar(style=3)
  for(i in 1:m){
    # collect the regression coefficient of marker effect
    y[1:n] = NA
    for(j in 1:E) y[(cVi[j]+1):cVi[j+1]] = ByEnv[[j]]$PolyTest[i,-c(1:12)]
    y = as.matrix(as.numeric(y),ncol=1)
    hh = ammi2(y,ammi)
    u = hh$SVD
    ld = hh$load
    df = hh$df
    # anova/LRT
    y[y==0]=NA
    l1a = lm(y~p+u,weights = W);
    l0a = lm(y~p,weights = W);
    l0 = logLik(l0a)[1]
    l1 = logLik(l1a)[1]
    lrt = -2*(l0-l1)
    pv = (-log(dchisq(lrt,df),base=10))
    #pv = -log(anova(l0a,l1a)$`Pr(>F)`[2])
    as = c(ld,l0,l1,lrt,df,pv)
    output[i,]=as
    setTxtProgressBar(pb,i/m)}
  close(pb)  
  FINAL_OUTPUT = list(PolyTest = data.frame(output),MAP = ByEnv[[1]]$MAP, SNPs = ByEnv[[1]]$SNPs, Method = "MetaGWA")
  class(FINAL_OUTPUT) = "NAM" 
  return(FINAL_OUTPUT)}
