mixed = function(y,random=NULL,fixed=NULL,data=NULL,X=list(),alg=emML,maxit=10,Deregress=FALSE,...){
  
  # Get y vector
  if(!is.null(data)) y = data[[deparse(substitute(y))]]
  Vy = var(y,na.rm=TRUE)
  
  # Remove missing
  if(anyNA(y)){
    wna = which(is.na(y)); y = y[-wna]
    data = droplevels.data.frame(data[-wna,])}
  
  # Random effects
  if(!is.null(random)){
    rnd = attr(terms(random),"term.labels")
    LMB0 = LMB = sapply(data[rnd], function(x) mean(table(x)) )
    df0 = sapply(data[rnd],function(x) ifelse(is.factor(x),length(unique(x)),1))
    RND = TRUE }else{ RND = FALSE }
  
  # Fixed effects or intercept
  if(!is.null(fixed)){
    fxd = attr(terms(fixed),"term.labels")
    df = sum(sapply(data[fxd],function(x) ifelse(is.factor(x),length(unique(x)),1)))
    FIX = TRUE }else{ FIX = FALSE; df=1 }
  B = H = A = list(); n = length(y)
  
  # Fixed Categorical & Continuous Term
  FCT = function(y,x){
    if(is.factor(x)){
      b = tapply(y,x,mean,na.rm=T);
      b[is.na(b)] = 0
      g = b[x]
      e = y-g
    }else{
      b = cov(y,x,use='pair')/var(x,na.rm=T);
      b = c(b);
      g = x*b
      e = y-g }
    return(list(g=g,b=b,e=e))}
  
  # Random Categorical & Continuous Term
  RCT = function(y0,x,lmb0){
    y00 = y0
    x00 = x
    if( anyNA(y0) ){
      x = x[!is.na(y0)]
      y0 = y0[!is.na(y0)]
      }
    Mean = function(x) sum(x)/(length(x)+lmb0)
    if(is.factor(x)){
      b = tapply(y0,x,Mean);
      if(anyNA(b)) b[is.na(b)] = 0
      g = b[x00]
    }else{
      b = cov(y0,x,use='pair')/(var(x,na.rm=T)+lmb0/length(y0));
      b = as.vector(b)
      g = b*x00
    }
    # Deregress
    if(Deregress){
      slope = c(crossprod(y0,g)/(crossprod(g)))
      xb = g*slope
      b0 = y0-xb
      g = b0+xb
    }
    return(list(h=g,b=b,e=y00-g))
  }
  
  # Structured random
  gws = function(e,i,...){
    x = data[[i]]
    if(iter==1){ e00 = e }else{ e00 = e + H[[i]]}
    e00 = e00 - mean(e00,na.rm = TRUE)
    e0 = tapply(e00,x,mean,na.rm=TRUE)[rownames(X[[i]])]
    # Fit WGR
    comn = intersect(names(e0),rownames(X[[i]]))
    h = alg(e0[comn],X[[i]][comn,],...)
    fit = c(h$hat)
    names(fit) = comn
    # Deregress
    if(Deregress){
      g = c(crossprod(fit,e0[comn])/(crossprod(fit)))
      xb = fit*g
      b0 = e0-xb
      fit = b0+fit*g
    }
    # Output
    hh = e0
    hh[comn] = fit
    hhh = hh[as.character(x)]
    if(anyNA(hhh)) hhh[is.na(hhh)] = 0
    res = e00 - hhh
    OUT = list(g=fit,h=hhh,b=h$b,e=res)
    return(OUT)}
  
  # Fit (loop)
  mu = mean(y,na.rm=T)
  e = yc = y-mu
  B[['Intercept']] = mu
  
  # Fit (loop)
  R2c = 0.5
  pb = txtProgressBar(style = 3)
  for(iter in 1:maxit){
    
    # Fixed effects 
    if(FIX){
      for(i in fxd){
        if(iter>1) e = e+H[[i]]
          go = FCT(e,data[[i]])
          e = go$e
          B[[i]] = go$b
          H[[i]] = go$g
      }
    }
    
    # Random effects
    if(RND){
      # Coefficients
      for(i in rnd){
        # Structured
        if(i%in%ls(X)){
          go = gws(e,i,...)
          e = go$e
          B[[i]] = go$g
          H[[i]] = go$h
          A[[i]] = go$b
        }else{
          # Unstructured 
          if(iter>=1) e+H[[i]]
            go = RCT(e,data[[i]],LMB[i])
            e = go$e
            B[[i]] = go$b
            H[[i]] = go$h
        }
      }
      
      # VarComp & Lambda
      Error = mean(yc*e,na.rm=T)
      SSa = sapply(H, function(h) mean(yc*h,na.rm=T) )
      Vg = c(SSa,Error=Error)
      Vg[Vg<(0.01*Vy)] = 0.01*Vy
      Va = Vg[rnd]
      Ve = Vg['Error']
      LMB = Ve/Va; names(LMB) = names(Va)
    }
    
    # Print R2 and check convergence based on Ve
    setTxtProgressBar(pb,iter/maxit)
    R2 = round(1-Ve/Vy,6)
    if(abs(R2c-R2)==0) break()
    R2c = R2
  }
  close(pb)
  
  # Fit model
  fit = rep(B[['Intercept']],length(y))
  for(i in 1:length(H)) fit = fit+H[[i]]
  
  # Wrapping up
  Fitness = list(obs=y,hat=fit,res=e,fits=H)
  OUT = list(Fitness=Fitness,Coefficients=B)
  if(RND){
    wVar = Vg/sum(Vg)
    VC = list( VarComponents = round(Vg,6),
               VarExplained = round(wVar,6) )
    OUT[['VarComp']] = VC;
    if(length(X)>0) OUT[['Structure']] = A
  }  
  
  class(OUT) = 'mixed'
  return(OUT)}

#############################################################################################################

mtmixed = function(resp, random=NULL, fixed=NULL, data, X=list(), maxit=10, init=10, regVC=FALSE){
  
  # Get y matrix
  k = length(resp) 
  rownames(data) = paste('Obs',1:nrow(data),sep='')
  y = data[,resp]
  Vy = apply(y,2,var,na.rm=TRUE)
  
  # Remove missing - (if missing for all traits)
  if(anyNA(y)){
    wna = which(rowMeans(is.na(y))==1);
    if(length(wna)>0){ y = y[-wna,]
    data = droplevels.data.frame(data[-wna,]) }}
  
  # Centralized y
  yc = apply(y,2,function(x) x-mean(x,na.rm=T))
  
  # Random effects
  if(!is.null(random)){
    rnd = attr(terms(random),"term.labels")
    LMB = apply(!is.na(y),2, function(q)  sapply(data[rnd], function(x) mean(table(droplevels(x[q])))) )
    LMB = matrix(LMB,nrow=length(rnd),ncol=k,dimnames=list(rnd,resp))
    df0 = apply(!is.na(y),2, function(q) sapply(data[rnd],function(x) ifelse(is.factor(x),length(unique(x[q])),1)) )
    # df0 = apply(!is.na(y),2, function(q) sapply(data[rnd],function(x) ifelse(is.factor(x), mean(table(x[q])) ,1)) )
    RND = TRUE }else{ RND = FALSE }
  
  # Fixed effects or intercept
  if(!is.null(fixed)){
    fxd = attr(terms(fixed),"term.labels")
    df = sum(sapply(data[fxd],function(x) ifelse(is.factor(x),length(unique(x)),1)))
    FIX = TRUE }else{ FIX = FALSE; df=1 }
  n = colSums(!is.na(y))
  B = H = list() # Store coefficients
  
  # Fixed Categorical & Continuous Term
  FCT = function(y,x){
    if(is.factor(x)){
      b = tapply(y,x,mean,na.rm=T);
      b[is.na(b)] = 0
      g = b[x]
      e = y-g
    }else{
      b = cov(y,x,use='pair')/var(x,na.rm=T);
      b = as.vector(b);
      g = x*b
      e = y-g }
    return(list(g=g,b=b,e=e))}
  
  # Random Categorical & Continuous Term
  RCT = function(yy,x,lmb){
    update = function(y0,x,lmb0){
      y00 = y0
      x00 = x
      if( anyNA(y0) ){
        x = x[!is.na(y0)]
        y0 = y0[!is.na(y0)] }
      Mean = function(x) sum(x)/(length(x)+lmb0)
      if(is.factor(x)){
        b = tapply(y0,x,Mean);
        if(anyNA(b)) b[is.na(b)] = 0
        g = b[x00]
      }else{
        b = cov(y0,x,use='pair')/(var(x,na.rm=T)+lmb0/length(y0));
        b = as.vector(b)
        g = b*x00 }
      return(list(h=g,b=b,e=y00-g))
    }
    outs = list()
    for(q in resp) outs[[q]] = update( y0=yy[,q],x,lmb0 = lmb[q])
    return(outs)
  }
  
  # MV parameters for structured terms
  MV = lapply(X, function(X){
    LIST = list()
    LIST[['b']] = matrix(0,ncol(X),k)
    rownames(LIST[['b']]) = colnames(X)
    colnames(LIST[['b']]) = resp
    kk = matrix(0,k,k,dimnames = list(resp))
    LIST[['iG']] = LIST[['vb']] = kk
    diag(LIST[['vb']]) = Vy/ncol(X)
    LIST[['iG']] = solve(LIST[['vb']])
    LIST[['ve']] = 0.5*Vy
    names(LIST[['ve']]) = resp
    return(LIST) } )
  
  # Structured random
  gws = function(e,i){
    x = data[[i]]
    if(iter==1){ e00 = e }else{ e00 = e + H[[i]][as.character(x),] }
    e0 = apply(e00,2,function(e00) tapply(e00,x,mean,na.rm=TRUE) )
    # Add for debugging
    if(any(!rownames(X[[i]])%in%rownames(e0))){
      ww = which(!rownames(X[[i]])%in%rownames(e0))
      ww = rownames(X[[i]])[ww]
      mm = matrix(NA,length(ww),k,dimnames=list(ww,resp))
      e0 = rbind(e0,mm)}
    e0 = e0[rownames(X[[i]]),]
    e0 = apply(e0,2,function(x) x-mean(x,na.rm=T))
    h = mtgsru(Y=e0,X=X[[i]],b=MV[[i]]$b,vb=MV[[i]]$vb,ve=MV[[i]]$ve,iG=MV[[i]]$iG,maxit=init)
    unlist(lapply(h,anyNA)); head(h$b)
    fit = h$hat
    rownames(fit) = rownames(X[[i]])
    comp = list(b=h$b,vb=h$vb,ve=h$ve,iG=h$iG,h2=h$h2,MSx=h$MSx)
    h = fit[as.character(x),]
    # dereg
    dereg = function(x,y){
      b = c(cov(x,y,use='p')/var(x,na.rm = T))
      xb = x*b
      mu = mean(y-xb,na.rm=T)
      ft = mu+xb
      return(ft)}
    for(drg in 1:k) h[,drg] = dereg(h[,drg],e00[,drg])               
    res = e00 - h
    colnames(fit) = resp
    OUT = list(g=fit,h=h,e=res,comp=comp)
    return(OUT)}
  
  # Get ready for loop
  e = y; R2c = 0.5;
  pb = txtProgressBar(style = 3)
  
  # Fit (loop)
  for(iter in 1:maxit){
    
    # Intercept
    if(iter==1){
      mu = apply(e,2,mean,na.rm=T)
      for(q in resp) e[,q] = e[,q]-mu[q]
      B[['Intercept']] = mu
    }else{
      mu = apply(e,2,mean,na.rm=T)
      for(q in resp) e[,q] = e[,q]-mu[q]
      B[['Intercept']] = mu+B[['Intercept']]
    }
    
    # Fixed effects 
    if(FIX){
      for(i in fxd){
        if(iter==1){
          go = apply(e,2,FCT,x=data[[i]])
          e = sapply(go, function(x) x$e )
          B[[i]] = sapply(go, function(x) x$b )
          H[[i]] = sapply(go, function(x) x$g )
        }else{
          E = e+H[[i]]
          go = apply(E,2,FCT,x=data[[i]])
          e = sapply(go, function(x) x$e )
          B[[i]] = sapply(go, function(x) x$b )
          H[[i]] = sapply(go, function(x) x$g )
        }}
    }
    
    # Random effects
    if(RND){
      
      # Coefficients
      for(i in rnd){
        
        # Structured
        if( i%in%ls(X) ){
          go = gws(e,i)
          e = go$e
          B[[i]] = go$g
          H[[i]] = go$h
          MV[[i]]$b = go$comp$b
          MV[[i]]$vb = go$comp$vb
          MV[[i]]$iG = go$comp$iG
          MV[[i]]$ve = go$comp$ve
          MV[[i]]$h2 = go$comp$h2
          MV[[i]]$MSx = go$comp$MSx
          rownames(MV[[i]]$iG) = resp
          colnames(MV[[i]]$vb) = resp
          
        }else{
          
          # Unstructured 
          if(iter==1){
            go = RCT( yy = e,x = data[[i]],lmb = LMB[i,])
            e = sapply(go, function(x) x$e )
            B[[i]] = sapply(go, function(x) x$b )
            H[[i]] = sapply(go, function(x) x$h )
          }else{
            E = e+H[[i]]
            go = RCT( yy = E,x = data[[i]],lmb = LMB[i,])
            e = sapply(go, function(x) x$e )
            B[[i]] = sapply(go, function(x) x$b )
            H[[i]] = sapply(go, function(x) x$h )
          }
        }
      }
      
      # VarComp & Lambda
      Error = colSums(e*e,na.rm=T)
      SSa = sapply(H, function(h) colSums(h*h,na.rm=T) )
      SS = cbind(SSa,Error)
      SS[SS<0] = 0.01
      if( regVC ) SS = sqrt(SS)
      wVar = t(apply(SS,1,function(x) x/sum(x) ))
      Vg = wVar*Vy
      Va = Vg[,rnd]/t(df0)*n
      Ve = Vg[,'Error'] * (n-df)/n
      LMB = t(Ve/Va)
      LMB = matrix(LMB,nrow=length(rnd),ncol=k,dimnames=list(rnd,resp))
      
      
    }
    
    setTxtProgressBar(pb,iter/maxit)
    # Print R2 and check convergence
    R2 = round(1-Ve/Vy,6)
    # cat('Iter ',iter,':',sep='')
    # cat(' R2 =',round(R2,3),'\n')
    if(sum(abs(R2c-R2))==0 & iter>=3){break(); setTxtProgressBar(pb, 1 )} 
    R2c = R2
    
    
  }
  close(pb)
  
  # Fit model
  fit = sapply(B[['Intercept']],function(x) rep(x,nrow(y)) )
  for(i in 1:length(H)) fit = fit+H[[i]]
  rownames(fit) = rownames(data)
  
  # Wrapping up
  Fitness = list(obs=y,hat=fit,res=e,fits=H)
  OUT = list(Fitness=Fitness,Coefficients=B)
  if(RND){
    VC = list( VarComponents = round(Vg,4),
               VarExplained = round(wVar,3) )
    OUT[['VarComp']] = VC;
    if(length(X)>0) OUT[['Structure']] = MV
  }
  
  # Return
  class(OUT) = 'mixed'
  return(OUT)}

Rcpp::cppFunction(
  'SEXP mtgsru(NumericMatrix Y, NumericMatrix X,
  NumericMatrix b, NumericMatrix vb, NumericVector ve,
  NumericMatrix iG, int maxit = 50){
  
  double tol = 10e-8;
  Rcpp::Environment base("package:base");
  Rcpp::Function solve = base["solve"];
  
  int k = Y.ncol(), p = X.ncol(), n0 = X.nrow();
  NumericMatrix fit(n0,k),o(n0,k),y(n0,k),e(n0,k);
  
  for(int i=0; i<k; i++){
  o(_,i) = ifelse(is_na(Y(_,i)),0,1);
  y(_,i) = ifelse(is_na(Y(_,i)),0,Y(_,i));}
  
  NumericVector n = colSums(o);
  NumericMatrix xx(p,k), vx(p,k); double tmp;
  for(int i=0; i<p; i++){ for(int j=0; j<k; j++){
  xx(i,j) = sum(X(_,i)*X(_,i)*o(_,j));
  tmp = sum(X(_,i)*o(_,j))/n(j);
  vx(i,j) = xx(i,j)/n(j)-tmp*tmp;}}
  
  NumericVector MSx = colSums(vx);
  NumericMatrix rho(k,k),LHS(k,k);
  NumericVector b0(k),b1(k),vy(k),RHS(k);
  
  for(int i=0; i<k; i++){ e(_,i) = y(_,i)+0;
  vy(i) = sum(e(_,i)*e(_,i))/n(i);}
  iG = solve(vb); NumericMatrix bc(p,k);
  
  int numit = 0; double cnv = 1; while(numit<maxit){
  bc = b+0; for(int j=0; j<p; j++){ b0 = b(j,_); LHS = iG+0;
  for(int i=0; i<k; i++){ LHS(i,i) = iG(i,i)+(xx(j,i)/ve(i));
  RHS(i) = (sum(e(_,i)*X(_,j))+xx(j,i)*b0(i))/ve(i);}
  b1 = solve(LHS, RHS);  b(j,_) = b1;
  
  for(int i=0; i<k; i++){ e(_,i) = (e(_,i)-X(_,j)*(b1(i)-b0(i)))*o(_,i);} }
  for(int i=0; i<k; i++){ ve(i) = (sum(e(_,i)*y(_,i)))/(n(i)-1);}
  for(int i=0; i<n0; i++){ for(int j=0; j<k; j++){fit(i,j) = sum(X(i,_)*b(_,j));}}
  for(int i=0; i<k; i++){ for(int j=0; j<k; j++){
  vb(i,j) = (sum(fit(_,i)*y(_,j))+sum(fit(_,j)*y(_,i))) / ((n(i)*MSx(i))+(n(j)*MSx(j))); }}
  for(int i=0; i<k; i++){ vb(i,i) = vb(i,i)*1.01;} // Ridging the diagonal of vb
  
  iG = solve(vb); ++numit; cnv = sum(abs(bc-b)); if( cnv<tol ){break;}}
  for(int i=0; i<n0; i++){for(int j=0; j<k; j++){fit(i,j) = sum(X(i,_)*b(_,j));}}
  
  NumericVector h2(k); h2 = 1-ve/vy;
  return List::create(Named("b")=b, Named("hat")=fit, Named("e")=fit, Named("MSx")=MSx,
  Named("vb")=vb, Named("ve")=ve, Named("iG")=iG, Named("h2")=h2);}')
