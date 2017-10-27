reml=function(y,X=NULL,Z=NULL,K=NULL){
  
  anyNA = function(x) any(is.na(x))
  if(!is.matrix(y)){
    # Dealing with random effect
    if(is.null(K)&is.null(Z))stop("Either Z or K must be specified")
    N = length(y)
    if(!is.null(Z)&!is.null(K)){
      model = 'Mixed'
    }else{
      if(is.null(Z)) model = 'Kernel'
      if(is.null(K)) model = 'Ridge'
    }
    # Dealing with fixed effect matrix
    if(is.null(X)){X = matrix(1,N,1)}else{
      if(is.matrix(X)){if(nrow(X)!=N) stop("Fixed effect does not match dimensions of response variable")
      }else{if(class(X)=="formula"){X=model.matrix(X)}}}
    if(model=='Ridge'){
      if(class(Z)=="formula"){
        Q = paste(Z); Q[2]=paste(Q[2],'- 1'); Z=as.formula(Q)
        Z=model.matrix(Z)
        }
      K=tcrossprod(Z)
      }
    if(model=='Kernel'){
      Z=diag(ncol(K))
      K0=K
      }
    if(model=='Mixed'){
      if(class(Z)=="formula"){
        Q = paste(Z); Q[2]=paste(Q[2],'- 1'); Z=as.formula(Q)
        Z=model.matrix(Z)
      }
      K0=K
      K = crossprod(t(Z),K); K=tcrossprod(K,Z)
      }
    # Function starts here
    m = which(is.na(y)) # missing values
    if(any(m)){
      y=y[-m]
      x=X[-m,]
      K=K[-m,-m]
      z=Z[-m,]
    }else{
      x=X;z=Z
    }
    x=as.matrix(x)
    # Defining log-REML
    loglike=function(theta){
      lambda=exp(theta)
      logdt=sum(log(lambda*delta+1))
      h=1/(lambda*delta+1)
      yy=sum(yu*h*yu);yx=matrix(0,q,1)
      xx=matrix(0,q,q);for(i in 1:q){
        yx[i]=sum(yu*h*xu[,i])
        for(j in 1:q){xx[i,j]=sum(xu[,i]*h*xu[,j])}}
      loglike = -0.5*logdt-0.5*(n-q)*log(yy-t(yx)%*%solve(xx)%*%yx)-0.5*log(det(xx))
      return(-loglike)}
    fixed=function(lambda){
      h=1/(lambda*delta+1)
      yy=sum(yu*h*yu)
      yx=timesVec(yu,h,xu,q)
      xx=timesMatrix(xu,h,xu,q,q)
      beta=qr.solve(xx,yx)
      sigma2=(yy-t(yx)%*%solve(xx)%*%yx)/(n-q)
      sigma2 = as.numeric(sigma2)
      var=diag((chol2inv(xx))*sigma2)
      stderr=sqrt(var)
      return(c(beta,stderr,sigma2))}
    # Eigendecomposition of K
    qq=eigen(as.matrix(K),symmetric=T)
    delta=qq[[1]]; uu=qq[[2]]; q=max(ncol(x),1)
    n=ncol(K); vp=var(y)
    yu=t(uu)%*%y; xu=t(uu)%*%x
    theta=0
    # Finding lambda through optimization
    parm=optim(par=theta,fn=loglike,method="L-BFGS-B",lower=-10,upper=10)
    lambda=exp(parm$par)
     # Variance components and fixed effect coefficient
    parmfix=fixed(lambda)
    beta=parmfix[1:q]
    sd=parmfix[(q+1):(2*q)]
    B = cbind(beta,sd)
    Ve=parmfix[2*q+1]
    Vg=lambda*Ve
    h2=Vg/(Vg+Ve); VC = data.frame(Vg,Ve,h2)
    # Random effect coefficient and prediction
    if(is.matrix(x)){re=y-x%*%beta}else{re =(y-tcrossprod(x,beta))}
    # RKHS solution
    if(model=='Kernel'){
      if(any(m)){
        qq=eigen(as.matrix(K0),symmetric=T)
        delta=qq[[1]]; uu=qq[[2]];
        dr = which(cumsum(delta)/sum(delta)>0.98)[1]
        g = qr.solve(uu[-m,],re)
        g = g/(Ve/(delta*Vg)+1)
        U = uu%*%g
      }else{
        dr = which(cumsum(delta)/sum(delta)>0.98)[1]
        uu = uu[,1:dr]
        delta = delta[1:dr]
        g = qr.solve(uu,re)
        g = g/(Ve/(delta*Vg)+1)
        U = uu%*%g
      }
    }
    if(model=='Ridge'){
      ZZ = crossprod(z)
      diag(ZZ) = diag(ZZ)+(Ve/Vg)
      Zy = crossprod(z,re)
      U = solve(ZZ,Zy)
    }
    if(model=='Mixed'){
      ZZ = crossprod(z)
      diag(K0)=diag(K0)+1e-8
      ZZ = ZZ + chol2inv(K0)*(Ve/Vg)
      Zy = crossprod(z,y)
      U = solve(ZZ,Zy)
    }
    REML = list("VC"=VC,"Fixed"=B,"EBV"=U)
    
  }else{
    N = nrow(y)
    if(is.null(X)){X = matrix(1,N,1)}else{
      if(is.matrix(X)){if(nrow(X)!=N) stop("Fixed effect does not match dimensions of response variable")
      }else{if(class(X)=="formula"){X=model.matrix(X)}}}
    if(is.null(K)&is.null(Z))stop("Either Z or K must be specified")
    if(is.null(K)){
      if(class(Z)=="formula"){Z=model.matrix(Z)-1}
      V=tcrossprod(Z)}
    if(is.null(Z)){V=K}
    if(is.null(Z)!=T&&is.null(K)!=T){
      if(class(Z)=="formula"){Z=model.matrix(Z)-1}
      V=crossprod(t(Z),K);V=tcrossprod(V,Z)};
    Z=diag(N); K=V
    if(any(is.na(Y))){
      hh = apply(Y,1,function(x)mean(is.na(x)))
      if(any(hh==1)){w = which(hh==1); Y = Y[-w,]; Z = Z[-w,]; X = X[-w,]}
      impY = function(Y){
        t = ncol(Y)
        NAs = which(is.na(Y))
        impX = function(x){x[is.na(x)]=mean(x,na.rm=TRUE);return(x)}
        Y2 = apply(Y,2,impX)
        rM=rowMeans(Y2,na.rm=TRUE)
        cM=colMeans(Y2,na.rm=TRUE)
        eM=tcrossprod(rM,cM)
        Y2 = Y; Y2[NAs]=eM[NAs]
        U = svd(x = Y2); Y3=Y
        for(i in 1:t){
          a = lm(Y[,i]~U$u)
          b = predict(a,data.frame(U$u))
          Y3[,i] = b}
        return(Y3)}
      Y = impY(Y)}
   ECM=function(Y,X,Z,K){Y=t(Y);X=t(X)
   ECM1=function(ytl, xtl, Vgt,Vet,Bt, deltal){Vlt=deltal*Vgt+Vet; invVlt=solve(Vlt+diag(1e-06,d))
   return(list(Vlt=Vlt, gtl=deltal*Vgt%*%invVlt%*%(ytl-Bt%*%xtl), Sigmalt=deltal*Vgt-deltal*Vgt%*%invVlt%*%(deltal*Vgt)))}
   wrapperECM1=function(l){ytl=Yt[,l]; xtl=Xt[,l]; deltal=eigZKZt$values[l]
   return( ECM1(ytl=ytl, xtl=xtl, Vgt=Vgt,Vet=Vet,Bt=Bt, deltal=deltal))}
   Vgfunc=function(l){Vgl=tcrossprod(outfromECM1[[l]]$gtl)
   return((1/n)*(1/eigZKZt$values[l])*(Vgl + outfromECM1[[l]]$Sigmalt))}
   Vefunc=function(l){etl = Yt[,l] - Bt%*%Xt[,l] - outfromECM1[[l]]$gtl
   return((1/n)*((tcrossprod(etl)+ outfromECM1[[l]]$Sigmalt)))}
   if (sum(is.na(Y))==0){ N=nrow(K); KZt=tcrossprod(K,Z)
   ZKZt=Z%*%KZt; eigZKZt = eigen(ZKZt); n=nrow(ZKZt); d=nrow(Y)
   Yt = Y%*%eigZKZt$vectors; Xt = X%*%eigZKZt$vectors
   Vgt =cov(t(Y))/2; Vet =cov(t(Y))/2; XttinvXtXtt=t(Xt)%*%solve(tcrossprod(Xt))
   Bt=Yt%*%XttinvXtXtt; Vetm1=Vet
   repeat{ outfromECM1=lapply(1:n, wrapperECM1)
         Vetm1=Vet; Gt=sapply(outfromECM1, function(x) {cbind(x$gtl)})
         Bt = (Yt - Gt) %*% XttinvXtXtt
         listVgts = lapply(1:n,Vgfunc); Vgt=Reduce('+', listVgts)
         listVets = lapply(1:n,Vefunc); Vet=Reduce('+', listVets)
         convnum=abs(sum(diag(Vet - Vetm1)))/abs(sum(diag(Vetm1)))
         convcond=tryCatch({convnum<1e-06}, error=function(e){return(FALSE)})
         if(convcond){break}}
 HobsInv=solve(kronecker(ZKZt,Vgt)+kronecker(diag(n),Vet)+diag(1e-06,d*n))
 ehat=matrix(Y - Bt%*%X,ncol=1, byrow=F)
 HobsInve=HobsInv%*%ehat; varvecG=kronecker(K,Vgt)
 gpred=varvecG%*%(kronecker(t(Z),diag(d)))%*%HobsInve
 Gpred=matrix(gpred, nrow=nrow(Y), byrow=F); 
 Fx=t(Bt); #rownames(Fx)=paste("beta",0:(nrow(Bt)-1),sep="")
 EBV = t(Gpred); colnames(Fx)=colnames(EBV)=paste("trait",1:d,sep="")
 VC = list("Vg"=Vgt,"Ve"=Vet,"h2"=diag(Vgt/(Vgt+Vet)),"GenCor"=cov2cor(Vgt))
 return(list("Fixed"=Fx,"VC"=VC,"EBV"=EBV))}}
 REML = ECM(Y=y,X=X,Z=Z,K=K)}

 class(REML) = "reml"
return(REML)}

MCreml = function(y,K,X=NULL,MC=300,samp=300){
  
  anyNA = function(x) any(is.na(x))
  
  if(nrow(K)!=length(y)){stop("Kinship and response variable have incompatible dimensions")}
  if(ncol(K)!=length(y)){stop("Kinship and response variable have incompatible dimensions")}
  n = MC; t = samp
  
  m = which(is.na(y)) # missing values
  if(any(m)){
    y=y[-m]
    X=X[-m,]
    K=K[-m,-m]
  }
  
  if(length(y)<t) t = round(0.5*length(y))
  
  moda=function (x){
    it=5;ny=length(x);k=ceiling(ny/2)-1; while(it>1){
      y=sort(x); inf=y[1:(ny-k)]; sup=y[(k+1):ny]
      diffs=sup-inf; i=min(which(diffs==min(diffs)))
      M=median(y[i:(i+k)]); it=it-1}; return(M)}
  
  h2 = c(); Vg = c(); Ve = c()
  for(i in 1:n){
    R = sample(1:length(y),t)
    fit = reml(y[R],X=X,K=K[R,R])
    Vg = c( fit$VC[1], Vg )
    Ve = c( fit$VC[2], Ve )
    h2 = c( fit$VC[3], h2 )
  }
  
  H = unlist(h2); G = unlist(Vg); E = unlist(Ve)
  samples = cbind(G,E,H); rownames(samples) = 1:MC
  mode.Vg = moda(G); mode.Ve = moda(E); mode.h2 = moda(H)
  VC = c(mode.Vg,mode.Ve,mode.h2); names(VC) = c("Vg","Ve","h2")
  result = list("samples"=samples,"modes"=VC)
  return(result)
}
