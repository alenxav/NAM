
# Loading data
data(met)
Y = Obs[,6]
DTA = Obs[,1:5]

# Fitting model
fit = gmm(y=Y,gen=Gen,dta=DTA,200,100)
SP = fit$sp

# 2013
w = which(DTA$Block==1)
M1 = t(as.matrix(Matrix::sparseMatrix(x=SP[w],i=DTA$Row[w],j=DTA$Col[w])))
M1[M1==0]=NA
w = which(DTA$Block==2)
M2 = t(as.matrix(Matrix::sparseMatrix(x=SP[w],i=DTA$Row[w],j=DTA$Col[w])))
M2[M2==0]=NA
F1 = rbind(M1,NA,M2)
#####
w = which(DTA$Block==3)
M3 = t(as.matrix(Matrix::sparseMatrix(x=SP[w],i=DTA$Row[w],j=DTA$Col[w])))
M3[M3==0]=NA
w = which(DTA$Block==4)
M4 = t(as.matrix(Matrix::sparseMatrix(x=SP[w],i=DTA$Row[w],j=DTA$Col[w])))
M4[M4==0]=NA
F2 = rbind(M3,NA,M4)
F2013 = cbind(F1,NA,NA,NA,F2)

# 2014
w = which(DTA$Block==5)
M1 = t(as.matrix(Matrix::sparseMatrix(x=SP[w],i=DTA$Row[w],j=DTA$Col[w])))
M1[M1==0]=NA
w = which(DTA$Block==6)
M2 = t(as.matrix(Matrix::sparseMatrix(x=SP[w],i=DTA$Row[w],j=DTA$Col[w])))
M2[M2==0]=NA
F1 = rbind(M1,NA,M2)
#####
w = which(DTA$Block==7)
M3 = t(as.matrix(Matrix::sparseMatrix(x=SP[w],i=DTA$Row[w],j=DTA$Col[w])))
M3[M3==0]=NA
w = which(DTA$Block==8)
M4 = t(as.matrix(Matrix::sparseMatrix(x=SP[w],i=DTA$Row[w],j=DTA$Col[w])))
M4[M4==0]=NA
F2 = rbind(M3,NA,M4)
F2014 = cbind(F1,NA,NA,NA,F2)

# Picturing field
par(mfrow=c(2,1))
image(F2013,col=heat.colors(1000),main='Field variation - SoyNAM 2013',xaxt='n',yaxt='n')
image(F2014,col=heat.colors(1000),main='Field variation - SoyNAM 2014',xaxt='n',yaxt='n')
