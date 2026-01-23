
# Loading data
data(met)
Y = Obs[,6]
DTA = Obs[,1:5]

# Fitting model
fit = gmm(y=Y,gen=Gen,dta=DTA,200,100)
SP = fit$sp

# 2013, block 1
w = which(DTA$Block==1)
tmp = data.frame(x=round(SP[w],2),i=DTA$Row[w],j=DTA$Col[w])
tmp$x[tmp$x==0]=NA

# reshape tmp with i as rows and j as columns, filled with x
X = reshape(tmp,direction = 'wide',idvar='j',timevar='i',v.names='x')
rownames(X) = X$j
X = data.matrix(X[order(X$j),-1])
colnames(X) = gsub('x\\.','',colnames(X))
X = X[,as.character(1:ncol(X))]

# Plot field variation
image(X,
      main='SoyNAM field variation, IN 2013 - Block 1',
      col=heat.colors(1000),xaxt='n',yaxt='n')

