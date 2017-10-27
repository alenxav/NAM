# Load data
data(met)
Obs = Obs[order(Obs$ID),]
Gen = data.matrix(Gen)
chr = data.frame(table(as.numeric(substring(colnames(Gen),3,4))))[,2]
fam = as.numeric(substring(rownames(Gen),6,7))

# Phenotype
Y = matrix(NA,nrow(Gen),3)
rownames(Y) = rownames(Gen)
colnames(Y) = 2013:2015

# 2013
y13 = Obs[which(Obs$Year==2013),'DTM']
names(y13) = as.character(Obs[which(Obs$Year==2013),'ID'])
Y[,1] = y13[rownames(Y)]
# 2014
y14 = Obs[which(Obs$Year==2014),'DTM']
names(y14) = as.character(Obs[which(Obs$Year==2014),'ID'])
Y[,2] = y14[rownames(Y)]
# 2015
y15 = Obs[which(Obs$Year==2015),'DTM']
names(y15) = as.character(Obs[which(Obs$Year==2015),'ID'])
Y[,3] = y15[rownames(Y)]

# Reduce the size of data for the example (skip for full-data example)
chr = round(chr/5)
subsample = seq(1,nrow(Y),5)
Gen = Gen[subsample,seq(1,ncol(Gen),5)]
Y = Y[subsample,]; fam = fam[subsample]

# Association
GxE_GWA = gwasGE(Y,Gen,fam,chr)

# Manhattan plot
par(mfrow=c(1,1))
plot(GxE_GWA,main="GWAS",pch=20)

# Plotting intercept and variance components
par(mfrow=c(2,2))
plot(GxE_GWA$PolyTest$vi,main="Allele Effect",pch=20,ylab="Mu")
plot(GxE_GWA$PolyTest$vi,main="GxE variance",pch=20,ylab="Var (GE)")
plot(GxE_GWA$PolyTest$vg,main="Genetic variance",pch=20,ylab="Var (G)")
plot(GxE_GWA$PolyTest$vg,main="Environmental variance",pch=20,ylab="Var (E)")
