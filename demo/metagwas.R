# Load data
data(met)
Obs = Obs[order(Obs$ID),]
chr = data.frame(table(as.numeric(substring(colnames(Gen),3,4))))[,2]
Gen = data.matrix(Gen)

# 2013
y13 = Obs[which(Obs$Year==2013),'DTM']
ids13 = as.character(Obs[which(Obs$Year==2013),'ID'])
names(y13) = ids13
n13 = length(y13)
fam13 = as.numeric(substring(ids13,6,7))

# 2014
y14 = Obs[which(Obs$Year==2014),'DTM']
ids14 = as.character(Obs[which(Obs$Year==2014),'ID'])
names(y14) = ids14
n14 = length(y14)
fam14 = as.numeric(substring(ids14,6,7))

# Reduce the size of data for the example (skip for full-data example)
chr = round(chr/5)
Gen = Gen[,seq(1,ncol(Gen),5)]
s13 = seq(1,n13,5)
y13 = y13[s13]; fam13 = fam13[s13]; ids13 = ids13[s13]
s14 = seq(1,n14,5)
y14 = y14[s14]; fam14 = fam14[s14]; ids14 = ids14[s14]
  
# Association
GWA13 = gwas3(y13,Gen[ids13,],fam13,chr)
GWA14 = gwas3(y14,Gen[ids14,],fam14,chr)

# Meta-analysis
Meta = meta3(list(GWA13,GWA14))

# Plot
par(mfrow=c(3,1))
plot(GWA13,main='2013')
plot(GWA14,main='2014')
plot(Meta,main='Meta-analysis')