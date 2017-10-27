# Load data
data(met)

# Response
Y = Obs[,6]

# ID, covariate (Year) and spatial parameters (Block-Row-Col)
DTA = Obs[,1:5]

# Fit model
fit = gmm(y=Y,gen=Gen,dta=DTA)

# Summary
cor(fit$hat,fit$obs)

# Plot
plot(fit$hat,fit$obs,col=as.numeric(Obs$Year)+1,pch=20)
legend('topleft',legend=2013:2015,col=2:4,pch=20)
