##############################################
#### 2. Simulation study comparing ols and gls
##############################################
library(MASS)

#### Set up the simulation study
length <- 10
x <- seq(0,1, length.out=length)
y <- seq(0,1, length.out=length)
grid <- expand.grid(x,y)
colnames(grid) <- c("X", "Y")
n <- nrow(grid)
Dist <- as.matrix(dist(grid))
Sigma <- exp(-Dist)
Sigma.inv <- solve(Sigma)
beta <- 10
X <- matrix(rep(1,n))
X.ols1 <- solve(t(X) %*% X) %*% t(X)
X.ols2 <- solve(t(X) %*% X)
X.ml1 <- solve(t(X) %*% Sigma.inv %*% X) %*% t(X) %*% Sigma.inv
X.ml2 <- solve(t(X) %*% Sigma.inv %*%  X)


#### Run the simulation study
n.simulation <- 1000
results <- array(NA, c(n.simulation, 4))
colnames(results) <- c("beta ols", "beta ml", "ci ols", "ci ml")

     for(i in 1:n.simulation)
     {
     ## Generate data
     z <- beta + mvrnorm(n=1, mu=rep(0, length^2), Sigma=Sigma)
          
     ## Estimate and CI using OLS
     beta.ols <- X.ols1 %*% z
     beta.sd.ols <- sqrt(X.ols2)
     ci.ols <- c(beta.ols - 1.96 * beta.sd.ols, beta.ols + 1.96 * beta.sd.ols)
     results[i, 1] <- beta.ols
     results[i,3] <- as.numeric(beta>ci.ols[1] &beta<ci.ols[2])
          
     ## Estimate and CI using GLS
     beta.ml <- X.ml1 %*% z
     beta.sd.ml <- sqrt(X.ml2)
     ci.ml <- c(beta.ml - 1.96 * beta.sd.ml, beta.ml + 1.96 * beta.sd.ml)
     results[i, 2] <- beta.ml
     results[i,4] <- as.numeric(beta>ci.ml[1] &beta<ci.ml[2])
     }


#### Summarise the results
## bias
apply(results[ ,1:2], 2, mean) - beta

## RMSE
sqrt(apply((results[ ,1:2]-beta)^2, 2, mean))

## Coverage
100 * apply(results[ ,3:4], 2, mean)










#### Plot realisations of exponential and Gaussian correlation functions
length <- 50
x <- seq(0,1, length.out=length)
y <- seq(0,1, length.out=length)
grid <- expand.grid(x,y)
par(mfrow=c(2,1))
exp.grf <- grf(grid=grid, cov.model="exponential", cov.pars=c(1, 1))
gau.grf <- grf(grid=grid, cov.model="gaussian", cov.pars=c(1, 1))
image(exp.grf)
image(gau.grf)


mat.grf <- grf(grid=grid, cov.model="matern", kappa=1.5, cov.pars=c(1, 1))
image(mat.grf)
