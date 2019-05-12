#-----------------------------Estimating Parameter of mu Beta Distribution--------------------------
set.seed(20)
# The package below is used to calculate the parameter of beta distribution given the quantiles.
# Uncomment if the package is not installed.
# install.packages("TeachBayes")
library("TeachBayes")

lower_quantile <- list(x = 45/56, p = .25)
upper_quantile <- list(x = 51/56, p = .75)

(ab <- beta.select(lower_quantile, upper_quantile))

#-----------------------------Estimating Parameter of sigma Beta Distribution-----------------------

set.seed(20)
# The package below is used to calculate the parameter of beta distribution given the quantiles.
# Uncomment if the package is not installed.
# install.packages("TeachBayes")
library("TeachBayes")

lower_quantile <- list(x = 60/600, p = .25)
upper_quantile <- list(x = 90/600, p = .75)

(ab <- beta.select(lower_quantile, upper_quantile))
p = seq(0,1, length=100)
plot(p*600, dbeta(p, ab[1], ab[2]), ylab="density", type ="l", col=4, xlab = "Time (s)")

#-----------------------------Plotting Beta distribution--------------------------------------------
par(mfrow=c(1,2))

set.seed(20)
p = seq(0,1, length=100)

#par(mfrow=c(2,1))

# Normalised Beta distribution
#plot(p, dbeta(p, 7.65, 1.23), ylab="density", type ="l", col=4)

# Unnormalised Beta distribution

plot(p*3360, dbeta(p, 16.03, 2.82)/sum(dbeta(p, 16.03, 2.82)), ylab="density", 
     type ="l", col="red", xlab = "Time (s)", main = "Stretch Beta(16.03,2.82)")
# LQ
abline(v=45*60, col="black")
# UQ
abline(v=51*60, col="black")

p = seq(0,1, length=100)

plot(p*3360,dunif(p, min = 0, max = 56*60), ylab = "density", 
     xlab = "Time (s)", main = "Uniform", type = "l", col="blue")

## Sigma (standard deviation)
par(mfrow=c(1,3))

set.seed(20)

p = seq(0,1, length=100)
plot(600*p, dbeta(p, 10, 68.92)/sum(dbeta(p, 10, 68.92)), ylab="density",
     type ="l", col="red", xlim = c(-1,600), xlab = "Time(s)",main = "Beta(10,68.92)")
abline(v=60, col="black")
abline(v=90, col="black")
plot(p*600,dunif(p, min = 0, max = 600), xlim=c(-10,610),type = "l", 
     xlab = "Time(s)", ylab = "density",main = "Uniform", col="blue")
plot(p*600, dbeta(p,2,2), xlim=c(-10,610),type = "l", xlab = "Time(s)", 
     ylab = "density", main = "Beta(2,2)", col="green")
#-----------------------------Pre-Posterior--------------------------------------------------------

N=10000

pre_posterior = NULL

for(i in 1:N){
  
  # Sampling from prior
  # Prior for mu
  A = rbeta(1,16.03, 2.82)
  mu = A*56*60
  
  # Prior for sigma
  B = rbeta(1,2,2)
  sigma= 9*B+1
  
  pre_posterior[i]= -1000
  
  while ((pre_posterior[i]>56*60)|(pre_posterior[i]<0)) {
    pre_posterior[i]= rnorm(1,mu,sigma)
  }
}

plot(density(pre_posterior), xlab="Time(s)", main = "Preposterior")

#----------------------------------Loading-Data---------------------------------------------------

data = c(2900,2810,2990,2980,2930,2990,3180,3200,3050,2980,2710,3150)

#----------------------------------Posterior------------------------------------------------------

#-----------Posterior-Function----------
# Log-Posterior of my prior
mypostfunction = function(m,s){
  (-12/(2*s^2)) - 
    (1/(2*s^2))*((2900-m)^2 + 
                  (2810-m)^2 +
                  (2990-m)^2 +
                  (2980-m)^2 +
                  (2930-m)^2 +
                  (2990-m)^2 + 
                  (3180-m)^2 +
                  (3200-m)^2 +
                  (3050-m)^2 +
                  (2980-m)^2 +
                  (2710-m)^2 +
                  (3150-m)^2) - 
    12*log(s*(pnorm(3360,m,s)-pnorm(0,m,s))) + 
    15.03*log(m/3360) + 
    1.82*log(1-m/3360) + 
    9*log(s/600)+
    67.92*log(1-s/600)
}

# Log-Posterior using naive prior
naivepostfunction = function(m,s){
  (-12/(2*s^2)) - 
    (1/(2*s^2))*((2900-m)^2 + 
                   (2810-m)^2 +
                   (2990-m)^2 +
                   (2980-m)^2 +
                   (2930-m)^2 +
                   (2990-m)^2 + 
                   (3180-m)^2 +
                   (3200-m)^2 +
                   (3050-m)^2 +
                   (2980-m)^2 +
                   (2710-m)^2 +
                   (3150-m)^2) - 
    12*log(s*(pnorm(3360,m,s)-pnorm(0,m,s)))
}

mypostfunction2 = function(m,s){
  (-12/(2*s^2)) - 
    (1/(2*s^2))*((2900-m)^2 + 
                   (2810-m)^2 +
                   (2990-m)^2 +
                   (2980-m)^2 +
                   (2930-m)^2 +
                   (2990-m)^2 + 
                   (3180-m)^2 +
                   (3200-m)^2 +
                   (3050-m)^2 +
                   (2980-m)^2 +
                   (2710-m)^2 +
                   (3150-m)^2) - 
    12*log(s*(pnorm(3360,m,s)-pnorm(0,m,s))) + 
    15.03*log(m/3360) + 
    1.82*log(1-m/3360) + 
    log(s/600)+
    log(1-s/600)
}


#---------------Metropolis-Hasting's-Sampling-------

# Metropolis Algorithm for Truncated Normal with mylogpostfunction
set.seed(20)

N = 100000
tpost = matrix(0,N,2)
#we need to specify two starting values (one for each parameter) 
tpost[1,1:2] = c(1000,70)
for (i in 2:N){
  #propose a new mu value and a sigma value
  mustar = rnorm(1,tpost[i-1,1],120)
  sigstar = rnorm(1,tpost[i-1,2],60)
  #we now calculate our (log) acceptance ratio, but we should only 
  #do this for possible values of mu and sigma...
  if (mustar < 3360 & mustar > 0 & sigstar > 60 & sigstar < 600){ 
  rtop = mypostfunction(mustar,sigstar)
  rbot = mypostfunction(tpost[i-1,1],tpost[i-1,2])
  U = runif(1,0,1)
  if (log(U)<(rtop-rbot)){
    tpost[i,1:2] = c(mustar,sigstar)
  } else { tpost[i,1:2] = tpost[i-1,1:2] } } else { tpost[i,1:2] = tpost[i-1,1:2] }
}

# Posterior using naive prior
set.seed(20)

N = 100000
tpost_naive = matrix(0,N,2)
#we need to specify two starting values (one for each parameter) 
tpost_naive[1,1:2] = c(1000,70)
for (i in 2:N){
  #propose a new mu value and a sigma value
  mustar = rnorm(1,tpost[i-1,1],120)
  sigstar = rnorm(1,tpost[i-1,2],60)
  #we now calculate our (log) acceptance ratio, but we should only #do this for possible values of mu and sigma...
  if (mustar < 3360 & mustar > 0 & sigstar > 60 & sigstar < 600){ rtop = naivepostfunction(mustar,sigstar)
  rbot = naivepostfunction(tpost_naive[i-1,1],tpost_naive[i-1,2])
  U = runif(1,0,1)
  if (log(U)<(rtop-rbot)){
    tpost_naive[i,1:2] = c(mustar,sigstar)
  } else { tpost_naive[i,1:2] = tpost_naive[i-1,1:2] } } else { tpost_naive[i,1:2] = tpost_naive[i-1,1:2] }
}

# Metropolis Algorithm for Truncated Normal with a bad prior (beta(2,2))
set.seed(20)

N = 100000
tpost2 = matrix(0,N,2)
#we need to specify two starting values (one for each parameter) 
tpost2[1,1:2] = c(1000,70)
for (i in 2:N){
  #propose a new mu value and a sigma value
  mustar = rnorm(1,tpost2[i-1,1],120)
  sigstar = rnorm(1,tpost2[i-1,2],60)
  #we now calculate our (log) acceptance ratio, but we should only 
  #do this for possible values of mu and sigma...
  if (mustar < 3360 & mustar > 0 & sigstar > 60 & sigstar < 600){ 
  rtop = mypostfunction2(mustar,sigstar)
  rbot = mypostfunction2(tpost2[i-1,1],tpost2[i-1,2])
  U = runif(1,0,1)
  if (log(U)<(rtop-rbot)){
    tpost2[i,1:2] = c(mustar,sigstar)
  } else { tpost2[i,1:2] = tpost2[i-1,1:2] } } else { tpost2[i,1:2] = tpost2[i-1,1:2] }
}

# Metropolis-Hasting's for Binomial likelihood
set.seed(20)

N = 100000
tpost3 = NULL
#we need to specify two starting values (one for each parameter) 
tpost3[1] = 0.1
alpha = 16.03+12*mean(data/10) - 1
beta = 2.82 + 12*(336-mean(data/10)) - 1
for (i in 2:N){
  tstar = rbeta(1,1,1)
  top = alpha*log(tstar)+beta*log(1-tstar)
  bot = alpha*log(tpost3[i-1])+beta*log(1-tpost3[i-1])
  qtop = dbeta(tpost3[i-1],2,2)
  qbot = dbeta(tstar,2,2)
  U = runif(1,0,1)
  if(log(U)<top+qtop-bot-qbot){
    tpost3[i]=tstar
  }else{tpost3[i]=tpost3[i-1]}
}

#---------------------------------------Plot of Posterior Data----------------------------------

#-------Plot of sample data----

par(mfrow=c(2,3))
plot(tpost[,2], type = "l",ylab = "sigma", main ="Beta(10,68.92) (my) prior for sigma")
plot(tpost_naive[,2], type = "l",ylab = "sigma", main = "Uniform (naive) Prior for Sigma")
plot(tpost2[,2],type = "l",ylab = "sigma",main = "Beta(2,2) prior for Sigma")
plot(tpost[,1],type = "l",ylab = "mu", main = "Beta(16.03,1.03) (my) for mu")
plot(tpost_naive[,1],type = "l",ylab = "mu", main = "Uniform (naive) for mu")

#-------ACF-------------------

par(mfrow=c(2,2))
acf(tpost[,1],main="acf for mu (my prior)")
acf(tpost[,2],main="acf for sigma (my prior)")
acf(tpost_naive[,1], main="acf for mu (naive prior)")
acf(tpost_naive[,2], main="acf for sigma (naive prior)")

#---Plot of how different priors affect posterior---
par(mfrow=c(2,2))

p = seq(0,1, length=100)
plot(p*600,dbeta(p, 10, 68.92)/sum(dbeta(p, 10, 68.92)),col="red",type = "l",
     xlab = "Time(s)", ylab = "density", main = "Prior distribution of sigma")
lines(p*600,dunif(p,min=0,max=600),col="blue",type = "l")
lines(p*600,0.05*dbeta(p,2,2),col="green", type = "l")
legend("topright",
       legend = c("Beta(10,68.92) Prior", "Uniform Prior", "Beta(2,2)"),
       col = c("red","blue","green"),
       lty=1
)

sig = density(tpost[1000:10000,2])
sig_naive = density(tpost_naive[1000:10000,2])
sig$y = sig$y/sum(sig$y)
sig_naive$y = sig_naive$y/sum(sig_naive$y)
sig2 = density(tpost2[1000:10000,2])
sig2$y = sig2$y/sum(sig2$y)
plot(sig, col="red",xlim=c(0,600),ylim=c(0,0.01),xlab = "Time(s)", 
     main = "Distribution of Sigma", ylab = "Density")
lines(sig_naive,col="blue")
lines(sig2,col="green")
legend("topright",
       legend = c("My Posterior","Naive Posterior", "beta(2,2) Posterior"),
       col = c("red","blue","green"),
       lty=1
)

plot(3360*p,0.01*dbeta(p, 16.03, 2.82)/sum(dbeta(p, 16.03, 2.82)),col="red",type = "l", 
     xlab = "Time(s)", ylab = "density", main = "Prior distribution of mu")
lines(3360*p,dunif(p,min = 0,max = 3360), col="blue", type = "l")
legend("topleft",
       legend = c("My Prior", "Naive Prior"),
       col = c("red","blue"),
       lty=1
)

mu = density(tpost[1000:10000,1])
mu_naive = density(tpost_naive[1000:10000,1])
mu_naive$y = mu_naive$y/sum(mu_naive$y)
plot(mu, col="red",xlim=c(0,3360),ylim=c(0,0.015), xlab = "Time(s)", 
     main = "Distribution of mu", ylab = "Density")
lines(mu_naive,col="blue")
legend("topleft",
       legend = c("My Posterior", "Naive Posterior"),
       col = c("red","blue"),
       lty=1
)

#---------------Predictive-Distribution------------------

# Posterior using my prior and truncated normal

M = length(2000:N)
preH = NULL
postmu = tpost[2000:N,1]
postsig = tpost[2000:N,2]
for(i in 1:M){
  mu = postmu[i]
  sig = postsig[i]
  preH[i]=-1000
  while((preH[i]<0)|(preH[i]>3360)){
    preH[i]=rnorm(1,mu,sig)
  }
}

# Posterior using naive prior and truncated normal
M = length(2000:N)
naivepreH = NULL
postmu = tpost_naive[2000:N,1]
postsig = tpost_naive[2000:N,2]
for(i in 1:M){
  mu = postmu[i]
  sig = postsig[i]
  naivepreH[i]=-1000
  while((naivepreH[i]<0)|(naivepreH[i]>3360)){
    naivepreH[i]=rnorm(1,mu,sig)
  }
}

# Predictive Posterior using my prior and binomial model

M = length(2000:N)
binompreH = NULL
postmu = tpost3[2000:N]
for(i in 1:M){
  mu = postmu[i]
  binompreH[i]=-1000
  while((binompreH[i]<0)|(binompreH[i]>3360)){
    binompreH[i]=rbinom(1,336,mu)
  }
}

#--------------Plotting-Predictive-Posterior-------------
mypred = density(preH)
mypred$y = mypred$y/sum(mypred$y)

naivepred = density(naivepreH)
naivepred$y = naivepred$y/sum(naivepred$y)

prepost = density(pre_posterior)
prepost$y = prepost$y/sum(prepost$y) 

plot(mypred, xlab = "Time(s)", ylab = "density", main = "Predictive Posterior", col="red", xlim=c(700,3360), ylim=c(0,0.008))
lines(naivepred, col="blue")
lines(prepost, col="black")
legend("topleft",
       legend = c("My Pre. Posterior", "Naive Pre. Posterior", "Preposterior"),
       col = c("red","blue", "black"),
       lty=1
)
#---------------Different-Likelihood---------------------

# Plot of posterior parameter deriving from truncated normal and binormial likelihood.
par(mfrow=c(1,2))
plot(tpost[,1],type="l", main = "plot of mu from truncated normal", ylab = "mu")
plot(tpost3,type="l", main = "plot of mu from binomial", ylab = "mu")

#--------Predictive Poster of Truncated Normal vs Binomial---
binormpred = density(binompreH)
binormpred$x = binormpred$x*10
binormpred$y = binormpred$y/sum(binormpred$y)
plot(mypred, xlab = "Time(s)", ylab = "density", main = "Predictive Posterior", col="red", xlim=c(2300,3360), ylim=c(0,0.008))
lines(binormpred, col="blue", type = "l")
legend("topleft",
       legend = c("Truncated Normal", "Binomial"),
       col = c("red","blue"),
       lty=1
)


#----------------------------------Credible-Interval-------------------------------

# Naive 95% credible interval, which forms a 
# highest density interval (HDI) for a unimodal distribution
Credible = function(density){
  for(i in 1:length(density$y)){
    if(sum(density$y[1:i])>0.025){
      LQ = density$x[i]
      break
    }
  }
  for(i in 1:length(density$y)){
    if(sum(density$y[1:i])>0.975){
      UQ = density$x[i]
      break
    }
  }
  return(c(LQ,UQ))
}

# Credible for posterior and predictive distribution for binomial data generating
theta = density(tpost3[2000:100000])
theta$y = theta$y/sum(theta$y)
Credible(theta)
mean(theta$x)
3360*Credible(theta)
3360*mean(theta$x)
# Credible for posterior and predictive distribution for Truncated Normal data generating

mu = density(tpost[1000:10000,1])
mu$y = mu$y/sum(mu$y)
sig = density(tpost[1000:10000,2])
sig$y = sig$y/sum(sig$y)
mypred$y = mypred$y/sum(mypred$y)
Credible(mu)
Credible(sig)
Credible(mypred)
mean(mu$x)
mean(sig$x)
mean(mypred$x)
