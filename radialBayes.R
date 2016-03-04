#########################
## Gamma model code
# for fitting radial sap flux density profiles with new data
# Aaron Berdanier, aaron.berdanier@gmail.com

### FUNCTIONS:
library(mvtnorm)
library(tmvtnorm)
library(msm)

a.n <- function(x) as.numeric(x)

gmaUpdate <- function(s,xx,yy,ab,dd){
  # yy are observed from group
  # ab are parameters from group
  # dd are random effects for individuals
  prop <- rtmvnorm(1, ab, jumps[,,s], c(-Inf,0,0))    
  
  vnow <- ab[1] + ab[2]*log(ab[3]*xx) - ab[3]*xx + dd
  vnew <- prop[1] + prop[2]*log(prop[3]*xx) - prop[3]*xx + dd
  pnow <- sum(dnorm(yy,vnow,sqrt(sg),log=T)) + dmvnorm(ab,rep(0,3),diag(rep(1000,3)),log=T)
  pnew <- sum(dnorm(yy,vnew,sqrt(sg),log=T)) + dmvnorm(prop,rep(0,3),diag(rep(1000,3)),log=T)
  
  nz   <- length(ab) # no. to accept
  a    <- exp(pnew - pnow) # acceptance PR
  z    <- runif(nz,0,1)
  keep <- which(z < a,arr.ind=T)
  ac <- 0
  if(any(keep)){
    ab[keep] <- prop[keep]           
    ac   <- length(keep)        
  } 
  temp <- list(x=ab,accept=ac)
  temp
}

dUpdate <- function(zhat,sg,tg){
  V <- ( as.numeric(table(ii))/sg + 1/tg )^(-1)
  vh <- zhat/sg
  Vv <- V*vh
  rnorm(ni,Vv,V)
}

updateSigmas <- function(y,mu,s1,s2){
  u1 <- s1 + length(y)/2
  u2 <- s2 + .5*sum( (y - mu)^2 )
  1/rgamma(1,u1,u2)
}
###

### DATA:
# v is instantaneous flux velocity (g/m2/s)
# x is the depth into the xylem (m)
# ii is an index for tree ID
# oi is the observation ID (unique)
# gi is the group ID
# ni is number of unique individuals

#
z <- log(v)

# initial values
abg <- matrix(rep(c(3,0.3,30),times=3),ncol=3,byrow=T)
dg <- rep(0, ni)
sg <- 2
tg <- 2
# initial metropolis jumps
jumps <- array(NA,c(3,3,3))
for(i in 1:3) jumps[,,i] <- diag(c(0.5,0.05,10))  

# priors
s1 <- 2 # higher value = tighter prior
s2 <- sg*(s1-1)
#hist(1/rgamma(10000,s1,s2))
t1 <- t2 <- 1
#hist(1/rgamma(10000,t1,t2),breaks="fd",xlim=c(0,20))


### MODEL:

ng <- 6000 # number of gibbs sampler iterations
burn <- 1000 # number of samples to burn

# data holders
abgib <- array(NA,c(ng,3,3)) # parameters
dgib <- matrix(0,ng,ni) # individual random effects
tgib <- sgib <- numeric(ng) # sgib = observation variance, tgib = variance on random effects
ygib <- matrix(0,ng,length(ii)) # predicted observations

prog <- txtProgressBar(min=0, max=ng, char="*", style=3)

for(g in 1:ng){
  
  # for each group, get slope parameters
  abg <- t(sapply(1:3,
                  function(s) gmaUpdate(s,
                                        xx=x[which(gi==s)],
                                        yy=z[which(gi==s)],
                                        ab=abg[s,],
                                        dd=dg[ii][which(gi==s)])$x
  ))
  abgib[g,,] <- c(t(abg))
  
  # for each individual, get random effects
  ytmp <- z - abg[gi,1] - abg[gi,2]*log(abg[gi,3]*x) + abg[gi,3]*x
  yin <- a.n(tapply(ytmp, ii, sum))
  dg <- dUpdate(yin, sg, tg)
  
  # for all observations, get uncertainty
  ytmp <- abg[gi,1] + abg[gi,2]*log(abg[gi,3]*x) - abg[gi,3]*x
  yhat <- ytmp+dg[ii]
  sg <- updateSigmas(z, yhat, s1, s2) 
  
  # update variance on random effects
  tg <- updateSigmas(dg, rep(0,length(dg)), t1, t2)
  
  dgib[g,] <- dg
  sgib[g] <- sg 
  tgib[g] <- tg
  
  ygib[g,] <- exp(yhat)  
  
  if(g %in% c(500,1000)) jumps[1:3,,] <- c(sapply(1:3,function(s) cov(abgib[1:g,,s])*1.9))
  
  setTxtProgressBar(prog, g)
}
close(prog)

# cut burns
abgib <- abgib[-c(1:burn),,] # posterior parameters
dgib <- dgib[-c(1:burn),] # posterior random effects
ygib <- ygib[-c(1:burn),] # posterior data predictions
sgib <- sgib[-c(1:burn)] # posterior observation variance
tgib <- tgib[-c(1:burn)] # posterior random effect variance