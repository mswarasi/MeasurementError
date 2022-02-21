
##########################################
# Measurement error models for gt.
# Try normal mixture latent models.
# Date: 08/28/2021
##########################################

rm( list=ls() )

##########################################

# Existing method:
loglik1 <- function(beta,Z,x,id,J,Se,Sp){
  p <- g( cbind(1,x)%*%beta )
  ll <- rep(-9,J)
  for(j in 1:J){
    pj <- p[(id[j]+1):id[j+1]]
    muj <- Se[j]-(Se[j]+Sp[j]-1)*prod(1-pj)
	ll[j] <- Z[j]*log(muj) + (1-Z[j])*log(1-muj)
  }
  return( -sum(ll) )
}

# Proposed method:
loglik.normal.known <- function(beta,mux,sdx,sdu,Z,w,id,J,Se,Sp){
  sd.xstar <- sdu*sdx/sqrt(sdu^2 + sdx^2)
  b <- beta[2]*sd.xstar
  ll <- rep(-9,J)
  for(j in 1:J){
    wj <- w[(id[j]+1):id[j+1]]
	mu.xstar <- (wj*sdx^2 + mux*sdu^2)/( sdx^2 + sdu^2 )
	a <- beta[1]+beta[2]*mu.xstar
	T1 <- (1-Se[j]+Z[j]*(2*Se[j]-1))*prod(dnorm(wj,mux,sqrt(sdx^2+sdu^2)))
	T2 <- (1-2*Z[j])*(Se[j]+Sp[j]-1)*prod(dnorm(wj,mux,sqrt(sdx^2+sdu^2))*(1-pnorm(a/sqrt(1+b^2))))
	ll[j] <- log( T1 + T2 )
  }
  return( -sum(ll) )
}  # Perfect!!!

# Proposed method:
loglik.mixture.known <- function(beta,mu1,sd1,mu2,sd2,r,sdu,Z,w,id,J,Se,Sp){
  sd.x1star <- sdu*sd1/sqrt(sdu^2 + sd1^2)
  sd.x2star <- sdu*sd2/sqrt(sdu^2 + sd2^2)
  b1 <- beta[2]*sd.x1star
  b2 <- beta[2]*sd.x2star
  ll <- rep(-9,J)
  for(j in 1:J){
    wj <- w[(id[j]+1):id[j+1]]
	mu.x1star <- (wj*sd1^2 + mu1*sdu^2)/( sd1^2 + sdu^2 )
	mu.x2star <- (wj*sd2^2 + mu2*sdu^2)/( sd2^2 + sdu^2 )
	a1 <- beta[1]+beta[2]*mu.x1star
	a2 <- beta[1]+beta[2]*mu.x2star
	K1 <- (1-Se[j]+Z[j]*(2*Se[j]-1))
	K2 <- (1-2*Z[j])*(Se[j]+Sp[j]-1)
	T1 <- K1*prod(r*dnorm(wj,mu1,sqrt(sd1^2+sdu^2))+(1-r)*dnorm(wj,mu2,sqrt(sd2^2+sdu^2)))
	T2 <- K2*prod(r*dnorm(wj,mu1,sqrt(sd1^2+sdu^2))*(1-pnorm(a1/sqrt(1+b1^2))) 
	        + (1-r)*dnorm(wj,mu2,sqrt(sd2^2+sdu^2))*(1-pnorm(a2/sqrt(1+b2^2))))
	ll[j] <- log(T1 + T2)
  }
  return( -sum(ll) )
}  # Perfect!!!

# Proposed method:
loglik.mixture.unknown <- function(param,r,sdu,Z,w,id,J,Se,Sp){
  beta <- param[1:2]
  mu1 <- param[3]
  sd1 <- param[4]
  mu2 <- param[5]
  sd2 <- param[6]

  sd.x1star <- sdu*sd1/sqrt(sdu^2 + sd1^2)
  sd.x2star <- sdu*sd2/sqrt(sdu^2 + sd2^2)
  b1 <- beta[2]*sd.x1star
  b2 <- beta[2]*sd.x2star
  ll <- rep(-9,J)
  for(j in 1:J){
    wj <- w[(id[j]+1):id[j+1]]
	mu.x1star <- (wj*sd1^2 + mu1*sdu^2)/( sd1^2 + sdu^2 )
	mu.x2star <- (wj*sd2^2 + mu2*sdu^2)/( sd2^2 + sdu^2 )
	a1 <- beta[1]+beta[2]*mu.x1star
	a2 <- beta[1]+beta[2]*mu.x2star
	K1 <- (1-Se[j]+Z[j]*(2*Se[j]-1))
	K2 <- (1-2*Z[j])*(Se[j]+Sp[j]-1)
	T1 <- K1*prod(r*dnorm(wj,mu1,sqrt(sd1^2+sdu^2))+(1-r)*dnorm(wj,mu2,sqrt(sd2^2+sdu^2)))
	T2 <- K2*prod(r*dnorm(wj,mu1,sqrt(sd1^2+sdu^2))*(1-pnorm(a1/sqrt(1+b1^2))) 
	        + (1-r)*dnorm(wj,mu2,sqrt(sd2^2+sdu^2))*(1-pnorm(a2/sqrt(1+b2^2))))
	ll[j] <- log(T1 + T2)
  }
  return( -sum(ll) )
}  # Perfect!!!

##########################################

# Proposed method:
loglik.mixture.unknown.all <- function(param,sdu,Z,w,id,J,Se,Sp){
  beta <- param[1:2]
  mu1 <- param[3]
  sd1 <- param[4]
  mu2 <- param[5]
  sd2 <- param[6]
  gma <- param[7]
  r <- exp(gma)/( 1 + exp(gma) )

  sd.x1star <- sdu*sd1/sqrt(sdu^2 + sd1^2)
  sd.x2star <- sdu*sd2/sqrt(sdu^2 + sd2^2)
  b1 <- beta[2]*sd.x1star
  b2 <- beta[2]*sd.x2star
  ll <- rep(-9,J)
  for(j in 1:J){
    wj <- w[(id[j]+1):id[j+1]]
	mu.x1star <- (wj*sd1^2 + mu1*sdu^2)/( sd1^2 + sdu^2 )
	mu.x2star <- (wj*sd2^2 + mu2*sdu^2)/( sd2^2 + sdu^2 )
	a1 <- beta[1]+beta[2]*mu.x1star
	a2 <- beta[1]+beta[2]*mu.x2star
	K1 <- (1-Se[j]+Z[j]*(2*Se[j]-1))
	K2 <- (1-2*Z[j])*(Se[j]+Sp[j]-1)
	T1 <- K1*prod(r*dnorm(wj,mu1,sqrt(sd1^2+sdu^2))+(1-r)*dnorm(wj,mu2,sqrt(sd2^2+sdu^2)))
	T2 <- K2*prod(r*dnorm(wj,mu1,sqrt(sd1^2+sdu^2))*(1-pnorm(a1/sqrt(1+b1^2))) 
	        + (1-r)*dnorm(wj,mu2,sqrt(sd2^2+sdu^2))*(1-pnorm(a2/sqrt(1+b2^2))))
	ll[j] <- log(T1 + T2)
  }
  return( -sum(ll) )
} ## Under inspection!


##########################################
# 3-component mixture.
# 9/9/2021
#
k.mix.loglik <- function(param,sdu,Z,w,id,J,Se,Sp){
  beta <- param[1:2]
  mu1 <- param[3]
  sd1 <- param[4]
  mu2 <- param[5]
  sd2 <- param[6]
  mu3 <- param[7]
  sd3 <- param[8]
  gma1 <- param[9]
  gma2 <- param[10]
  r1 <- exp(gma1)/( 1 + exp(gma1) )
  r2 <- exp(gma2)/( 1 + exp(gma2) )
  r3 <- 1 - r1 - r2

  sd.x1star <- sdu*sd1/sqrt(sdu^2 + sd1^2)
  sd.x2star <- sdu*sd2/sqrt(sdu^2 + sd2^2)
  sd.x3star <- sdu*sd3/sqrt(sdu^2 + sd3^2)

  b1 <- beta[2]*sd.x1star
  b2 <- beta[2]*sd.x2star
  b3 <- beta[2]*sd.x3star

  ll <- rep(-9,J)
  for(j in 1:J){
    wj <- w[(id[j]+1):id[j+1]]
	mu.x1star <- (wj*sd1^2 + mu1*sdu^2)/( sd1^2 + sdu^2 )
	mu.x2star <- (wj*sd2^2 + mu2*sdu^2)/( sd2^2 + sdu^2 )
	mu.x3star <- (wj*sd3^2 + mu3*sdu^2)/( sd3^2 + sdu^2 )

	a1 <- beta[1] + beta[2]*mu.x1star
	a2 <- beta[1] + beta[2]*mu.x2star
	a3 <- beta[1] + beta[2]*mu.x3star

	K1 <- (1-Se[j]+Z[j]*(2*Se[j]-1))
	K2 <- (1-2*Z[j])*(Se[j]+Sp[j]-1)

	T1 <- K1*prod(  r1*dnorm(wj,mu1,sqrt(sd1^2+sdu^2))
			      + r2*dnorm(wj,mu2,sqrt(sd2^2+sdu^2))
	              + r3*dnorm(wj,mu3,sqrt(sd3^2+sdu^2))
				  )

	T2 <- K2*prod( r1*dnorm(wj,mu1,sqrt(sd1^2+sdu^2))*(1-pnorm(a1/sqrt(1+b1^2))) 
	             + r2*dnorm(wj,mu2,sqrt(sd2^2+sdu^2))*(1-pnorm(a2/sqrt(1+b2^2))) 
				 + r3*dnorm(wj,mu3,sqrt(sd3^2+sdu^2))*(1-pnorm(a3/sqrt(1+b3^2)))
				 )

	ll[j] <- log(T1 + T2)
  }
  return( -sum(ll) )
} ## Under inspection!

##########################################

# Probit model:
g <- function(t) pnorm(t)

##########################################

## Simulation

J <- 400
cvec <- rep(5,J)
N <- sum(cvec)
id <- cumsum( c(0,cvec) )
beta.t <- c(-2, 1)
Se <- rep( .95, J)
Sp <- rep( .98, J)

sdu <- sqrt(.25)
mu1 <- 2.35
sd1 <- sqrt(.41)
mu2 <- -.26
sd2 <- sqrt(.38)
mu3 <- 0
sd3 <- 0.1

gma <-  -2.2
gma1 <- -2.2
gma2 <- -4.2

r <-  exp(gma)/( 1 + exp(gma) )  ## r = 0.1
r1 <- exp(gma1)/( 1 + exp(gma1) )  ## r = 0.1
r2 <- exp(gma2)/( 1 + exp(gma2) )  ## r = 0.1
r3 <- 1 - r1 - r2

param.t <- c(beta.t,mu1,sd1,mu2,sd2,mu3,sd3,gma1,gma2)
param0 <- param.t + .3

sims <- 1
se1 <- res1 <- matrix(-9,sims,2)
se2 <- res2 <- matrix(-9,sims,2)
se3 <- res3 <- matrix(-9,sims,length(param0))
se4 <- res4 <- matrix(-9,sims,length(param0))
conv3 <- conv4 <- rep(-9,sims)

## set.seed(321)

for(s in 1:sims){

  # Data simulation:
  u <- rnorm(N, 0, sdu)
# x <- rnorm(N, mux, sdx)

# Generate data from mixture:
  ru <- runif(N)
  x <- rep(-9,N)
  for(i in 1:N){
    if(ru[i] < r) x[i] <- rnorm(1,mu1,sd1)
    else x[i] <- rnorm(1,mu2,sd2)
  }
#  mean(x); sd(x)
#  hist(x)

  w <- x + u
  p.t <- g( cbind(1,x)%*%beta.t )
  Ytil <- rbinom(N,1,p.t)
  Z <- rep(-9, J)
  for(j in 1:J){
    Yj <- Ytil[(id[j]+1):id[j+1]]
    Z[j] <- rbinom(1,1,ifelse(sum(Yj)>0,Se[j],1-Sp[j]))
  }

# Model fitting:
  out1 <- optim( par=beta.t,fn=loglik1,Z=Z,x=w,id=id,J=J,Se=Se,Sp=Sp,hessian=FALSE,method="BFGS" )
  res1[s, ] <- out1$par
#  se1[s, ] <- sqrt(diag(solve(out1$hessian)))

#  out3 <- optim( par=param0,fn=loglik.mixture.unknown.all,sdu=sdu,Z=Z,w=w,id=id,J=J,Se=Se,Sp=Sp,hessian=TRUE,method="BFGS" )
#  res3[s, ] <- out3$par
#  se3[s, ] <- sqrt(diag(solve(out3$hessian)))
#  conv3[s] <- out3$convergence

  out4 <- optim( par=param0,fn=k.mix.loglik.temp,sdu=sdu,Z=Z,w=w,id=id,J=J,Se=Se,Sp=Sp,hessian=FALSE,method="BFGS" )
  res4[s, ] <- out4$par
#  se4[s, ] <- sqrt(diag(solve(out4$hessian)))
  conv4[s] <- out4$convergence

  print(s)
}

# Existing:
colMeans(res1)
# colMeans(se1)
apply(res1, 2, sd)

# Proposed:
colMeans(res4)
# colMeans(se4)
apply(res4, 2, sd)

## True
param.t


