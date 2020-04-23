#Exercise 1:  Generate realizations for x using M and MH sampling for uniform (0,1) target density
#Method 1: Metropolis sampling

k = .2
lag=1
burnin=0
N = 5000*lag + burnin
acceptx1 = 0
xvec1 = NULL

#generate initial x value
x = runif(1,0,1)

for(i in 1:N){
#propose new value from uniform proposal (x-k,x+k)
xstar= runif(1,x-k,x+k)
while(xstar>1 | xstar<0){
	xstar= runif(1,x-k,x+k)
	}

#compute alpha probability
alphastar = min(1, (1/(2*k))/(1/(2*k)))

#accept xstar with prob alphastar and set x = xstar if so before storing value; otherwise store current x
	if(runif(1,0,1)<alphastar){
		x = xstar
		}
		
		if(i>burnin & i%%lag==0){
			xvec1 = c(xvec1,x)
			if(x == xstar){
			acceptx1 = acceptx1 +1
			}
		}	
}


#Method 2
#set constants
k = .2
lag=1
burnin=0
N = 5000*lag+burnin
acceptx2 = 0
xvec2 = NULL


#generate initial x value
x = runif(1,0,1)

for(i in 1:N){
#propose new value from uniform proposal (x-k,x+k)
xstar= runif(1,x-k,x+k)
while(xstar>1 | xstar<0){
	xstar= runif(1,x-k,x+k)
	}

#compute MH adjustment factor
num= max(1/(2*k),1/(xstar+k),1/(1-xstar+k))
denom = max(1/(2*k),1/(x+k),1/(1-x+k))

#compute alpha probability
alphastar = min(1, (1/(2*k))/(1/(2*k))*(num/denom))

#accept xstar with prob alphastar and set x = xstar if so before storing value; otherwise store current x
if(runif(1,0,1)<alphastar){
		x = xstar
		}
		if(i>burnin & i%%lag==0){
			xvec2 = c(xvec2,x)
			if(x == xstar){
			acceptx2 = acceptx2 +1
			}
		 }	
}


#results for both methods
results = function(xvec1,acceptx1, xvec2,acceptx2){
#acceptance prob for x
probx1 =acceptx1/length(xvec1)
probx2 = acceptx2/length(xvec2)

#corresponding acf plot and histograms into pdf
pdf("MHplots.pdf")
par(mfrow= c(2,2))

acf(xvec1, lag.max = 200)
hist(xvec1, probability = T)
acf(xvec2, lag.max = 200)
hist(xvec2, probability = T)

dev.off()

return(info = list("probx1"=probx1,"probx2"=probx2))
}

#call function
results(xvec1,acceptx1, xvec2,acceptx2)

#inital acceptance prob with lags 60, 75 respectively
$probx1
[1] 1

$probx2
[1] 0.9512
#these probs should be the same before and after implementing the necessary lag correct?

#******************************************************************************************
#Exercise 2  Generate realizations for mu and sigma using MH sampling for normal target and proposal densities
#read in Y's from data file
y = read.table(file.choose(), header = F)

#find the sample mean and variance
vary = var(y$V1)
meany = mean(y$V1)


n = 155
b= .07
c= .07
acceptmu = 0
acceptsig = 0
sigvec = NULL
muvec = NULL
lag = 50  #based of conservative mu=60 sigma = 80 lag values
burnin= 0
N= 10000*lag + burnin
muzero = meany #initial value for mu
sigzero = vary #initial value for sigma

for(i in 1:N){
#use proposal density for mu to generate realization of mustar
mu = rnorm(1, muzero, b)

#generate mu ratio
jointmu = exp(-1/(2*sigzero)*(n*(meany-mu)^2+mu^2/10-n*(meany-muzero)^2-muzero^2/10))

alphamu = min(jointmu,1)

#accept mustar if with prob alphamu and set muzero = mu if so before storing value; otherwise store current muzero
if(runif(1,0,1)<alphamu){
		muzero = mu
		}
		if(i>burnin & i%%lag==0){
			muvec = c(muvec,muzero)
			if(mu == muzero){
				acceptmu = acceptmu +1
			}
		}

#use proposal density for sig to generate realization for sigstar
sig = rnorm(1, sigzero, c)
while(sig<0){
	sig = rnorm(1,sigzero,c)
	}
	
#generate sig ratio
astar = log((1/sig)^(n/2+4))
bstar = 1/(2*sig)*((n-1)*vary+n*(meany-mu)^2+mu^2/10+2)
azero = log((1/sigzero)^(n/2+4))
bzero = 1/(2*sigzero)*((n-1)*vary+n*(meany-mu)^2+mu^2/10+2)
jointsig = exp(astar-bstar-azero+bzero)

MHcor = dnorm(sigzero,sig,c)/dnorm(sig,sigzero,c)*((1-pnorm(0,sigzero,c))/(1-pnorm(0,sig,c)))

alphasig = min(jointsig*MHcor,1)

#accept sigstar if with prob alphasig and set sigzero = sig if so before storing value; otherwise store current sigzero
if(runif(1,0,1)<alphasig){
		sigzero = sig
		}
		if(i>burnin & i%%lag==0){
			sigvec = c(sigvec,sigzero)
			if(sigzero == sig){
			acceptsig = acceptsig +1
			}
			}
}

results = function(muvec, sigvec,acceptmu, acceptsig){
#acceptance prob for mu
pmu =acceptmu/10000

#acceptance for sigma
psig =acceptsig/10000

#corresponding trace plot and histograms into pdf (plots may need adjusted between case 1 and 2)
pdf("MHplots2.pdf")
par(mfrow= c(2,2))

acf(muvec, lag.max = 200)
acf(sigvec, lag.max = 200)
ts.plot(muvec, xlab = "Iteration")
ts.plot(sigvec, xlab = "Iteration")
hist(muvec, probability = T)
hist(sigvec, probability = T)
yy = seq(0.3,0.8, length = 250)
alpha = 155/2-1
beta = (155-1)*vary/2
lines(yy, beta^alpha/gamma(alpha)*yy^(-alpha-1)*exp(-beta/yy))

dev.off()

return(info = list("muprob" =pmu, "SigProb"= psig))
}

results(muvec, sigvec,acceptmu, acceptsig)

#accept prob before lag
$muprob
[1] 0.6646

$SigProb
[1] 0.6697

#accept prob after lag
$muprob
[1] 0.6605

$SigProb
[1] 0.6599
