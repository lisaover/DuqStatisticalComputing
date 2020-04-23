
#printPDF function receives one vector and a filename. It prints three plots for the vector to the specified filename: ts plot, hist, acf.
 
printPDF <- function(v1,filename) {
	
   pdf(filename)
    #PLOTS to PDF
  	par(mfrow=c(3,1)) #split plotting window into 3 rows
	ts.plot(v1,xlab="Iterations") 
	hist(v1,probability=T, cex.lab=1.5, cex.axis=1.5)
  	acf(v1, lag.max=1000)
   dev.off()
	
}


#check_alpha function receives three parameters "alpha.star" for the parameter probability, "param.star" for the next generated parameter, and "param0" for the previously accepted parameter.

check_alpha <- function(alpha.star, param.star, param0) {
	
	#if alpha.star is less than 1, draw a random uniform to compare
	#if random uniform draw is less than alpha star, return param.star and count=1
	#else return param0 and count=0
	if(alpha.star < 1) {
		if(runif(1,0,1) < alpha.star) {
			results <- list("param" = param.star, "count" = 1)
  			return(results) 
		}
		else {
			results <- list("param" = param0, "count" = 0)
  			return(results)
		}
	}
	#check if alpha star for parameter is equal to one, return parameter star and count=1
	else {
		results <- list("param" = param.star, "count" = 1)
  		return(results) 
	}
	
}

#The metro_hastings function receives six parameters: "x0" for the initial sample value, "k" for the interval, "N" for number of independent realizations, "lag" for determining how many realizations to skip between saves, and "burnin" for detemrining how many realizations to skip before starting to save.

metro_hastings <- function(x0,k,N,lag,burnin) {
		
	#Set N to be N*lag+burnin
	N <- N*lag + burnin

	#Initialize vector to hold x
	x.s <- NULL
	
	#store the acceptance rate for x in 'x.cnt'
	x.cnt = 0
	
  for(i in 1:N) {
	
	#generate x from uniform proposal density with (k) as the interval limits
	x.star <- runif(1,x0-k,x0+k)
	while((x.star < 0) | (x.star > 1)) {
		x.star <- runif(1,x0-k,x0+k)
	}
		
	#alpha star for x using uniform density as target density and hastings correction (proposal.ratio) to make density integrate to 1
	target.ratio = (1/(2*k))/(1/(2*k))
	proposal.ratio = max((1/(2*k)), (1/(x.star+k)), (1/(1-x.star+k)))/max((1/(2*k)), (1/(x0+k)), (1/(1-x0+k)))
	alpha.star = min(target.ratio*proposal.ratio, 1)
	
	#use check_alpha function to determine if x.star should be accepted
	results = check_alpha(alpha.star,x.star,x0)
	
	#if i is greater than burnin and if i is a multiple of the lag, store value
	#add results$count to total x.cnt--results$count will be 0 if x.star was not accepted and 1 ow
	#add accepted value to vector - x0 or x.star was returned as results$param
  if(i > burnin) {
  	if(i %% lag == 0) {
		x.cnt = x.cnt + results$count
		x.s = c(x.s, results$param)
    }
   }
   #set x0 equal to whatever check_alpha returned for next iteration
   x0 = results$param
	
   }
  filename = sprintf("Documents/R-FILES/HW8-hCorrect-lag-U-k-%s.pdf",k)
  printPDF(x.s,filename)

  vector <- list("x" = x.s, "xCount" = x.cnt)
  return(vector)
 }
 
#Set burnin=20, lag=100, x0=(random uniform value), k=0.2; run Metro with N=5000
N = 5000
lag = 100
burnin = 20
k = 0.2
x0 = runif(1,0,1)

vector = metro_hastings(x0,k,N,lag,burnin)

vector$xCount/N
#[1] 1 for metropolis without hastings correction
#[1] 0.9462 for metropolis-hastings

#########################################################################

#printPDF function receives two vectors and a filename. It prints three plots for each vector to the specified filename: ts plot, hist, acf.
 
printPDF <- function(v1,v2,alpha,lambda,filename) {
	
   pdf(filename)
    #PLOTS to PDF
  	par(mfrow=c(3,2)) #split plotting window into 3 rows and 2 columns
	ts.plot(v1,xlab="Iterations") 
	ts.plot(v2,xlab="Iterations") 
	hist(v1,probability=T, cex.lab=1.5, cex.axis=1.5)
	hist(v2,probability=T, cex.lab=1.5, cex.axis=1.5)
	yy = seq(10,70,by=.001)
	lines(yy, ((lambda^alpha)/gamma(alpha))*((1/yy)^(alpha+1))*exp(-lambda/yy))
  	acf(v1)
    acf(v2) 
   dev.off()
	
}

#metro_hastings function receives seven parameters: "data" for the data values of interest, "m0" for the population mean, "s0" for the population variance, "b" for the mu interval, "c" for the sig sqr interval, "N" for number of independent random normal realizations, "lag" for determining how many realizations to skip between saves, and "burnin" for detemrining how many realizations to skip before starting to save.

metro_hastings <- function(data,m0,s0,b,c,N,lag,burnin) {
	
	#obtain mean (y.bar), variance (s.sqr), and length (n) of data
	y.bar = mean(data)
    s.sqr = var(data)
    n = length(data)
	
	#Set N to be N*lag+burnin
	N <- N*lag + burnin

	#Initialize vectors to hold mu and sigma sqr
	m.s <- NULL
	s.s <- NULL
	
	#store the acceptance rate for mu in 'mu.cnt' and sigma sqr in 'sig.cnt'
	mu.cnt = 0
	sig.cnt = 0
	
  for(i in 1:N) {
  	
#Generate a mu star 'm.star' from the proposal density (Normal) using the population mu (m0) and (b) as the interval limits
	m.star = rnorm(1,m0,b)
	
#Compute alpha star for mu using mean (y.bar), variance (s.sqr), length (n) of data, and population variance (s0) -- Use target.ratio which is ratio of target density given m.star to the target density given m0 (target density is Normal)
	target.ratio = exp(-(1/(2*s0))*(n*(y.bar-m.star)^2 + m.star^2/10 - n*(y.bar-m0)^2 - m0^2/10))
	
#Use target.ratio*proposal.ratio to determine if s.star should be accepted
	save.param = 0
	accept.code = 0
		if(runif(1,0,1) < target.ratio) {
			save.param = m.star
			accept.code = 1
		}
		else {
			save.param = m0
			accept.code = 0
		}

	
#if i is greater than burnin and if i is a multiple of the lag, store value
#add accept.code to total mu.cnt - accept.code will be 0 if m.star was not accepted and 1 ow
#add accepted value - m0 or m.star as save.param
	if(i > burnin) {
  	if(i %% lag == 0) {
		mu.cnt = mu.cnt + accept.code
		m.s = c(m.s, save.param)
    }
   }
   #set m0 equal to save.param for next iteration
   m0 = save.param
	
#Generate a sigma sqr star 's.star' from the proposal density using the population variance (s0) and (c) as the interval limits
	s.star = rnorm(1,s0,c)
	while(s.star < 0) {
		s.star = rnorm(1,s0,c)
	}
	
	print("s0")
	print(s0)
	print("s.star")
	print(s.star)
		
#Compute alpha star for sig sqr using mean (y.bar), variance (s.sqr), length (n) of data, population variance (s0) and mean (m0)

	#Use target.ratio which is ratio of target density given s.star to the target density given s0 (target density is inverse gamma)
	A.star = (1/s.star)^((n/2)+4)
	B.star = (1/(2*s.star))*((n - 1)*s.sqr + n*(y.bar-m0)^2 + m0^2/10 + 2)
	A.zero = (1/s0)^((n/2)+4)
	B.zero = (1/(2*s0))*((n - 1)*s.sqr + n*(y.bar-m0)^2 + m0^2/10 + 2)
	target.ratio = exp(log(A.star) - B.star - log(A.zero) + B.zero)
	
	print("A.star")
	print(A.star)
	print("A.zero")
	print(A.zero)
	
	#Multiply target.ratio by Hastings correction (proposal.ratio) which is the ratio of the proposal s.star density given s0 to the proposal s0 density given s.star -- proposal density is Normal -- Divide proposal densities by proportion of curve that appears after the boundary (above zero) so that the density equals one when integrated
	density.star = dnorm(s0, s.star, c)/(1-pnorm(0, s.star, c))
	density0 = dnorm(s.star, s0, c)/(1-pnorm(0, s0, c))
	proposal.ratio = density.star/density0
	
	#Use target.ratio*proposal.ratio to determine if s.star should be accepted
	save.param = 0
	accept.code = 0
		if(runif(1,0,1) < target.ratio*proposal.ratio) {
			save.param = s.star
			accept.code = 1
		}
		else {
			save.param = s0
			accept.code = 0
		}

#if i is greater than burnin and if i is a multiple of the lag, store value
	#add accept.code to total sig.cnt - accept.code will be 0 if s.star was not accepted and 1 ow
	#add accepted value - s0 or s.star as save.param
  if(i > burnin) {
  	if(i %% lag == 0) {
		sig.cnt = sig.cnt + accept.code
		s.s = c(s.s, save.param)
    }
   }
   #set s0 equal to results$param for next iteration
   s0 = save.param
   
   }
  
  #IG parameters to superpose a line over sigma sqr distribution
  alpha = length(data)/2 + 3.5
  lambda = ((length(data)-1)*var(data))/2 + 2
  #alpha = length(data)/2 - 1
  #lambda = ((length(data)-1)*var(data))/2
  
  filename = sprintf("Documents/R-FILES/HW8-hastings-N-b-%s-c-%s.pdf",b,c)
  printPDF(m.s,s.s,alpha,lambda,filename)

  vectors <- list("mu" = m.s, "muCount" = mu.cnt, "sigsqr" = s.s, "sigCount" = sig.cnt)
  return(vectors)
 
 }
 
data2 = read.table(file.choose(), header=F)
attach(data2)
data2 = data2$V1

#Set burnin, lag, m0 of data, and s0; run metro to obtain N realizations
N = 10000
lag = 30
burnin = 0
b = .5
c = .5
m0 = mean(data2)
s0 = var(data2)

vectors = metro_hastings(data2,m0,s0,b,c,N,lag,burnin)

vectors$sigCount/N
vectors$muCount/N

#Hastings
vectors$sigCount/N
#[1] 0.427
vectors$muCount/N
#[1] 0.1095

yy = seq(0,20,by=.001)
alpha = length(data2)/2 + 3.5
lambda = ((length(data2)-1)*var(data2))/2 + 2

yy = seq(10,70,by=.001)
alpha = length(data2)/2 - 1
lambda = ((length(data2)-1)*var(data2))/2

hist(vectors$sigsqr,probability=T, cex.lab=1.5, cex.axis=1.5)
lines(yy, ((lambda^alpha)/gamma(alpha))*((1/yy)^(alpha+1))*exp(-lambda/yy))

xx = seq(.01,50,length=250)
plot(xx,(1/xx)^80)

#########################

[1] 0
[1] "astar"
[1] 2.146938e+80
[1] "azero"
[1] Inf
[1] "astar"
[1] 2.953445e-20
[1] "azero"
[1] Inf
[1] "astar"
[1] Inf
[1] "azero"
[1] Inf
Error in if (target.ratio * proposal.ratio < 1) { : 
  missing value where TRUE/FALSE needed
  
  log(Inf)
[1] Inf
exp(Inf - Inf) = "NaN"
Inf - Inf = "NaN"
  exp(log(Inf) - log(Inf))
[1] NaN
> exp(log(Inf) - log(8))
[1] Inf
> exp(log(8) - log(Inf))
