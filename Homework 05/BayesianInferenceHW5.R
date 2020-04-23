#Gibbs function receives five parameters: "data" for the data values of interest, "m" for the population mean, "N" for number of realizations, "lag" for determining how many realizations to skip between saves, and "burnin" for detemrining how many realizations to skip before starting to save. Gibbs generates N independent random normal values.
 
gibbs <- function(data,m,N,lag,burnin) {
	
	#obtain mean (y.bar), variance (s.sqr), and length (n) of data
	y.bar = mean(data)
    s.sqr = var(data)
    n = length(data)
	
	#Set N to be N*lag+burnin
	N <- N*lag + burnin

	#Initialize vectors to hold mu and sigma sqr
	m.s <- NULL
	s.s <- NULL
	
  for(i in 1:N) {
	
	#generate a sigma sqr s using mean (y.bar), variance (s.sqr), and length (n) of data
	s <- 1/rgamma(1,n/2,0.5*((n-1)*s.sqr+n*(y.bar-m)^2))
	#generate a mu m using mean (y.bar) and length (n) of data and sigma sqr (s) from last step
	m <- rnorm(1, y.bar, sqrt(s/n))
	
  #if i is greater than burnin and if i is a multiple of the lag, store m and s
  if(i > burnin) {
  	if(i %% lag == 0) {
		m.s <- c(m.s,m)
		s.s <- c(s.s,s)
    }
   }
  }
  vectors <- list("mu" = m.s, "sigsqr" = s.s)
  return(vectors)
 }
 
data = read.table(file.choose(), header=F)
attach(data)
data = data$V1

#Set burnin to 10 and leave lag=1 and run Gibbs with N=2500
N = 2500
lag = 1
burnin = 10 
m = 80

vectors = gibbs(data,m,N,lag,burnin)

#From Gibbs with N=500, lag=1 and burnin=0, run ts.plot and autocorrelation (acf) to determine values for lag and burnin
#vectors$mu represents Gibbs sampler realizations for mu
#vectors$sigsqr represents Gibbs sampler realizations for sigma sqr
par(mfrow=c(2,2)) #split plotting window into 2 rows and 2 columns
ts.plot(vectors$mu,xlab="Iterations") 
ts.plot(vectors$sigsqr,xlab="Iterations") 
hist(vectors$mu,probability=T, cex.lab=1.5, cex.axis=1.5)
hist(vectors$sigsqr,probability=T, cex.lab=1.5, cex.axis=1.5)
yy = seq(10,70,by=.001)
alpha = length(data)/2 -1
lambda = ((length(data)-1)*var(data))/2
lines(yy, ((lambda^alpha)/gamma(alpha))*((1/yy)^(alpha+1))*exp(-lambda/yy))

par(mfrow=c(2,1))
acf(vectors$mu)
acf(vectors$sigsqr) 

mt.s = NULL
vt.s = NULL
N = 2500
n = length(data)
s.sqr = var(data)
y.bar = mean(data)

for(i in 1:N) {
	v = 1/rgamma(1,(n/2-1),(n-1)*s.sqr/2)
    m = rnorm(1,y.bar,sqrt(v/n))
	mt.s = c(mt.s,m)
	vt.s = c(vt.s,v)
}

par(mfrow=c(2,2)) #split plotting window into 2 rows and 2 columns
ts.plot(mt.s,xlab="Iterations") 
ts.plot(vt.s,xlab="Iterations")
hist(mt.s,probability=T, cex.lab=1.5, cex.axis=1.5)
hist(vt.s,probability=T, cex.lab=1.5, cex.axis=1.5)
yy = seq(10,70,by=.001)
alpha = length(data)/2 -1
lambda = ((n-1)*s.sqr)/2
lines(yy, ((lambda^alpha)/gamma(alpha))*((1/yy)^(alpha+1))*exp(-lambda/yy))

quantile(vectors$mu, 0.025)
#9.2466
quantile(vectors$mu, 0.975)
#12.4249
quantile(vectors$sigsqr, 0.025)
#17.7869
quantile(vectors$sigsqr, 0.975)
#41.4715

quantile(vectors$mu, 0.005)
#8.5788
quantile(vectors$mu, 0.995)
#13.0917
quantile(vectors$sigsqr, 0.005)
#15.9805
quantile(vectors$sigsqr, 0.995)
#50.6254

quantile(mt.s, 0.025)
#9.18
quantile(mt.s, 0.975)
#12.3514
quantile(vt.s, 0.025)
#17.9174
quantile(vt.s, 0.975)
#43.9736

quantile(mt.s, 0.005)
#8.6776
quantile(mt.s, 0.995)
#12.9029
quantile(vt.s, 0.005)
#15.9959
quantile(vt.s, 0.995)
#50.4055 
