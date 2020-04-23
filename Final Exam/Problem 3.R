#Select file and read data
data = read.csv(file.choose(), header=TRUE)
attach(data)

x = X4
p = 0.5
b = 1
m = mean(x)
Y = m + b
N = 25000
lag = 1
burnin = 0

#Gibbs function receives six parameters: "x" for the sample data values, "p" for the probability that data value, xi, comes from distribution 1, "Y" for the parameter of distribution 1 (b+m with b=1), "N" for number of realizations, "lag" for determining how many realizations to skip between saves, and "burnin" for detemrining how many realizations to skip before starting to save. Gibbs generates N independent z vectors of indicator variables, N independent means from a Poisson distribution, and N independent probabilities that are each a probability of drawing from distribution 1.
 
gibbs <- function(x,p,Y,N,lag,burnin) {
	
	#obtain length of x
	n = length(x)
    	
	#Set N to be N*lag+burnin
	N <- N*lag + burnin
	
	#Initialize vectors to hold the z vector realizations, m2 realizations, and p realizations
	zs = NULL
	Ys = NULL
	ps = NULL
	
	for(i in 1:N) {
	
	#Calculate the vector of probabilities where each value represents the probabiliy that the corresponding value in the data vector y was drawn from distribution 2
	probz1 = (1-p)*exp(-1)/factorial(x)
	probz2 = p*Y^x*exp(-Y)/factorial(x)
	probz = probz1 + probz2
	pvec = probz1/probz
	
	#Calculate the z vector of latent indicator values where each value is either 0 or 1 to indicate which distribution the corresponding value in the data vector x is from - "0" indicates distribution 1 and "1" indicates distribution 2
	z = rbinom(n,1,pvec)
	
	#Compute a mean for the mixed model from a gamma distribution using the indicator variable vector z and the data vector x.
	xbarz = sum(x*z)/sum(z)
	Y = rgamma(1,(sum(x)-sum(x*z)+1),(n-sum(z)+1/10))
	
	#Compute a p for the probability that the data value comes from distribution 1
	p = rbeta(1,n-sum(z)+1,sum(z)+1)
	
	#if i is greater than burnin and if i is a multiple of the lag, store z, Y, and p
    if(i > burnin) {
  	  if(i %% lag == 0) {
		zs <- c(zs,z)
		Ys <- c(Ys,Y)
		ps <- c(ps,p)
      }
     }
}

    vectors <- list("zs" = zs, "Ys" = Ys, "ps" = ps)
    return(vectors)

}

vectors = gibbs(x,p,Y,N,lag,burnin)

#vectors$zbars represents Gibbs sampler realizations for the mean of the z vector
#vectors$Ys represents Gibbs sampler realizations for the mean+1 of distribution 2
#vectors$ps represents Gibbs sampler realizations for the probability that an x value comes from distribution 1
zs = vectors$zs
Ys = vectors$Ys
ps = vectors$ps

#Since Y = m + 1, subtract 1 from Ys values to create a vector of m values
ms = Ys - 1

par(mfrow=c(3,2)) #split plotting window into 2 rows and 2 columns
ts.plot(ms,xlab="Iterations") 
ts.plot(ps,xlab="Iterations") 
hist(ms,probability=T, cex.lab=1.5, cex.axis=1.5)
hist(ps,probability=T, cex.lab=1.5, cex.axis=1.5)
acf(ms,lag.max=500)
acf(ps, lag.max=500) 

#Convert z to an n column matrix
zmat = matrix(zs, ncol=length(x), byrow=TRUE)
#Get column means of matrix zmat
zimeans = colMeans(zmat)

plot(x,zimeans)

mean(ms)
#[1] 6.259425
mean(ps)
#[1] 0.3321361
quantile(ms, 0.025)
#5.092297
quantile(ms, 0.975)
#7.501003
quantile(ps, 0.025)
#0.2256441
quantile(ps, 0.975)
#0.450701
