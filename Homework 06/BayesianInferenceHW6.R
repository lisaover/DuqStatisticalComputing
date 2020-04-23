#printPDF function receives two vectors and a filename. It prints three plots for each vector to the specified filename: ts plot, hist, acf.
 
printPDF <- function(m.s,s.s,filename) {
	
   pdf(filename)
    #PLOTS to PDF
  	par(mfrow=c(3,2)) #split plotting window into 3 rows and 2 columns
	ts.plot(m.s,xlab="Iterations") 
	ts.plot(s.s,xlab="Iterations") 
	hist(m.s,probability=T, cex.lab=1.5, cex.axis=1.5)
	hist(s.s,probability=T, cex.lab=1.5, cex.axis=1.5)
  	acf(m.s)
    acf(s.s) 
   dev.off()
	
}

#Gibbs function receives five parameters: "data" for the data values of interest, "m" for the population mean, "m0","s0","n0","k0" as prior parameters (mu, sigma sqr, nu, kappa), "N" for number of realizations, "lag" for determining how many realizations to skip between saves, and "burnin" for detemrining how many realizations to skip before starting to save. Gibbs generates N independent random normal values.

gibbs <- function(data,m,m0,s0,n0,k0,N,lag,burnin) {
	
	#obtain mean (y.bar), variance (s.sqr), and length (n) of data
	y.bar = mean(data)
    s.sqr = var(data)
    n = length(data)
	
	#Set N to be N*lag+burnin
	N <- N*lag + burnin

	#Initialize vectors to hold mu and sigma sqr
	m.s <- NULL
	s.s <- NULL
	#Vector to hold y values from normal distribution with each mu, sigma sqr as parameters
	y.s <- NULL
	
  for(i in 1:N) {
	
	#generate a sigma sqr s using mean (y.bar), variance (s.sqr), length (n) of data, and prior parameters
	s <- 1/rgamma(1,(n+1+n0)/2,((n-1)*s.sqr + n*(y.bar-m)^2 + ((m-m0)^2)*k0 + n0*s0)/2)
	#generate a mu m using mean (y.bar) and length (n) of data and sigma sqr (s) from last step
	m <- rnorm(1, (y.bar*n + m0*k0)/(n + k0), (sqrt((s/n)/(n + k0))))
	
	y.s <- c(y.s,rnorm(1,m,sqrt(s)))
	
  #if i is greater than burnin and if i is a multiple of the lag, store m and s
  if(i > burnin) {
  	if(i %% lag == 0) {
		m.s <- c(m.s,m)
		s.s <- c(s.s,s)
    }
   }
  }
  #filename = sprintf("Documents/R-FILES/HW6-%s-%s-%s.pdf",k0,n0,s0)
  #printPDF(m.s,s.s,filename)

  #vectors <- list("mu" = m.s, "sigsqr" = s.s)
  #return(vectors)
  
  return(y.s)
 }
 
data = read.table(file.choose(), header=F)
attach(data)
data = data$V1

#Set burnin to 10 and leave lag=1 and run Gibbs with N=2500
N = 25000
lag = 1
burnin = 10 
m = 80
#FIX m0 at 20
m0 = 20
#CHOOSE s0,n0,k0
k0 = 1
s0 = 1
n0 = 1

vector = gibbs(data,m,m0,s0,n0,k0,N,lag,burnin)
mean(vectors$mu)
quantile(vectors$mu, 0.025)
quantile(vectors$mu, 0.975)
mean(vectors$sigsqr)
quantile(vectors$sigsqr, 0.025)
quantile(vectors$sigsqr, 0.975)

y.10 = vector > 10
prob = length(y.10[y.10=="TRUE"])/N
#0.5758

quantile(vector, 0.025)
#0.9
quantile(vector, 0.975)
#21.34




###WRONG APPROACH BUT NICE CODE### 
#For each parameter, (k0, s0, n0), create two vectors of mean sigma sqr values and mean mu values. Start each set with parameter = 1. Continue running Gibbs 100 times each time making the parameter both larger (multiply by 5) and smaller (divide by 5).
large = 1
small = 1
small.m = NULL #stores mean mu values as parameter (k0, s0, n0) gets smaller
large.m = NULL #stores mean mu values as parameter (k0, s0, n0) gets larger
small.s = NULL #stores mean sigsqr values as (k0, s0, n0) parameter gets smaller
large.s = NULL #stores mean sigsqr values as (k0, s0, n0) parameter gets larger

# loop takes parameter, (k0, s0, n0), as small from #0^1 to #0^-100 and as large from #0^1 to #0^100
for(k in 1:1000) {
	#run Gibbs with a small parameter
	vectors = gibbs(data,m,m0,s0,n0,small,N,lag,burnin)
	#store the mean sigsqr and mean mu value in respective vectors 
	small.s = c(small.s,mean(vectors$sigsqr))
	small.m = c(small.m,mean(vectors$mu))
	
	#run Gibbs with a large parameter
	vectors = gibbs(data,m,m0,s0,n0,large,N,lag,burnin)
	large.s = c(large.s,mean(vectors$sigsqr))
	large.m = c(large.m,mean(vectors$mu))
	
	large = large+1
	small = small/2
}

quantile(large.s)
quantile(small.s)
quantile(large.m)
quantile(small.m)