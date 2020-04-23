#metro function receives nine parameters: two vectors "x.s" and “y.s” for the data values of interest, "a0" for the initial alpha, "b0" for the initial beta, "k" for the alpha interval, "c" for the beta interval, "N" for number of independent random normal realizations, "lag" for determining how many realizations to skip between saves, and "burnin" for determining how many realizations to skip before starting to save.

metro <- function(x.s,y.s,a0,b0,k,c,N,lag,burnin) {
	
	#Set N to be N*lag+burnin
	N <- N*lag + burnin

	#Initialize vectors to hold alpha (a.v) and beta (b.v) values
	a.s <- NULL
	b.s <- NULL
	
	#store the acceptance rate for alpha (a.cnt) and beta (b.cnt)
	a.cnt = 0
	b.cnt = 0
	
  for(i in 1:N) {
  	
#Generate an alpha star 'a.star' from the proposal density (Normal) using a0 as mu and k as std dev
	a.star = rnorm(1,a0,k)
	
#Compute probability for alpha using Normal priors with mean (0), variance (100) -- Use full conditional (from joint posterior which is likelihood*prior(alpha)*prior(beta)) -- compute ratio of full conditional given a.star to the full conditional given a0
	numerator = sum(y.s*a.star) - sum(log(1 + exp(a.star + b0*x.s))) - a.star^2/200
	denominator = sum(y.s*a0) - sum(log(1 + exp(a0 + b0*x.s))) - a0^2/200
	target.ratio = exp(numerator - denominator)
	
	#Use target.ratio to determine if a.star should be accepted
	param = 0
	accept.code = 0
	if(target.ratio < 1) {
		if(runif(1,0,1) < target.ratio) {
			param = a.star
			accept.code = 1
		}
		else {
			param = a0
			accept.code = 0
		}
	}
	else {
		param = a.star
		accept.code = 1
	}

#if i is greater than burnin and if i is a multiple of the lag, store value
#add accept.code to total a.cnt - accept.code will be 1 if a.star was accepted and 0 ow
#add accepted value to vector - a0 or a.star - stored as param
	if(i > burnin) {
  	if(i %% lag == 0) {
		a.cnt = a.cnt + accept.code
		a.s = c(a.s, param)
		#set a0 equal to param (parameter that was stored in vector -- a.star or a0) for next iteration
   		a0 = param
    }
   }
	
#Reset target.ratio, numerator, denominator
numerator = 0
denominator = 0
target.ratio = 0
#Generate a beta star 'b.star' from the proposal density (Normal) using b0 as mu and c as std dev
	b.star = rnorm(1,b0,c)
	
#Compute probability for beta using Normal priors with mean (0), variance (100) -- Use full conditional (from joint posterior which is likelihood*prior(alpha)*prior(beta)) -- compute ratio of full conditional given b.star to the full conditional given b0
	numerator = sum(y.s*b.star*x.s) - sum(log(1 + exp(a.star + b.star*x.s))) - b.star^2/200
	denominator = sum(y.s*b0*x.s) - sum(log(1 + exp(a.star + b0*x.s))) - b0^2/200
	target.ratio = exp(numerator - denominator)
	
	#Use target.ratio to determine if b.star should be accepted
	param = 0
	accept.code = 0
	if(target.ratio < 1) {
		if(runif(1,0,1) < target.ratio) {
			param = b.star
			accept.code = 1
		}
		else {
			param = b0
			accept.code = 0
		}
	}
	else {
		param = b.star
		accept.code = 1
	}
	
#if i is greater than burnin and if i is a multiple of the lag, store value
#add accept.code to total b.cnt - accept.code will be 1 if b.star was accepted and 0 ow
#add accepted value to vector - b0 or b.star - stored as param
	if(i > burnin) {
  	if(i %% lag == 0) {
		b.cnt = b.cnt + accept.code
		b.s = c(b.s, param)
		#set b0 equal to param (parameter that was stored in vector -- b.star or b0) for next iteration
   		b0 = param
    }
   }
   

  }

  vectors <- list("alpha" = a.s, "aCount" = a.cnt, "beta" = b.s, "bCount" = b.cnt)
  return(vectors)
 }
 
#Convert to dataframe and attach
df = data.frame(state.x77)
attach(df)

#Income per capita - create y.s vector or 1s and 0s where 1 is if Income is greater than income per capita and 0 otherwise
perCapita = 4445
y.s = as.numeric(Income > perCapita)
x.s = HS.Grad

#Obtain alpha and beta from data to use for a0 and b0
output = summary(glm(y.s~x.s,family=binomial))

#Set metro parameters
N = 10000
lag = 1000
burnin = 20
k=.35
c=.007
a0 = output$coef[1,1]
b0 = output$coef[2,1]

vectors = metro(x.s,y.s,a0,b0,k,c,N,lag,burnin)

aCount = vectors$aCount
bCount = vectors$bCount
alpha = vectors$alpha
beta = vectors$beta
 
par(mfrow=c(3,2)) #split plotting window into 3 rows and 2 columns
	ts.plot(alpha,xlab="Iterations") 
	ts.plot(beta,xlab="Iterations") 
	hist(alpha,probability=T, cex.lab=1.5, cex.axis=1.5)
	hist(beta,probability=T, cex.lab=1.5, cex.axis=1.5)
  	acf(alpha, lag.max=1000)
    acf(beta, lag.max=1000) 
    
alpha.accept = aCount/N
beta.accept = bCount/N
alpha.accept
#0.6963
beta.accept
#0.6824

mean(alpha)
mean(beta)
var(alpha)
var(beta)

quantile(vectors$beta, 0.025)
quantile(vectors$beta, 0.975)




###########################################

> mean(vectors$alpha)
[1] -7.959136
> mean(vectors$beta)
[1] 0.1575588
> var(vectors$alpha)
[1] 9.823369
> var(vectors$beta)
[1] 0.003514773


> quantile(vectors$alpha, 0.025)
     2.5% 
-14.40186 
> quantile(vectors$alpha, 0.975)
    97.5% 
-3.044593 
> quantile(vectors$beta, 0.025)
      2.5% 
0.06520356 
> quantile(vectors$beta, 0.975)
    97.5% 
0.2811901 

detach(df)

> output

Call:
glm(formula = y.s ~ x.s, family = binomial)

Deviance Residuals: 
    Min       1Q   Median       3Q      Max  
-2.3269  -0.6611   0.5067   0.9300   1.4732  

Coefficients:
            Estimate Std. Error z value Pr(>|z|)   
(Intercept) -8.02370    2.63976  -3.040  0.00237 **
x.s          0.15842    0.04976   3.184  0.00145 **
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

(Dispersion parameter for binomial family taken to be 1)

    Null deviance: 68.029  on 49  degrees of freedom
Residual deviance: 53.848  on 48  degrees of freedom
AIC: 57.848

Number of Fisher Scoring iterations: 4

#DATA TO SAS

 
#Convert to dataframe and attach
df = data.frame(state.x77)
attach(df)

#Income per capita - create y.s vector or 1s and 0s where 1 is if Income is greater than income per capita and 0 otherwise
perCapita = 4445
y.s = as.numeric(Income > perCapita)
x.s = HS.Grad

# Create dataframe of y.s and x.s for export
thedata = data.frame(x.s, y.s)

# export data frame to text file 
 # write out text datafile and
 # a SAS program to read it
 library(foreign)
 write.foreign(thedata, "c:/Users/TA00/My Documents/thedata.txt", "c:/Users/TA00/My Documents/thedata.sas",   package="SAS") 
 