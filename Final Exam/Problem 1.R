#The proposal density is uniform over a 4 unit interval so the height of the proposal desnity is 0.25. The maximum height of the target density is 0.5 so multiply 0.25 by 2 (A=2) to raise the proposal density above the target density. 
A = 2
real.s = NULL
N = 25000

for(i in 1:N) {
	#Draw a uniform independent random variable and transform it to be within the domain of the given density, 2 <= x <= 6.
	x = runif(1) * 4 + 2
		
	#Determine what part of the piecewise function to use based on the value of x and set the numerator equal to the appropriate function and the denominator equal to A times the height of the uniform density over the four unit interval.
	if (x <= 3) {
		num = (x - 2)/2
		denom = A*0.25
		}
	else {
		num = (2 - (x/3))/2
		denom = A*0.25
		}
	#Calculate the probability that the previously generated uniform is from the target density
	prob <- num/denom
	
	#Generate a uniform independent random variable to use as a comparison probability: if the probability, prob, that the transformed x uniform could be from the given distribution is greater than the probability represented by this newly generated uniform random variable, select the x transformed uniform and store in real.s
	if(runif(1) < prob) {
		real.s <- c(real.s,x)
	}
}
hist(real.s)
mean(real.s)
#[1] 3.664185


