/*8.1*/
data one;
	mu = 50;
	sd = 20;
	seed = 46327;
	do i = 1 to 1000;
		z = rannor(seed);
		x = (z*sd) + mu;
		output;
	end;
run;

goptions cpattern = black htext = 1.5;

proc gchart;
vbar x/ space = 0 midpoints = -10 to 110 by 10 width = 8;
title 'Random Observations from a Normal Distribution with mu=50 and sigma=20';
run;

/*8.2*/
data one;
	mu = 10;
	sd = 10;
	N = 5000;
	seed = 46327;
	do i = 1 to N;
		z = rannor(seed);
		x = (z*sd) + mu;
		output;
	end;
run;

goptions cpattern = black htext = 1.5;
proc univariate;
var x;
histogram x/ normal endpoints = -20 to 40 by 5;
title1 'Distribution of the Mean of 5000 Samples';
title2 'from a N(10,10) Distribution';
run;

/*8.3*/
data one;
	lambda = 7;
	N = 1000;
	seed = 46327;
	do i = 1 to N;
		z = ranexp(seed);
		x = z/lambda;
		output;
	end;
run;

goptions cpattern = black htext = 1.5;

proc gchart;
vbar x/ space = 0 midpoints = 0 to 1 by .1 width = 5;
title 'Random Observations from an Exponential Distribution with lambda=7';
run;

/*8.4*/
data one;
	mu = 5;
	N = 700;
	seed = 46327;
	do i = 1 to N;
		x = ranpoi(seed,mu);
		output;
	end;
run;

goptions cpattern = black htext = 1.5;

proc gchart;
vbar x/ space = 0;
title 'Random Observations from a Poisson Distribution with mu=5';
run;

/*8.5*/
data one;
	N = 500;
	n = 40;
	p = 0.2;
	seed = 4491;
	do i to N;
		x = ranbin(seed, n, p);
		output;
	end;
	run;

	goptions cpattern = black htext = 1.5;

proc gchart;
vbar x/ space = 0;
title 'Random Observations from a Binomial Distribution with n=40 and p=0.2';
run;

/*8.6*/
data one;
	lambda = 7;
	N = 1000;
	seed = 46327;
	do i = 1 to N;
		z1 = ranexp(seed);
		x1 = z1/lambda;
		z2 = ranexp(seed);
		x2 = z2/lambda;
		z3 = ranexp(seed);
		x3 = z3/lambda;
		z4 = ranexp(seed);
		x4 = z4/lambda;
		z5 = ranexp(seed);
		x5 = z5/lambda;
		z6 = ranexp(seed);
		x6 = z6/lambda;
		z7 = ranexp(seed);
		x7 = z7/lambda;
		z8 = ranexp(seed);
		x8 = z8/lambda;
		z9 = ranexp(seed);
		x9 = z9/lambda;
		z10 = ranexp(seed);
		x10 = z10/lambda;
		avg_x = (x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)/10;
		output;
	end;
run;

goptions cpattern = black htext = 1.5;

proc univariate;
var avg_x;
histogram avg_x/ normal;
title1 'Distribution of the Mean of 1000 Samples';
title2 'from an Exponential(7) Distribution';
run;

/*8.7*/
data one;
	N = 1000;
	a = 10;
	b = 20;
	seed = 4491;
	do i to N;
		z1 = ranuni(seed);
		x1 = a+(b-a)*z1;
		z2 = ranuni(seed);
		x2 = a+(b-a)*z2;
		z1 = ranuni(seed);
		x1 = a+(b-a)*z1;
		z3 = ranuni(seed);
		x3 = a+(b-a)*z3;
		z4 = ranuni(seed);
		x4 = a+(b-a)*z4;
		z5 = ranuni(seed);
		x5 = a+(b-a)*z5;
		z6 = ranuni(seed);
		x6 = a+(b-a)*z6;
		z7 = ranuni(seed);
		x7 = a+(b-a)*z7;
		z8 = ranuni(seed);
		x8 = a+(b-a)*z8;
		z9 = ranuni(seed);
		x9 = a+(b-a)*z9;
		z10 = ranuni(seed);
		x10 = a+(b-a)*z10;
		avg_x = (x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)/10;
		output;
	end;
run;

goptions cpattern = black htext = 1.5;

proc univariate;
var avg_x;
histogram avg_x/ normal;
title1 'Distribution of the Mean of 1000 Samples';
title2 'from a Uniform(10,20) Distribution';
run;
