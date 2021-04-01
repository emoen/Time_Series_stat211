'''
For n = 5 and n = 1000, generate a random sample (e.g. using the R-command runif)
from the unif(0; 1) distribution. 
Perform the two tests in c) having size alpha = 0.05 on your datasets. 
For which of the tests do you reject H0? Discuss your fndings.
'''

min_n=5
max_n=1000
x_min = runif(min_n, min = 0, max = 1) #The Uniform Distribution
x_max = runif(max_n, min = 0, max = 1) #runif (open) will not generate either of the extreme values

#The likelihood ratio statistics rejects H0 if max(x) > 1 which it isnt from U(0,1)
# theta_MoM = 2*X_bar

theta_MoM_min = 2* mean(x_min)
theta_MoM_max = 2* mean(x_max)

sample.mean = 2*mean(x_min)
print(sample.mean)
sample.n = min_n
sample.sd = sd(x_min)
sample.se = sample.sd/sqrt(sample.n)
bound = ci_for_u(sample)

#2*sample_mean [1] 0.5812125
#t-score [1] 2.776445
#confidence interval [1] 0.2281778 0.9342472

sample.mean = 2* mean(x_max)
print(sample.mean)
sample.n = max_n
sample.sd = sd(x_max)
sample.se = sample.sd/sqrt(sample.n)
bound = ci_for_u(sample)
#2*sample_mean [1] 0.9991154
#t-score [1] 1.962341
#confidence interval [1] 0.9814277 1.0168030

# discuss: likelihood ratio test is better.
# theta_Mom contains level alpha=0.05 test with theta > 1 when we know U(0,1)

ci_for_u <- function(sample) {
	alpha = 0.05
	degrees.freedom = sample.n - 1
	t.score = qt(p=alpha/2, df=degrees.freedom,lower.tail=F)
	print(t.score)
	margin.error = t.score * sample.se
	bound.lower = sample.mean - margin.error
	bound.upper = sample.mean + margin.error
	print(c(bound.lower,bound.upper))
	return (c(bound.lower,bound.upper))
}


'''
 here is a little simulation for the MSE of the estimators of the mean,
 showing that while the MLE if we do not know the lower bound is zero is unbiased,
 the MSEs for the two variants are identical, suggesting that the estimator which 
 incorporates knowledge of the lower bound reduces variability.
 
https://stats.stackexchange.com/questions/252129/is-there-an-example-where-mle-produces-a-biased-estimate-of-the-mean 
'''

theta <- 1
mean <- theta/2
reps <- 500000
n <- 5
mse <- bias <- matrix(NA, nrow = reps, ncol = 2)

for (i in 1:reps){
  x <- runif(n, min = 0, max = theta)
  mle.knownlowerbound <- max(x)/2
  mle.unknownlowerbound <- (max(x)+min(x))/2
  mse[i,1] <- (mle.knownlowerbound-mean)^2
  mse[i,2] <- (mle.unknownlowerbound-mean)^2
  bias[i,1] <- mle.knownlowerbound-mean
  bias[i,2] <- mle.unknownlowerbound-mean

}

> colMeans(mse)
[1] 0.01194837 0.01194413

> colMeans(bias)
[1] -0.083464968 -0.000121968