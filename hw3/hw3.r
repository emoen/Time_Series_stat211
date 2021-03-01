'''
In this homework, you are supposed to learn how to simulate and estimate AR/MA models
in R, and how to use ACF and PACF to identify models
'''

# 3.1
'''
Consider an AR(1) process f(t) = ø*f(t-1) + z(t), z~N(0,4)   |ø| < 1
Simulate with ø = 0.7, t=1..n, n=10000=10^4

Note: ACF of AR(1): p = ø^h, h >= 1
partial autocorrelation function alpha(h) = 1 if h=1 else 0
'''

#A: 1) Plot the sample autocorrelation function ρ_hat sample partial autocorrelation function α_hat of the simulated process.
#2) how are these consistent with an AR(1) process? Can you see the relation between ø_a and ρ_hat(1)


set.seed(1234)
y <- arima.sim(
	model = list(ar = 0.7),
	n = 10000,
	innov = rnorm(10000, sd = 2)
)

par(mfrow=c(1,2)) #create a window with 1 row and 2 columns
acf(x)
AutoCorrelation = acf(x, plot = FALSE)
plot(AutoCorrelation, main = "Auto Correlation Function of AR(1)")
PartialAutoCorr = pacf(x, plot=FALSE)
plot(PartialAutoCorr, main="Partial Auto Correlation Function")