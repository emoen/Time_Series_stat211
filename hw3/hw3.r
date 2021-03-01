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

#A: Plot the sample autocorrelation function ρ_hat sample partial autocorrelation function α_hat of the simulated process.
#How are these consistent with an AR(1) process? Can you see the relation between ø and ρ_hat(1)


set.seed(1234)
y <- arima.sim(
	model = list(ar = 0.7),
	n = 10000,
	innov = rnorm(10000, sd = 2)
)

par(mfrow=c(1,2)) #create a window with 1 row and 2 columns
AutoCorrelation = acf(x, plot = FALSE)
plot(AutoCorrelation, main = "Auto Corr. Func. AR(1)")
PartialAutoCorr = pacf(x, plot=FALSE)
plot(PartialAutoCorr, main="Partial Auto Corr. Func.")

#The ACF "tails off" to || x_n - x_m || < Epsilon after t>=10
#The PACF "cuts off" at t=2

'''
A\B  |      AR(P)     | MA(q)           | ARMA(p,q)
----------------------------------------------
ACF  |Tails off       | cuts off at t=q | Tails off
PACF |cuts off at t=p | Tails off       | Tails off
'''

alpha_1 =pacf(x, plot=FALSE)$acf[1]

#[1] 0.7011247
# alpha(1) approx ø

# B: Make a plot comparing the sample autocorrelation with the theoretical autocorrelation. Do the same for the partial autocorrelation

phi = .7
acf_theo = phi^(0:40) # Numerical Theoretical ACF with up to 40 lags - length 40
pacf_theo = c(phi, rep(0,40)) # Theoretical PACF with up to 40 lags - length 40
# Plotting ACF and PACF
par(mfrow=c(1,2))

# length acf(x, plot=FALSE)$acf = 41
plot(0:40, acf(x, plot=FALSE)$acf, type = "h", ylab = "ACF", xlab = "h")
lines(0:40, acf_theo , col=2, type = "l", lty=2)
legend("topright", col = 1:2,
legend = c(expression(hat(rho)(h)),expression(rho(h))), lty=1:2)
# length( pacf(x,plot=FALSE)$acf ) = 40
plot(1:40, pacf(x, plot=FALSE)$acf, type = "h", ylab = "PACF", xlab ="h")
lines(1:40, pacf_theo, col=2, type = "l",lty=2)
legend("topright", col = 1:2, legend = c(expression(hat(alpha)(h)),expression(alpha(h))), lty=1:2)

# length [acf = pacf +1]
